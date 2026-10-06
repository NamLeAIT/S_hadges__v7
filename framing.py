from __future__ import annotations

import bz2
import gzip
import io
import lzma
import mimetypes
import re
import struct
import zipfile
import zlib
from dataclasses import dataclass, asdict
from typing import Callable, Iterable

import storage_container

PAYLOAD_BITS = 103
HEDGES_FILE_MAGIC = b"HDF1"
HEDGES_FILE_VERSION = 1
HEDGES_FILE_HEADER = struct.Struct(">4sBBBII")
DOMAIN_CODES = {"binary": 0, "text": 1, "image": 2, "audio": 3, "video": 4, "document": 5, "archive": 6}
CODE_DOMAINS = {value: key for key, value in DOMAIN_CODES.items()}


@dataclass(frozen=True)
class ReconstructionCandidate:
    hypothesis: str
    bits: str
    data: bytes
    format: str
    extension: str
    mime: str
    strong_valid: bool
    validation: str
    padding_bits: int
    filename: str | None = None
    domain: str | None = None

    def to_dict(self) -> dict:
        d = asdict(self)
        d.pop("bits", None)
        d.pop("data", None)
        return d


def bytes_to_bits(data: bytes) -> str:
    return "".join(f"{b:08b}" for b in bytes(data or b""))


def bits_to_bytes_exact(bits: str) -> bytes:
    if any(c not in "01" for c in bits):
        raise ValueError("Bitstream must contain only 0/1")
    if len(bits) % 8:
        raise ValueError("Bitstream is not byte aligned")
    return bytes(int(bits[i:i+8], 2) for i in range(0, len(bits), 8))


def domain_from_mime(mime: str | None) -> str:
    value = str(mime or "").lower()
    if value.startswith("text/"): return "text"
    if value.startswith("image/"): return "image"
    if value.startswith("audio/"): return "audio"
    if value.startswith("video/"): return "video"
    if value in {"application/pdf", "application/msword"} or "officedocument" in value: return "document"
    if value in {"application/zip", "application/gzip", "application/x-xz", "application/x-bzip2"}: return "archive"
    return "binary"


def pack_hedges_file(data: bytes, filename: str, mime: str | None = None) -> bytes:
    """Embed file length, name/domain and CRC inside the HEDGES bitstream."""
    data = bytes(data or b"")
    name = str(filename or "stored.bin").rsplit("/", 1)[-1].rsplit("\\", 1)[-1] or "stored.bin"
    name_b = name.encode("utf-8")
    if len(name_b) > 255:
        raise ValueError("Stored file name is too long for HEDGES framing")
    if len(data) > 0xFFFFFFFF:
        raise ValueError("Stored file is too large for HEDGES framing")
    domain = domain_from_mime(mime or mimetypes.guess_type(name)[0])
    header = HEDGES_FILE_HEADER.pack(
        HEDGES_FILE_MAGIC,
        HEDGES_FILE_VERSION,
        DOMAIN_CODES[domain],
        len(name_b),
        len(data),
        zlib.crc32(data) & 0xFFFFFFFF,
    )
    return header + name_b + data


def hedges_file_overhead(filename: str) -> int:
    name = str(filename or "stored.bin").rsplit("/", 1)[-1].rsplit("\\", 1)[-1] or "stored.bin"
    return HEDGES_FILE_HEADER.size + len(name.encode("utf-8"))


def _hedges_envelope_candidate(full_bits: str) -> ReconstructionCandidate | None:
    """Recover a self-contained HEDGES file without external size/domain hints."""
    if len(full_bits) < HEDGES_FILE_HEADER.size * 8 or any(c not in "01" for c in full_bits):
        return None
    prefix = bits_to_bytes_exact(full_bits[:HEDGES_FILE_HEADER.size * 8])
    magic, version, domain_code, name_n, data_n, checksum = HEDGES_FILE_HEADER.unpack(prefix)
    if magic != HEDGES_FILE_MAGIC:
        return None
    if version != HEDGES_FILE_VERSION or domain_code not in CODE_DOMAINS:
        raise ValueError("HEDGES file envelope version/domain is invalid")
    envelope_n = HEDGES_FILE_HEADER.size + int(name_n) + int(data_n)
    used_bits = envelope_n * 8
    if used_bits > len(full_bits):
        raise ValueError("HEDGES file envelope is truncated")
    tail = full_bits[used_bits:]
    if tail and tail != "1" + "0" * (len(tail) - 1):
        raise ValueError("HEDGES padding delimiter is invalid")
    envelope = bits_to_bytes_exact(full_bits[:used_bits])
    name_b = envelope[HEDGES_FILE_HEADER.size:HEDGES_FILE_HEADER.size + int(name_n)]
    data = envelope[HEDGES_FILE_HEADER.size + int(name_n):]
    if (zlib.crc32(data) & 0xFFFFFFFF) != int(checksum):
        raise ValueError("HEDGES reconstructed file CRC32 mismatch")
    name = name_b.decode("utf-8", errors="replace") or "recovered.bin"
    fmt, detected_ext, detected_mime, _strong, note = detect_and_validate_exact(data)
    suffix = mimetypes.guess_type(name)[0]
    ext = ("." + name.rsplit(".", 1)[-1]) if "." in name else detected_ext
    mime = suffix or detected_mime or "application/octet-stream"
    domain = CODE_DOMAINS[int(domain_code)]
    return ReconstructionCandidate(
        "H2 · embedded length/name/CRC", bytes_to_bits(data), data, fmt, ext or ".bin", mime,
        True, f"Embedded length + CRC32 valid; {note}", len(tail), filename=name, domain=domain,
    )


def frame_103(bits: str) -> tuple[list[str], int, str]:
    """Self-delimiting tail framing without adding an extra payload block.

    p = 0: B
    p > 0: B || 1 || 0^(p-1)
    """
    if any(c not in "01" for c in bits):
        raise ValueError("bits must contain only 0/1")
    if not bits:
        return [], 0, ""
    p = (-len(bits)) % PAYLOAD_BITS
    framed = bits if p == 0 else bits + "1" + "0" * (p - 1)
    if len(framed) % PAYLOAD_BITS:
        raise AssertionError("103-bit framing failed")
    return [framed[i:i+PAYLOAD_BITS] for i in range(0, len(framed), PAYLOAD_BITS)], p, framed


def delimiter_candidate(full_bits: str) -> tuple[str, int] | None:
    """Return candidate before 1 0* delimiter and removed-bit count."""
    if not full_bits or any(c not in "01" for c in full_bits):
        return None
    i = len(full_bits) - 1
    while i >= 0 and full_bits[i] == "0":
        i -= 1
    if i < 0 or full_bits[i] != "1":
        return None
    candidate = full_bits[:i]
    removed = len(full_bits) - i
    return candidate, removed


def _png_end(data: bytes) -> tuple[int | None, str]:
    sig = b"\x89PNG\r\n\x1a\n"
    if not data.startswith(sig):
        return None, "not PNG"
    pos = 8
    seen_ihdr = False
    while pos + 12 <= len(data):
        n = int.from_bytes(data[pos:pos+4], "big")
        typ = data[pos+4:pos+8]
        end = pos + 12 + n
        if n > 100_000_000 or end > len(data):
            return None, "truncated PNG chunk"
        body = data[pos+8:pos+8+n]
        stored = int.from_bytes(data[pos+8+n:end], "big")
        crc = zlib.crc32(typ)
        crc = zlib.crc32(body, crc) & 0xFFFFFFFF
        if stored != crc:
            return None, "PNG CRC mismatch"
        if typ == b"IHDR":
            if pos != 8 or seen_ihdr:
                return None, "invalid PNG IHDR"
            seen_ihdr = True
        if typ == b"IEND":
            if n != 0 or not seen_ihdr:
                return None, "invalid PNG IEND"
            return end, "PNG chunks + CRC + IEND valid"
        pos = end
    return None, "PNG IEND not found"


def _jpeg_end(data: bytes) -> tuple[int | None, str]:
    if len(data) < 4 or not data.startswith(b"\xff\xd8"):
        return None, "not JPEG"
    # A full JPEG parser is unnecessary for endpoint framing, but require a terminal EOI.
    end = data.rfind(b"\xff\xd9")
    if end < 2:
        return None, "JPEG EOI not found"
    return end + 2, "JPEG SOI/EOI valid"


def _riff_end(data: bytes) -> tuple[int | None, str, str, str]:
    if len(data) < 12 or data[:4] != b"RIFF":
        return None, "not RIFF", "", ""
    end = int.from_bytes(data[4:8], "little") + 8
    if end < 12 or end > len(data):
        return None, "invalid RIFF size", "", ""
    form = data[8:12]
    if form == b"WEBP":
        return end, "RIFF WEBP size valid", ".webp", "image/webp"
    if form == b"WAVE":
        return end, "RIFF WAVE size valid", ".wav", "audio/wav"
    if form == b"AVI ":
        return end, "RIFF AVI size valid", ".avi", "video/x-msvideo"
    return end, f"RIFF {form!r} size valid", ".riff", "application/octet-stream"


def _bmp_end(data: bytes) -> tuple[int | None, str]:
    if len(data) < 14 or data[:2] != b"BM":
        return None, "not BMP"
    end = int.from_bytes(data[2:6], "little")
    if end < 14 or end > len(data):
        return None, "invalid BMP file size"
    return end, "BMP header size valid"


def _gif_end(data: bytes) -> tuple[int | None, str]:
    if not (data.startswith(b"GIF87a") or data.startswith(b"GIF89a")):
        return None, "not GIF"
    # Trailer byte 0x3B terminates the GIF data stream.
    end = data.rfind(b"\x3b")
    if end < 6:
        return None, "GIF trailer not found"
    return end + 1, "GIF header/trailer valid"


def _zip_end(data: bytes) -> tuple[int | None, str]:
    if not data.startswith((b"PK\x03\x04", b"PK\x05\x06", b"PK\x07\x08")):
        return None, "not ZIP"
    # EOCD: signature + 18-byte fixed fields + 2-byte comment length => 22 bytes.
    start = max(0, len(data) - 65557 - 256)
    pos = data.rfind(b"PK\x05\x06", start)
    if pos < 0 or pos + 22 > len(data):
        return None, "ZIP EOCD not found"
    comment_len = int.from_bytes(data[pos+20:pos+22], "little")
    end = pos + 22 + comment_len
    if end > len(data):
        return None, "truncated ZIP EOCD/comment"
    try:
        with zipfile.ZipFile(io.BytesIO(data[:end]), "r") as zf:
            bad = zf.testzip()
            if bad is not None:
                return None, f"ZIP CRC failed at {bad}"
    except Exception as exc:
        return None, f"ZIP validation failed: {exc}"
    return end, "ZIP EOCD + CRC validation passed"


def _pdf_end(data: bytes) -> tuple[int | None, str]:
    if not data.startswith(b"%PDF-"):
        return None, "not PDF"
    pos = data.rfind(b"%%EOF")
    if pos < 0:
        return None, "PDF EOF not found"
    end = pos + 5
    while end < len(data) and data[end] in b"\x00\x09\x0a\x0c\x0d\x20":
        end += 1
    return end, "PDF header/EOF valid"


def _mp4_end(data: bytes) -> tuple[int | None, str]:
    if len(data) < 8:
        return None, "too short for ISO-BMFF"
    pos = 0
    seen_ftyp = False
    boxes = 0
    while pos + 8 <= len(data):
        size = int.from_bytes(data[pos:pos+4], "big")
        typ = data[pos+4:pos+8]
        header = 8
        if size == 1:
            if pos + 16 > len(data):
                return None, "truncated extended MP4 box"
            size = int.from_bytes(data[pos+8:pos+16], "big")
            header = 16
        elif size == 0:
            size = len(data) - pos
        if size < header or pos + size > len(data):
            return None, "invalid MP4 box length"
        if typ == b"ftyp":
            seen_ftyp = True
        pos += size
        boxes += 1
        if pos == len(data):
            break
    if seen_ftyp and boxes > 0:
        return pos, "ISO-BMFF boxes + ftyp valid"
    return None, "MP4 ftyp not found"


def _zlib_end(data: bytes) -> tuple[int | None, str]:
    """Validate one complete RFC 1950 zlib stream and recover its exact endpoint."""
    if len(data) < 6:
        return None, "too short for zlib"
    cmf, flg = data[0], data[1]
    if (cmf & 0x0F) != 8 or ((cmf << 8) + flg) % 31 != 0:
        return None, "not zlib"
    try:
        obj = zlib.decompressobj(zlib.MAX_WBITS)
        obj.decompress(data)
        obj.flush()
        if not obj.eof:
            return None, "incomplete zlib stream"
        end = len(data) - len(obj.unused_data)
        return end, "zlib stream + Adler-32 validated"
    except Exception as exc:
        return None, f"zlib validation failed: {exc}"


def _gzip_end(data: bytes) -> tuple[int | None, str]:
    if not data.startswith(b"\x1f\x8b"):
        return None, "not gzip"
    try:
        obj = zlib.decompressobj(16 + zlib.MAX_WBITS)
        obj.decompress(data)
        if not obj.eof:
            return None, "incomplete gzip stream"
        end = len(data) - len(obj.unused_data)
        return end, "gzip stream validated"
    except Exception as exc:
        return None, f"gzip validation failed: {exc}"


def _bz2_end(data: bytes) -> tuple[int | None, str]:
    if not data.startswith(b"BZh"):
        return None, "not bzip2"
    try:
        d = bz2.BZ2Decompressor()
        d.decompress(data)
        if not d.eof:
            return None, "incomplete bzip2 stream"
        end = len(data) - len(d.unused_data)
        return end, "bzip2 stream validated"
    except Exception as exc:
        return None, f"bzip2 validation failed: {exc}"


def _xz_end(data: bytes) -> tuple[int | None, str]:
    if not data.startswith(b"\xfd7zXZ\x00"):
        return None, "not xz"
    try:
        d = lzma.LZMADecompressor(format=lzma.FORMAT_XZ)
        d.decompress(data)
        if not d.eof:
            return None, "incomplete xz stream"
        end = len(data) - len(d.unused_data)
        return end, "XZ stream validated"
    except Exception as exc:
        return None, f"XZ validation failed: {exc}"


def _ddsc_end(data: bytes) -> tuple[int | None, str]:
    if not data.startswith(storage_container.MAGIC):
        return None, "not DDSC"
    try:
        restored = storage_container.unpack(data)
        return len(data), f"DDSC length + CRC32 valid ({restored['algorithm']})"
    except Exception as exc:
        return None, f"DDSC validation failed: {exc}"


def _flac_end(data: bytes) -> tuple[int | None, str]:
    if not data.startswith(b"fLaC"):
        return None, "not FLAC"
    # FLAC has no mandatory byte-count footer. Signature alone is not a strong endpoint proof.
    return len(data), "FLAC signature valid; endpoint weak"


def detect_and_validate_exact(data: bytes) -> tuple[str, str, str, bool, str]:
    """Return format, extension, mime, strong_valid, validation.

    Strong means the format itself establishes an endpoint equal to len(data), not
    merely that a recognizable prefix exists.
    """
    checks = []
    end, note = _ddsc_end(data); checks.append((end, "DDSC", ".ddsc", "application/vnd.dna-storage.compressed", note))
    end, note = _png_end(data); checks.append((end, "PNG", ".png", "image/png", note))
    end, note = _jpeg_end(data); checks.append((end, "JPEG", ".jpg", "image/jpeg", note))
    r_end, r_note, r_ext, r_mime = _riff_end(data)
    if r_end is not None:
        r_fmt = {".webp":"WebP", ".wav":"WAV", ".avi":"AVI"}.get(r_ext, "RIFF")
        checks.append((r_end, r_fmt, r_ext, r_mime, r_note))
    end, note = _bmp_end(data); checks.append((end, "BMP", ".bmp", "image/bmp", note))
    end, note = _gif_end(data); checks.append((end, "GIF", ".gif", "image/gif", note))
    end, note = _zip_end(data); checks.append((end, "ZIP", ".zip", "application/zip", note))
    end, note = _pdf_end(data); checks.append((end, "PDF", ".pdf", "application/pdf", note))
    end, note = _mp4_end(data); checks.append((end, "MP4/ISO-BMFF", ".mp4", "video/mp4", note))
    end, note = _zlib_end(data); checks.append((end, "ZLIB", ".zlib", "application/zlib", note))
    end, note = _gzip_end(data); checks.append((end, "GZIP", ".gz", "application/gzip", note))
    end, note = _bz2_end(data); checks.append((end, "BZIP2", ".bz2", "application/x-bzip2", note))
    end, note = _xz_end(data); checks.append((end, "XZ", ".xz", "application/x-xz", note))

    for end, fmt, ext, mime, note in checks:
        if end is not None and end == len(data):
            return fmt, ext, mime, True, note

    # Text is deliberately weak: UTF-8 validity does not prove an endpoint.
    if data:
        try:
            text = data.decode("utf-8")
            printable = sum(ch.isprintable() or ch in "\r\n\t" for ch in text) / max(1, len(text))
            if printable >= 0.98:
                return "UTF-8 text", ".txt", "text/plain", False, "UTF-8 text candidate; endpoint is not intrinsically self-describing"
        except Exception:
            pass
    return "Unknown binary", ".bin", "application/octet-stream", False, "No strong self-describing container endpoint detected"


def restore_lossless_payload(data: bytes) -> dict | None:
    """Deterministically expose the payload inside a recovered lossless wrapper.

    This function uses only the recovered bytes themselves. It never accepts the
    original file, compression method, filename, SHA, payload bits, or session
    metadata. The exact recovered wrapper is preserved by the caller; this helper
    only returns an additional inner-payload view when unwrapping is deterministic.
    """
    raw = bytes(data or b"")
    if not raw:
        return None

    fmt, ext, mime, strong, note = detect_and_validate_exact(raw)
    if not strong:
        return None

    inner = None
    inner_name = None
    method = None
    try:
        if fmt == "DDSC":
            restored = storage_container.unpack(raw)
            inner = restored["data"]
            inner_name = restored["name"]
            method = restored["algorithm"]
        elif fmt == "ZLIB":
            inner = zlib.decompress(raw)
            method = "zlib"
        elif fmt == "GZIP":
            inner = gzip.decompress(raw)
            method = "gzip"
        elif fmt == "BZIP2":
            inner = bz2.decompress(raw)
            method = "bzip2"
        elif fmt == "XZ":
            inner = lzma.decompress(raw, format=lzma.FORMAT_XZ)
            method = "xz"
        elif fmt == "ZIP":
            with zipfile.ZipFile(io.BytesIO(raw), "r") as zf:
                members = [x for x in zf.infolist() if not x.is_dir()]
                if len(members) != 1:
                    return None
                member = members[0]
                # CRC is checked by ZipFile.read(). One-member archives are
                # deterministic to unwrap; the recovered ZIP bytes remain available.
                inner = zf.read(member)
                inner_name = member.filename.rsplit("/", 1)[-1] or None
                method = "zip single-member"
        else:
            return None
    except Exception:
        return None

    if inner is None:
        return None

    inner_fmt, inner_ext, inner_mime, inner_strong, inner_note = detect_and_validate_exact(inner)
    if inner_name:
        guessed_mime = None
        try:
            import mimetypes as _mimetypes
            guessed_mime = _mimetypes.guess_type(inner_name)[0]
        except Exception:
            guessed_mime = None
        if guessed_mime:
            inner_mime = guessed_mime
        suffix = "." + inner_name.rsplit(".", 1)[-1] if "." in inner_name else ""
        if suffix:
            inner_ext = suffix
    else:
        inner_name = "recovered" + (inner_ext or ".bin")

    return {
        "data": inner,
        "name": inner_name,
        "format": inner_fmt,
        "extension": inner_ext or ".bin",
        "mime": inner_mime or "application/octet-stream",
        "strong_valid": bool(inner_strong),
        "validation": inner_note,
        "wrapper_format": fmt,
        "wrapper_validation": note,
        "method": method,
    }


def reconstruction_hypotheses(full_bits: str) -> list[tuple[str, str, int]]:
    if not full_bits or any(c not in "01" for c in full_bits):
        return []
    out: list[tuple[str, str, int]] = []
    # H0 only if the complete recovered stream can itself be a byte file.
    if len(full_bits) % 8 == 0:
        out.append(("H0 · no padding", full_bits, 0))
    d = delimiter_candidate(full_bits)
    if d is not None:
        bits, removed = d
        if len(bits) % 8 == 0:
            out.append(("H1 · self-delimiting padding", bits, removed))
    # De-duplicate identical bit candidates while retaining the first label.
    seen = set(); unique = []
    for row in out:
        if row[1] not in seen:
            seen.add(row[1]); unique.append(row)
    return unique


def autonomous_reconstruct(full_bits: str) -> tuple[str, ReconstructionCandidate | None, list[ReconstructionCandidate]]:
    """Fail-closed autonomous file reconstruction from decoded payload bits.

    No file size, MIME, filename, reference bytes, or padding count is accepted.
    """
    # New files carry a self-contained envelope.  Its byte length and CRC remove
    # padding without using uploaded-file/session metadata.  Legacy streams fall
    # through to the original format-based hypotheses below.
    if full_bits.startswith("".join(f"{b:08b}" for b in HEDGES_FILE_MAGIC)):
        try:
            embedded = _hedges_envelope_candidate(full_bits)
        except Exception:
            return "FAIL", None, []
        if embedded is not None:
            return "OK", embedded, [embedded]

    candidates: list[ReconstructionCandidate] = []
    for hyp, bits, p in reconstruction_hypotheses(full_bits):
        try:
            data = bits_to_bytes_exact(bits)
        except Exception:
            continue
        fmt, ext, mime, strong, note = detect_and_validate_exact(data)
        candidates.append(ReconstructionCandidate(hyp, bits, data, fmt, ext, mime, strong, note, p))

    strong = [c for c in candidates if c.strong_valid]
    if len(strong) == 1:
        return "OK", strong[0], candidates
    if len(strong) > 1:
        return "AMBIGUOUS", None, candidates

    # If there is exactly one byte-aligned hypothesis, framing itself is unique.
    # This safely covers unknown binary for p>0 when H0 is byte-impossible.
    if len(candidates) == 1:
        return "OK_WEAK", candidates[0], candidates
    if len(candidates) == 0:
        return "FAIL", None, candidates
    return "AMBIGUOUS", None, candidates
