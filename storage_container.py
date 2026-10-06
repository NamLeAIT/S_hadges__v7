"""Self-describing lossless-compression container used before DNA encoding.

The compressed stream carries its algorithm, original file name, MIME type,
uncompressed byte length and CRC32.  Nothing from Streamlit session state is
needed to restore the uploaded file after DNA decoding.
"""
from __future__ import annotations

import bz2
import lzma
import mimetypes
import struct
import zlib
from pathlib import Path

MAGIC = b"DDSC"
VERSION = 1
HEADER = struct.Struct(">4sBBHHQQI")

ALGORITHM_IDS = {
    "zstandard": 1,
    "brotli": 2,
    "xz": 3,
    "zlib": 4,
    "bzip2": 5,
}
ID_ALGORITHMS = {value: key for key, value in ALGORITHM_IDS.items()}


def _compress(algorithm: str, raw: bytes, level: int | None = None) -> bytes:
    algorithm = str(algorithm).lower()
    if algorithm == "zstandard":
        import zstandard as zstd
        return zstd.ZstdCompressor(level=int(level if level is not None else 19)).compress(raw)
    if algorithm == "brotli":
        import brotli
        return brotli.compress(raw, quality=int(level if level is not None else 11), mode=brotli.MODE_TEXT)
    if algorithm == "xz":
        return lzma.compress(raw, format=lzma.FORMAT_XZ, preset=int(level if level is not None else 9))
    if algorithm == "zlib":
        return zlib.compress(raw, level=int(level if level is not None else 9))
    if algorithm == "bzip2":
        return bz2.compress(raw, compresslevel=int(level if level is not None else 9))
    raise ValueError(f"Unsupported storage compressor: {algorithm}")


def _decompress(algorithm: str, payload: bytes) -> bytes:
    if algorithm == "zstandard":
        import zstandard as zstd
        return zstd.ZstdDecompressor().decompress(payload)
    if algorithm == "brotli":
        import brotli
        return brotli.decompress(payload)
    if algorithm == "xz":
        return lzma.decompress(payload, format=lzma.FORMAT_XZ)
    if algorithm == "zlib":
        return zlib.decompress(payload)
    if algorithm == "bzip2":
        return bz2.decompress(payload)
    raise ValueError(f"Unsupported storage compressor: {algorithm}")


def pack(raw: bytes, algorithm: str, name: str, mime: str | None = None, level: int | None = None) -> bytes:
    raw = bytes(raw or b"")
    algorithm = str(algorithm).lower()
    if algorithm not in ALGORITHM_IDS:
        raise ValueError(f"Unsupported storage compressor: {algorithm}")
    clean_name = Path(str(name or "payload.bin")).name or "payload.bin"
    clean_mime = str(mime or mimetypes.guess_type(clean_name)[0] or "application/octet-stream")
    name_b = clean_name.encode("utf-8")
    mime_b = clean_mime.encode("ascii", errors="replace")
    if len(name_b) > 65535 or len(mime_b) > 65535:
        raise ValueError("File name or MIME type is too long for the storage container")
    compressed = _compress(algorithm, raw, level)
    header = HEADER.pack(
        MAGIC,
        VERSION,
        ALGORITHM_IDS[algorithm],
        len(name_b),
        len(mime_b),
        len(raw),
        len(compressed),
        zlib.crc32(raw) & 0xFFFFFFFF,
    )
    return header + name_b + mime_b + compressed


def unpack(blob: bytes) -> dict:
    blob = bytes(blob or b"")
    if len(blob) < HEADER.size:
        raise ValueError("Storage container is truncated")
    magic, version, algorithm_id, name_n, mime_n, raw_n, payload_n, checksum = HEADER.unpack_from(blob)
    if magic != MAGIC or version != VERSION:
        raise ValueError("Storage container magic/version mismatch")
    algorithm = ID_ALGORITHMS.get(int(algorithm_id))
    if algorithm is None:
        raise ValueError("Unknown storage compression algorithm")
    metadata_end = HEADER.size + int(name_n) + int(mime_n)
    expected_end = metadata_end + int(payload_n)
    if expected_end != len(blob):
        raise ValueError("Storage container length mismatch")
    name_b = blob[HEADER.size:HEADER.size + int(name_n)]
    mime_b = blob[HEADER.size + int(name_n):metadata_end]
    payload = blob[metadata_end:expected_end]
    raw = _decompress(algorithm, payload)
    if len(raw) != int(raw_n):
        raise ValueError("Restored file length mismatch")
    if (zlib.crc32(raw) & 0xFFFFFFFF) != int(checksum):
        raise ValueError("Restored file CRC32 mismatch")
    name = name_b.decode("utf-8", errors="replace") or "recovered.bin"
    mime = mime_b.decode("ascii", errors="replace") or mimetypes.guess_type(name)[0] or "application/octet-stream"
    return {
        "data": raw,
        "name": Path(name).name or "recovered.bin",
        "mime": mime,
        "algorithm": algorithm,
        "compressed_bytes": len(payload),
        "container_bytes": len(blob),
        "crc32": f"{int(checksum):08x}",
    }


def is_container(blob: bytes) -> bool:
    return bytes(blob or b"").startswith(MAGIC)
