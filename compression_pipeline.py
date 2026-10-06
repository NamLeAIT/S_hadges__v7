# compression_pipeline.py
# Ported from accepted App1_v3(app2) compression pipeline.
# Headerless self-describing compression benchmark
# ------------------------------------------------
# New policy:
#   - No Lossy/Lossless UI dependency.
#   - Compression mode benchmarks all available valid compressors.
#   - Only self-describing outputs are accepted, based on detect_magic().
#   - The smallest byte stream is selected for DNA encoding.
#
# Requirements from existing project:
#   - utils_bits_v2.py: detect_magic, safe_basename
#   - ui_helpers.py: fmt_bytes, get_domain
#   - config.py: SELF_DESCRIBING_KINDS
#
# Optional:
#   - Pillow for image codecs
#   - ffmpeg on PATH for audio/video codecs

from __future__ import annotations

import gzip
import hashlib
import io
import bz2
import lzma
import mimetypes
import os
import subprocess
import tempfile
import zipfile
import zlib
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple

try:
    from PIL import Image
except Exception:
    Image = None

from utils_bits_v2 import detect_magic, safe_basename
from config import SELF_DESCRIBING_KINDS
from ui_helpers import fmt_bytes, get_domain
import storage_container


# ============================================================
# Data model
# ============================================================

@dataclass
class CompressionCandidate:
    rank: int
    method: str
    data: bytes
    ext: str
    kind: str
    mime: str
    lossy: bool
    size_bytes: int
    compression_ratio: float
    saving_pct: float
    estimated_dna_nt: int
    note: str = ""
    codec_family: str = ""
    quality_tier: str = ""

    def public_row(self) -> dict:
        return {
            "Rank": self.rank,
            "Method": self.method,
            
            "File extension": self.ext,
            
            "Size": fmt_bytes(self.size_bytes),
            "Ratio": f"{self.compression_ratio:.2f}",
           
            "Estimate DNA (nt)": f"{self.estimated_dna_nt:,}"
        }


def _candidate_from_bytes(
    method: str,
    data: bytes,
    original_size: int,
    lossy: bool,
    note: str = "",
    codec_family: str = "",
    quality_tier: str = "",
) -> Optional[CompressionCandidate]:
    """Create a candidate only if bytes are self-describing by magic signature."""
    if not data:
        return None

    m = detect_magic(data)
    if not m:
        return None

    if m.kind not in SELF_DESCRIBING_KINDS:
        return None

    size = len(data)
    ratio = (original_size / size) if size else 0.0
    saving = (1.0 - (size / original_size)) * 100.0 if original_size else 0.0

    return CompressionCandidate(
        rank=0,
        method=method,
        data=data,
        ext=m.ext,
        kind=m.kind,
        mime=m.mime,
        lossy=lossy,
        size_bytes=size,
        compression_ratio=ratio,
        saving_pct=saving,
        estimated_dna_nt=size * 4,  # exact for simple 2-bit mapping; estimate for rule-based mapping
        note=note or getattr(m, "note", ""),
        codec_family=codec_family or method.split("_", 1)[0],
        quality_tier=quality_tier,
    )


def _lossless_container_candidate(
    raw: bytes,
    original_name: str,
    original_size: int,
    algorithm: str,
    level: int,
    method: str,
    note: str,
) -> Optional[CompressionCandidate]:
    """Create a self-contained compressor result that restores name + MIME."""
    try:
        packed = storage_container.pack(
            raw,
            algorithm,
            original_name,
            mimetypes.guess_type(original_name)[0] or "application/octet-stream",
            level,
        )
    except Exception:
        return None
    return _candidate_from_bytes(
        method, packed, original_size, lossy=False, note=note,
        codec_family=algorithm, quality_tier="Exact",
    )


# ============================================================
# Generic self-describing containers
# ============================================================

def zip_store_bytes(name: str, data: bytes) -> bytes:
    """ZIP container with no compression; useful to force self-describing bytes."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", compression=zipfile.ZIP_STORED) as zf:
        zf.writestr(safe_basename(name), data)
    return buf.getvalue()


def zip_deflate_bytes(name: str, data: bytes, level: int = 6) -> bytes:
    """ZIP container with DEFLATE compression."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=int(level)) as zf:
        zf.writestr(safe_basename(name), data)
    return buf.getvalue()


def generic_container_candidates(raw: bytes, original_name: str, original_size: int) -> List[CompressionCandidate]:
    """
    Generic lossless self-describing candidates.

    Used for:
      - text/json/csv
      - binary/unknown
      - fallback for image/audio/video/document/archive
    """
    out: List[CompressionCandidate] = []

    for algorithm, level, method in (
        ("zstandard", 19, "zstd_level19"),
        ("xz", 9, "xz_preset9"),
        ("zlib", 9, "zlib_level9"),
    ):
        c = _lossless_container_candidate(
            raw, original_name, original_size, algorithm, level, method,
            "general-purpose exact compression",
        )
        if c:
            out.append(c)

    return out


def text_compression_candidate(raw: bytes, original_name: str, original_size: int) -> List[CompressionCandidate]:
    """High-ratio exact candidates for text and structured text."""
    out: List[CompressionCandidate] = []
    for algorithm, level, method in (
        ("brotli", 11, "brotli_quality11"),
        ("zstandard", 19, "zstd_level19"),
        ("xz", 9, "xz_preset9"),
    ):
        c = _lossless_container_candidate(
            raw, original_name, original_size, algorithm, level, method,
            "text/structured-data exact compression",
        )
        if c:
            out.append(c)
    return out


def document_compression_candidates(raw: bytes, original_name: str, original_size: int) -> List[CompressionCandidate]:
    """Lossless candidates for PDF and Office-like documents."""
    return generic_container_candidates(raw, original_name, original_size)


def fast_generic_container_candidates(raw: bytes, original_name: str, original_size: int) -> List[CompressionCandidate]:
    """Small fallback set for media files, where ffmpeg transcodes dominate runtime."""
    return generic_container_candidates(raw, original_name, original_size)


# ============================================================
# Image candidates
# ============================================================

def image_candidates(raw: bytes, original_name: str, original_size: int, quality_mode: str = "Compression") -> List[CompressionCandidate]:
    """
    Image-native candidates.

    In the new pipeline, Compression mode always tries both:
      - PNG and WebP lossless representatives
      - JPEG/WebP lossy: light, medium and strong compression

    quality_mode is kept only for backward compatibility and is not used.
    """
    out: List[CompressionCandidate] = []
    if Image is None:
        return out

    try:
        img = Image.open(io.BytesIO(raw))
        img.load()
    except Exception:
        return out
    # A single-frame preview metric cannot represent an animated GIF/WebP or a
    # multipage TIFF.  Keep those files on the general-purpose exact path.
    if int(getattr(img, "n_frames", 1)) != 1:
        return out

    png_img = img
    if png_img.mode not in ("RGB", "RGBA", "L", "P"):
        png_img = png_img.convert("RGBA")
    try:
        buf = io.BytesIO()
        png_img.save(buf, format="PNG", optimize=True, compress_level=9)
        c = _candidate_from_bytes(
            "png_optimized", buf.getvalue(), original_size, lossy=False,
            note="image-native exact-pixel compression", codec_family="PNG", quality_tier="Exact",
        )
        if c:
            out.append(c)
    except Exception:
        pass

    try:
        buf = io.BytesIO()
        img.save(buf, format="WEBP", lossless=True, method=6)
        c = _candidate_from_bytes(
            "webp_lossless", buf.getvalue(), original_size, lossy=False,
            note="image-native exact-pixel compression", codec_family="WebP", quality_tier="Exact",
        )
        if c:
            out.append(c)
    except Exception:
        pass

    jpeg_img = img if img.mode == "RGB" else img.convert("RGB")
    has_transparency = img.mode in ("RGBA", "LA") or (img.mode == "P" and "transparency" in img.info)
    for q, tier in ((90, "Light"), (70, "Medium"), (10, "Strong")):
        if has_transparency:
            break
        try:
            buf = io.BytesIO()
            jpeg_img.save(buf, format="JPEG", quality=int(q), optimize=True)
            c = _candidate_from_bytes(
                f"jpeg_{tier.lower()}_q{q}", buf.getvalue(), original_size, lossy=True,
                note="image JPEG lossy", codec_family="JPEG", quality_tier=tier,
            )
            if c:
                out.append(c)
        except Exception:
            pass

    for q, tier in ((90, "Light"), (70, "Medium"), (10, "Strong")):
        try:
            buf = io.BytesIO()
            img.save(buf, format="WEBP", quality=int(q), lossless=False)
            c = _candidate_from_bytes(
                f"webp_{tier.lower()}_q{q}", buf.getvalue(), original_size, lossy=True,
                note="image WebP lossy", codec_family="WebP", quality_tier=tier,
            )
            if c:
                out.append(c)
        except Exception:
            pass

    return out


# ============================================================
# ffmpeg helpers
# ============================================================

def has_ffmpeg() -> bool:
    """Return True if ffmpeg is available on PATH."""
    try:
        p = subprocess.run(
            ["ffmpeg", "-version"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=5,
        )
        return p.returncode == 0
    except Exception:
        return False


def run_ffmpeg_transcode(input_path: str, output_ext: str, args: List[str], timeout_sec: int = 120) -> Optional[bytes]:
    """
    Run ffmpeg and return output bytes.
    Returns None if ffmpeg is unavailable, codec fails, or timeout occurs.
    """
    try:
        with tempfile.TemporaryDirectory() as td:
            out_path = Path(td) / f"out{output_ext}"
            cmd = [
                "ffmpeg",
                "-y",
                "-hide_banner",
                "-loglevel",
                "error",
                "-i",
                input_path,
            ] + args + [str(out_path)]

            p = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=int(timeout_sec),
            )

            if p.returncode != 0 or not out_path.exists():
                return None

            return out_path.read_bytes()
    except Exception:
        return None


# ============================================================
# Audio candidates
# ============================================================

def audio_compression_candidates(
    input_path: str,
    raw: bytes,
    original_size: int,
    quality_mode: str = "Compression",
) -> List[CompressionCandidate]:
    """
    Audio-native candidates.

    New pipeline: always try audio-native candidates in Compression mode.

    Requires ffmpeg.  Keep the default set short so the Streamlit app remains
    responsive while still showing representative lossless and lossy choices.
    """
    out: List[CompressionCandidate] = []

    if not has_ffmpeg():
        return out

    # FLAC lossless
    flac = run_ffmpeg_transcode(input_path, ".flac", ["-vn", "-c:a", "flac"])
    if flac:
        c = _candidate_from_bytes(
            "flac_lossless", flac, original_size, lossy=False,
            note="audio-native exact-sample compression", codec_family="FLAC", quality_tier="Exact",
        )
        if c:
            out.append(c)

    # OGG/Opus lossy levels.
    for br, tier in ((128, "Light"), (64, "Medium"), (32, "Strong")):
        ogg = run_ffmpeg_transcode(
            input_path,
            ".ogg",
            ["-vn", "-c:a", "libopus", "-b:a", f"{br}k"],
        )
        if ogg:
            c = _candidate_from_bytes(
                f"opus_{tier.lower()}_{br}k", ogg, original_size, lossy=True,
                note="audio Opus lossy", codec_family="Opus", quality_tier=tier,
            )
            if c:
                out.append(c)

    return out


# ============================================================
# Video candidates
# ============================================================

def video_compression_candidates(
    input_path: str,
    raw: bytes,
    original_size: int,
    quality_mode: str = "Compression",
) -> List[CompressionCandidate]:
    """
    Video-native candidates.

    New pipeline: always try video-native candidates in Compression mode.

    Requires ffmpeg.  VP9 is intentionally omitted from the default benchmark
    because it is much slower than H.264 in an interactive app.
    """
    out: List[CompressionCandidate] = []

    if not has_ffmpeg():
        return out

    # H.264 MP4 levels. Higher CRF means stronger compression.
    for crf, tier in ((20, "Light"), (28, "Medium"), (36, "Strong")):
        mp4 = run_ffmpeg_transcode(
            input_path,
            ".mp4",
            [
                "-c:v", "libx264",
                "-preset", "veryfast",
                "-crf", str(crf),
                "-pix_fmt", "yuv420p",
                "-movflags", "+faststart",
                "-c:a", "aac",
                "-b:a", "128k",
            ],
            timeout_sec=180,
        )
        if mp4:
            c = _candidate_from_bytes(
                f"h264_{tier.lower()}_crf{crf}", mp4, original_size, lossy=True,
                note="video H.264 lossy", codec_family="H.264", quality_tier=tier,
            )
            if c:
                out.append(c)

    # H.265/HEVC MP4 levels. More efficient than H.264, but slower and not
    # guaranteed to be previewable in every browser.
    for crf, tier in ((24, "Light"), (30, "Medium"), (36, "Strong")):
        mp4 = run_ffmpeg_transcode(
            input_path,
            ".mp4",
            [
                "-c:v", "libx265",
                "-preset", "ultrafast",
                "-crf", str(crf),
                "-pix_fmt", "yuv420p",
                "-tag:v", "hvc1",
                "-movflags", "+faststart",
                "-c:a", "aac",
                "-b:a", "128k",
            ],
            timeout_sec=180,
        )
        if mp4:
            c = _candidate_from_bytes(
                f"h265_{tier.lower()}_crf{crf}", mp4, original_size, lossy=True,
                note="video H.265 lossy", codec_family="H.265", quality_tier=tier,
            )
            if c:
                out.append(c)

    return out


# ============================================================
# Benchmark orchestration
# ============================================================

def _deduplicate_candidates(candidates: List[CompressionCandidate]) -> List[CompressionCandidate]:
    """Deduplicate by method, size, and sha1 prefix."""
    seen = set()
    unique: List[CompressionCandidate] = []

    for c in candidates:
        key = (c.method, c.size_bytes, hashlib.sha1(c.data[:4096]).hexdigest())
        if key in seen:
            continue
        seen.add(key)
        unique.append(c)

    return unique


def run_compression_benchmark(
    input_path: str,
    raw: bytes,
    quality_mode: str = "Compression",
) -> Tuple[CompressionCandidate, List[CompressionCandidate]]:
    """
    Build and rank all valid self-describing compression candidates.

    The smallest candidate by byte size is selected.

    quality_mode is retained for backward compatibility but should no longer
    control compressor inclusion. Compression mode tries all available candidates.
    """
    original_name = os.path.basename(input_path) or "input.bin"
    original_size = len(raw)
    domain = get_domain(input_path, raw)

    candidates: List[CompressionCandidate] = []

    # Domain-specific candidates.
    if domain == "image":
        candidates.extend(image_candidates(raw, original_name, original_size, quality_mode=quality_mode))
        candidates.extend(generic_container_candidates(raw, original_name, original_size))

    elif domain == "audio":
        candidates.extend(audio_compression_candidates(input_path, raw, original_size, quality_mode=quality_mode))
        candidates.extend(generic_container_candidates(raw, original_name, original_size))

    elif domain == "video":
        candidates.extend(video_compression_candidates(input_path, raw, original_size, quality_mode=quality_mode))
        candidates.extend(generic_container_candidates(raw, original_name, original_size))

    elif domain == "text":
        candidates.extend(text_compression_candidate(raw, original_name, original_size))

    elif domain in {"other", "unknown", "binary"}:
        candidates.extend(generic_container_candidates(raw, original_name, original_size))

    elif domain in {"archive", "document"}:
        candidates.extend(document_compression_candidates(raw, original_name, original_size))

    else:
        candidates.extend(generic_container_candidates(raw, original_name, original_size))

    # If a domain-specific encoder is unavailable, keep the app usable with a
    # small generic fallback set instead of showing no compression options.
    if not candidates:
        candidates.extend(fast_generic_container_candidates(raw, original_name, original_size))

    unique = _deduplicate_candidates(candidates)

    # Last-resort ZIP store: ensures one self-describing payload even for raw/unknown bytes.
    if not unique:
        b = zip_store_bytes(original_name, raw)
        c = _candidate_from_bytes("zip_store_fallback", b, original_size, lossy=False, note="fallback wrapper")
        if not c:
            raise RuntimeError("Could not create any self-describing compression candidate.")
        unique = [c]

    # Rank by smallest stored output; this is the byte stream that becomes DNA.
    unique.sort(key=lambda x: (x.size_bytes, x.method))

    for i, c in enumerate(unique, start=1):
        c.rank = i

    return unique[0], unique
