from __future__ import annotations
from pathlib import Path
from utils_bits_v2 import detect_magic
from config import IMAGE_KINDS, AUDIO_KINDS, VIDEO_KINDS

def fmt_bytes(n):
    if n is None: return "—"
    x=float(n)
    for unit in ["B","KB","MB","GB","TB"]:
        if x < 1024 or unit=="TB":
            return f"{int(x)} B" if unit=="B" else f"{x:.2f} {unit}"
        x/=1024

def get_domain(path: str, data: bytes) -> str:
    ext=Path(path).suffix.lower()
    document_exts={".pdf",".doc",".docx",".ppt",".pptx",".xls",".xlsx",".odt",".ods",".odp"}
    video_exts={".mp4",".mov",".m4v",".webm",".avi",".mkv",".wmv",".flv",".mpeg",".mpg",".3gp",".ts",".mts",".m2ts"}
    audio_exts={".wav",".mp3",".ogg",".opus",".flac",".m4a",".aac"}
    m=detect_magic(data)
    if not m:
        if ext in document_exts: return "document"
        if ext in video_exts: return "video"
        if ext in audio_exts: return "audio"
        try:
            data.decode("utf-8")
            return "text"
        except Exception:
            return "unknown"
    if m.kind in IMAGE_KINDS: return "image"
    if m.kind in AUDIO_KINDS: return "audio"
    if m.kind in VIDEO_KINDS or ext in video_exts: return "video"
    if m.kind in {"pdf","docx","pptx","xlsx","epub"}: return "document"
    if m.kind in {"zip","gzip","xz","bz2","zlib"}: return "archive"
    if m.kind == "text": return "text"
    if ext in document_exts: return "document"
    return "other"
