from __future__ import annotations
from exports import collect

import csv
import base64
import gzip
import hashlib
import html
import io
import json
import math
import mimetypes
import os
import random
import re
import subprocess
import tempfile
import zipfile
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List, Tuple

import pandas as pd
import numpy as np
import streamlit as st

try:
    from PIL import Image
except Exception:
    Image = None

import legacy_codec
import hedges_bridge as hb
import fastq_ecc
import randomness_ui
import compression_pipeline as app2_compression
from framing import bytes_to_bits, detect_and_validate_exact, restore_lossless_payload

APP_TITLE = "DNA Data Storage System"
DESIGNS = ["TVH-V5.2 HEDges"]
DESIGN_DISPLAY = {
    "Simple Mapping": "SM",
    "R0_B9": "R0-B9 (paper)",
    "R1_B12": "R1-B12 (paper)",
    "R2_B15": "R2-B15 (paper)",
    # Public name is concise; the manifest/internal ID retains shuffled-v1 so
    # legacy RINF_B16 files can never be decoded with the wrong rule ordering.
    "RINF_SHUFFLED_B16": "R∞",
    "TVH-V5.2 HEDges": "TVH-V5.2 HEDges",
}
FBR = hb.FBR
RBR = hb.RBR
TARGET_STRAND_NT = 125
BODY_NT = 76
LEGACY_SI_NT = 8

# -----------------------------------------------------------------------------
# Visual system — keeps the accepted App1 six-panel, bordered, two-column grammar.
# -----------------------------------------------------------------------------

def apply_style() -> None:
    st.markdown(
        """
<style>
.block-container {max-width: 1500px; padding-top: 1.1rem; padding-bottom: 4rem;}
html, body, [class*="css"] {font-family: Inter, "Source Sans 3", -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;}
.hero-card {border:1px solid rgba(15,23,42,.10); border-radius:18px; padding:1.05rem 1.25rem; margin-bottom:.8rem;
  background:linear-gradient(135deg, rgba(219,234,254,.45), rgba(240,253,250,.48));}
.hero-title {font-size:1.55rem; font-weight:800; letter-spacing:-.025em; color:#0f172a;}
.hero-subtitle {margin-top:.22rem; color:#475569; font-size:.96rem;}
.pipeline-steps {display:grid; grid-template-columns:repeat(6,1fr); gap:.45rem; margin:.55rem 0 1.0rem 0;}
.pipeline-step {border:1px solid rgba(100,116,139,.22); border-radius:12px; padding:.55rem .65rem; background:rgba(248,250,252,.7);}
.pipeline-step.done {background:rgba(220,252,231,.65); border-color:rgba(34,197,94,.28);}
.pipeline-step.current {background:rgba(219,234,254,.7); border-color:rgba(59,130,246,.32);}
.step-num {display:inline-grid; place-items:center; width:1.35rem; height:1.35rem; border-radius:999px; background:#e2e8f0; font-weight:800; font-size:.76rem; margin-right:.38rem;}
.step-name {font-weight:700; font-size:.78rem; color:#334155;}
.step-state {font-size:.68rem; color:#64748b; margin-top:.18rem;}
.seq-box {border:1px solid rgba(100,116,139,.24); border-radius:10px; background:#f8fafc; padding:.65rem .75rem; max-height:220px; overflow:auto;
  font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace; font-size:.77rem; line-height:1.48; word-break:break-all; white-space:pre-wrap;}
.segment {display:inline-block; border-radius:4px; padding:1px 3px; margin:1px; font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace; font-size:.77rem; word-break:break-all;}
.seg-fbr{background:#dbeafe;color:#1e3a8a}.seg-si{background:#fef3c7;color:#78350f}.seg-body{background:#dcfce7;color:#14532d}.seg-fill{background:#f1f5f9;color:#475569}.seg-rbr{background:#ede9fe;color:#4c1d95}
.err-sub{background:#fecaca;color:#7f1d1d;font-weight:800;border-radius:4px;padding:1px 3px}.err-ins{background:#e9d5ff;color:#581c87;font-weight:800;border-radius:4px;padding:1px 3px}.err-del{background:#fed7aa;color:#7c2d12;font-weight:800;border-radius:4px;padding:1px 3px;text-decoration:line-through}.err-legend{display:inline-block;border-radius:5px;padding:2px 7px;margin-right:6px;font-size:.76rem;font-weight:700}
.small-note {font-size:.82rem; color:#64748b;}
div[data-testid="stMetric"] {border:1px solid rgba(100,116,139,.17); border-radius:11px; padding:.7rem .78rem;}
@media (max-width:900px){.pipeline-steps{grid-template-columns:repeat(2,1fr)}}
</style>
""",
        unsafe_allow_html=True,
    )


def step_header(no: int, title: str) -> None:
    st.markdown(f"### {no}. {title}")


def _step_done(no: int) -> bool:
    checks = {
        1: bool(st.session_state.get("input_bytes")),
        2: bool(st.session_state.get("stored_bytes")),
        3: bool(st.session_state.get("encoded_ready")),
        4: bool(st.session_state.get("strand_rows")),
        5: bool(st.session_state.get("decode_result")),
        6: bool(st.session_state.get("decode_result")),
    }
    return bool(checks.get(no))


def render_stepper() -> None:
    labels = ["Input", "Compression", "Encoding", "Strand Preparation", "Decoding", "Summarization"]
    parts = ['<div class="pipeline-steps">']
    for i, label in enumerate(labels, 1):
        done = _step_done(i)
        prior = all(_step_done(x) for x in range(1, i)) if i > 1 else True
        css = "done" if done else ("current" if prior else "")
        state = "Done" if done else ("Next" if prior else "Waiting")
        parts.append(f'<div class="pipeline-step {css}"><div><span class="step-num">{i}</span><span class="step-name">{label}</span></div><div class="step-state">{state}</div></div>')
    parts.append("</div>")
    st.markdown("".join(parts), unsafe_allow_html=True)


def seq_box(text: str, limit: int = 8000) -> None:
    text = str(text or "")
    shown = text[:limit] + ("\n… [preview truncated]" if len(text) > limit else "")
    st.markdown(f'<div class="seq-box">{html.escape(shown)}</div>', unsafe_allow_html=True)


def property_table(rows: List[Dict[str, Any]]) -> None:
    if rows:
        frame = pd.DataFrame(rows)
        for column in frame.select_dtypes(include="object"):
            frame[column] = frame[column].map(lambda value: "—" if value is None else str(value))
        st.dataframe(frame, width="stretch", hide_index=True)


def _df_csv(rows: List[Dict[str, Any]]) -> bytes:
    if not rows:
        return b""
    buf = io.StringIO()
    writer = csv.DictWriter(buf, fieldnames=list(rows[0].keys()))
    writer.writeheader(); writer.writerows(rows)
    return buf.getvalue().encode("utf-8")


def _reset_after_input() -> None:
    for key in ('process_run','sim_decode','sim_public_signature'):
        st.session_state.pop(key,None)
    for k in ["stored_bytes","stored_name","stored_mime","compression_method","compression_rows","compression_analysis_rows","compression_selected_fingerprint","compression_version_choice","compression_version_selector","encoded_ready","dna_design","legacy_dna","legacy_bits","legacy_meta","hedges_strands","hedges_meta","strand_rows","direct_error_rows","direct_error_stats","fastq_bytes","fastq_sim_stats","fastq_origin","consensus_map","consensus_rows","consensus_position_rows","consensus_stats","consensus_eval_rows","consensus_eval_positions","consensus_eval_stats","consensus_sid_mode","decode_result","encoding_rand_overall","encoding_rand_rows","encoding_rand_hash"]:
        st.session_state.pop(k, None)


def _reset_after_compression() -> None:
    for key in ('process_run','sim_decode','sim_public_signature'):
        st.session_state.pop(key,None)
    for k in ["encoded_ready","dna_design","legacy_dna","legacy_bits","legacy_meta","hedges_strands","hedges_meta","strand_rows","direct_error_rows","direct_error_stats","fastq_bytes","fastq_sim_stats","fastq_origin","consensus_map","consensus_rows","consensus_position_rows","consensus_stats","consensus_eval_rows","consensus_eval_positions","consensus_eval_stats","consensus_sid_mode","decode_result","encoding_rand_overall","encoding_rand_rows","encoding_rand_hash"]:
        st.session_state.pop(k, None)


def _reset_after_encoding() -> None:
    for key in ('process_run','sim_decode','sim_public_signature'):
        st.session_state.pop(key,None)
    for k in ["strand_rows","direct_error_rows","direct_error_stats","fastq_bytes","fastq_sim_stats","fastq_origin","consensus_map","consensus_rows","consensus_position_rows","consensus_stats","consensus_eval_rows","consensus_eval_positions","consensus_eval_stats","consensus_sid_mode","decode_result","encoding_rand_overall","encoding_rand_rows","encoding_rand_hash"]:
        st.session_state.pop(k, None)


def _detect_type(data: bytes, name: str = "") -> tuple[str, str, str]:
    fmt, ext, mime, strong, note = detect_and_validate_exact(bytes(data or b""))
    if fmt != "Unknown binary":
        return fmt, ext, mime
    guessed = mimetypes.guess_type(name)[0] or "application/octet-stream"
    return ("Text" if guessed.startswith("text/") else "Binary"), Path(name).suffix or ".bin", guessed


def _ooxml_text_preview(data: bytes, name: str) -> str | None:
    ext=Path(name).suffix.lower()
    if ext not in {".docx",".pptx",".xlsx",".odt",".ods",".odp"}:
        return None
    try:
        with zipfile.ZipFile(io.BytesIO(data),"r") as archive:
            members=archive.namelist()
            if ext==".docx": selected=[x for x in members if x=="word/document.xml"]
            elif ext==".pptx": selected=sorted(x for x in members if re.fullmatch(r"ppt/slides/slide\d+\.xml",x))[:20]
            elif ext==".xlsx": selected=[x for x in members if x=="xl/sharedStrings.xml"]+sorted(x for x in members if x.startswith("xl/worksheets/sheet"))[:5]
            else: selected=[x for x in members if x=="content.xml"]
            values=[]
            for member in selected:
                root=ET.fromstring(archive.read(member))
                values.extend((node.text or "").strip() for node in root.iter() if node.tag.rsplit("}",1)[-1] in {"t","p"} and (node.text or "").strip())
            return "\n".join(values[:1000]) if values else None
    except Exception:
        return None


def _preview_bytes(data: bytes, name: str, mime: str | None = None, title: str = "Preview", key: str | None = None) -> None:
    st.markdown(f"#### {title}")
    mime = mime or mimetypes.guess_type(name)[0] or "application/octet-stream"
    try:
        if mime.startswith("image/"):
            st.image(data, caption=name, width="stretch")
        elif mime.startswith("audio/"):
            st.audio(data, format=mime)
        elif mime.startswith("video/"):
            st.video(data, format=mime)
        elif mime == "application/pdf" or str(name).lower().endswith(".pdf"):
            if hasattr(st, "pdf"):
                try:
                    st.pdf(data, height=520)
                    return
                except Exception:
                    pass
            encoded = base64.b64encode(data).decode("ascii")
            import streamlit.components.v1 as components
            components.html(
                f'<iframe src="data:application/pdf;base64,{encoded}" width="100%" height="520" style="border:0"></iframe>',
                height=530,
            )
        elif Path(str(name)).suffix.lower() in {".docx",".pptx",".xlsx",".odt",".ods",".odp"}:
            extracted=_ooxml_text_preview(data,name)
            if extracted:
                st.text_area("Document content preview",extracted[:12000],height=260,disabled=True,key=key or f"office_{hashlib.sha1(data[:512]).hexdigest()[:12]}")
            else:
                st.info("Document structure is valid, but no browser-safe text preview was available.")
        elif mime.startswith("text/"):
            try:
                txt = data.decode("utf-8")
            except Exception:
                txt = data.decode("utf-8", errors="replace")
            if key is None:
                token = hashlib.sha1(
                    (str(title) + "|" + str(name) + "|" + str(mime) + "|" + str(len(data))).encode("utf-8")
                    + bytes(data[:128])
                ).hexdigest()[:16]
                key = f"text_preview_{token}"
            st.text_area("Text preview", txt[:10000], height=220, disabled=True, key=key)
        else:
            st.info("Binary/container file. Preview is available after reconstruction when the container is supported by Streamlit.")
    except Exception as exc:
        st.warning(f"Preview unavailable: {exc}")


# -----------------------------------------------------------------------------
# Compression
# -----------------------------------------------------------------------------

def _zip_bytes(data: bytes, filename: str) -> bytes:
    bio = io.BytesIO()
    with zipfile.ZipFile(bio, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as zf:
        zf.writestr(Path(filename).name or "payload.bin", data)
    return bio.getvalue()


def _image_convert(data: bytes, fmt: str, quality: int | None = None) -> bytes | None:
    if Image is None:
        return None
    try:
        im = Image.open(io.BytesIO(data))
        out = io.BytesIO()
        if fmt == "PNG":
            im.save(out, format="PNG", compress_level=9)
        elif fmt == "JPEG":
            if im.mode not in ("RGB", "L"):
                im = im.convert("RGB")
            im.save(out, format="JPEG", quality=int(quality or 90), optimize=True)
        elif fmt == "WEBP":
            im.save(out, format="WEBP", quality=int(quality or 90), method=6)
        else:
            return None
        return out.getvalue()
    except Exception:
        return None


@st.cache_data(show_spinner=False)
def compression_candidates(data: bytes, filename: str) -> List[dict]:
    """Use the complete accepted App2 domain-specific compression benchmark.

    No-compression is kept as an explicit baseline. All other candidates are
    generated by the original App2 compressor families and settings.
    """
    rows = [{
        "Method": "No compression",
        "data": data,
        "name": filename,
        "mime": mimetypes.guess_type(filename)[0] or "application/octet-stream",
        "lossy": False,
        "note": "baseline",
        "codec": "Original",
        "tier": "Exact",
    }]
    with tempfile.TemporaryDirectory() as td:
        original_name = Path(filename).name or "payload.bin"
        tmp_path = str(Path(td) / original_name)
        Path(tmp_path).write_bytes(data)
        _best, candidates = app2_compression.run_compression_benchmark(tmp_path, data)
        stem = Path(filename).stem or "stored"
        for c in candidates:
            rows.append({
                "Method": c.method,
                "data": c.data,
                "name": stem + c.ext,
                "mime": c.mime,
                "lossy": bool(c.lossy),
                "note": c.note,
                "kind": c.kind,
                "codec": c.codec_family,
                "tier": c.quality_tier,
            })
    # Avoid duplicate byte-identical/method-identical entries while preserving App2 order.
    seen=set(); out=[]
    for r in rows:
        key=(r["Method"], hashlib.sha1(r["data"]).hexdigest())
        if key in seen: continue
        seen.add(key); out.append(r)
    return out


def _format_method_name(method: str) -> str:
    names = {
        "No compression": "Original file",
        "brotli_quality11": "Brotli · exact",
        "zstd_level19": "Zstandard · exact",
        "xz_preset9": "XZ/LZMA2 · exact",
        "zlib_level9": "zlib · exact",
        "png_optimized": "PNG optimized · exact pixels",
        "webp_lossless": "WebP · exact pixels",
        "flac_lossless": "FLAC · exact samples",
    }
    value = str(method or "No compression")
    if value in names:
        return names[value]
    return value.replace("_", " ").replace("light", "Light").replace("medium", "Medium").replace("strong", "Strong")


def _candidate_preview(cand: dict) -> tuple[bytes, str, str]:
    """Return the human-viewable representation, unwrapping exact compressors."""
    data = bytes(cand.get("data") or b"")
    restored = restore_lossless_payload(data)
    if restored is not None:
        return restored["data"], restored["name"], restored["mime"]
    return data, str(cand.get("name") or "preview.bin"), str(cand.get("mime") or "application/octet-stream")


def _chart_lollipop(rows: List[Dict[str, Any]], label_col: str, value_col: str, *, key: str) -> None:
    if not rows:
        return
    df = pd.DataFrame(rows).copy()
    if label_col not in df or value_col not in df:
        return
    df["Value"] = pd.to_numeric(df[value_col], errors="coerce")
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["Value"])
    if df.empty:
        return
    df["Order"] = range(len(df)); df["Zero"] = 0.0
    df["Label"] = df["Value"].map(lambda x: f"{float(x):.3f}")
    spec = {
        "height": max(140, min(310, 35 * len(df) + 30)),
        "encoding": {
            "y": {"field": label_col, "type": "nominal", "sort": {"field":"Order","order":"ascending"}, "axis":{"labelLimit":220,"title":None}},
            "x": {"field":"Value","type":"quantitative","title":value_col,"scale":{"zero":True},"axis":{"grid":True}},
        },
        "layer": [
            {"mark":{"type":"rule","strokeWidth":2,"color":"#94a3b8"},"encoding":{"x":{"field":"Zero","type":"quantitative"},"x2":{"field":"Value"}}},
            {"mark":{"type":"point","filled":True,"size":75,"color":"#2563eb"}},
            {"mark":{"type":"text","align":"left","baseline":"middle","dx":8,"fontSize":11,"color":"#0f172a"},"encoding":{"text":{"field":"Label","type":"nominal"}}},
        ],
    }
    st.vega_lite_chart(df, spec, width="stretch", key=key)


@st.cache_data(show_spinner=False)
def _image_quality_metrics(original: bytes, candidate: bytes) -> Dict[str, float] | None:
    if Image is None:
        return None
    try:
        ref = Image.open(io.BytesIO(original)).convert("RGB")
        got = Image.open(io.BytesIO(candidate)).convert("RGB")
        if got.size != ref.size:
            return None  # Never resize a candidate silently to improve its score.
        a = np.asarray(ref, dtype=np.float64)
        b = np.asarray(got, dtype=np.float64)
        mae = float(np.mean(np.abs(a-b)))
        mse = float(np.mean((a-b)**2))
        psnr = float("inf") if mse == 0 else float(20.0*np.log10(255.0/np.sqrt(mse)))
        window=min(7,a.shape[0],a.shape[1])
        if window % 2 == 0:window-=1
        ssim=None
        if window>=3:
            from scipy.ndimage import uniform_filter
            vals=[]; C1=(0.01*255)**2; C2=(0.03*255)**2
            pad=window//2; covariance_scale=window**2/(window**2-1)
            for ch in range(3):
                x=a[:,:,ch]; y=b[:,:,ch]
                mux=uniform_filter(x,size=window); muy=uniform_filter(y,size=window)
                vx=covariance_scale*(uniform_filter(x*x,size=window)-mux*mux)
                vy=covariance_scale*(uniform_filter(y*y,size=window)-muy*muy)
                cov=covariance_scale*(uniform_filter(x*y,size=window)-mux*muy)
                score=((2*mux*muy+C1)*(2*cov+C2))/((mux*mux+muy*muy+C1)*(vx+vy+C2))
                vals.append(float(score[pad:-pad,pad:-pad].mean()))
            ssim=float(np.mean(vals))
        return {"PSNR":psnr,"SSIM":ssim,"MAE":mae}
    except Exception:
        return None


def _safe_media_suffix(name: str, fallback: str) -> str:
    suffix = Path(str(name or "")).suffix.lower()
    return suffix if re.fullmatch(r"\.[a-z0-9]{1,8}", suffix or "") else fallback


def _decode_audio_pcm(data: bytes, suffix: str) -> np.ndarray | None:
    try:
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / ("audio" + _safe_media_suffix(suffix, ".bin"))
            path.write_bytes(data)
            proc = subprocess.run(
                ["ffmpeg", "-hide_banner", "-loglevel", "error", "-i", str(path), "-vn", "-ac", "1", "-ar", "48000", "-f", "f32le", "-"],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=90,
            )
            if proc.returncode != 0 or not proc.stdout:
                return None
            return np.frombuffer(proc.stdout, dtype="<f4").astype(np.float64)
    except Exception:
        return None


@st.cache_data(show_spinner=False)
def _audio_quality_metrics(original: bytes, original_name: str, candidate: bytes, candidate_name: str) -> Dict[str, float] | None:
    ref = _decode_audio_pcm(original, original_name)
    got = _decode_audio_pcm(candidate, candidate_name)
    if ref is None or got is None:
        return None
    n = min(len(ref), len(got))
    if n < 32:
        return None
    ref = ref[:n]; got = got[:n]
    error = ref - got
    signal_rms = float(np.sqrt(np.mean(ref * ref)))
    error_rms = float(np.sqrt(np.mean(error * error)))
    snr = float("inf") if error_rms <= 1e-15 else float(20 * np.log10(max(signal_rms, 1e-15) / error_rms))
    denom = float(np.linalg.norm(ref) * np.linalg.norm(got))
    correlation = 1.0 if np.array_equal(ref, got) else (float(np.dot(ref, got) / denom) if denom > 0 else 0.0)
    nrmse = error_rms / max(signal_rms, 1e-15)
    return {"Audio SNR": snr, "Audio correlation": correlation, "Audio NRMSE": nrmse}


def _ffmpeg_video_metric(original: bytes, original_name: str, candidate: bytes, candidate_name: str, metric: str) -> float | None:
    try:
        with tempfile.TemporaryDirectory() as td:
            ref = Path(td) / ("reference" + _safe_media_suffix(original_name, ".mp4"))
            got = Path(td) / ("candidate" + _safe_media_suffix(candidate_name, ".mp4"))
            ref.write_bytes(original); got.write_bytes(candidate)
            proc = subprocess.run(
                ["ffmpeg", "-hide_banner", "-i", str(ref), "-i", str(got), "-t", "20", "-lavfi", f"[0:v][1:v]{metric}", "-f", "null", "-"],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=120,
            )
            log = proc.stderr or ""
            if metric == "ssim":
                hits = re.findall(r"All:([0-9.]+)", log)
            else:
                hits = re.findall(r"average:([0-9.]+|inf)", log)
            if not hits:
                return None
            return float("inf") if hits[-1] == "inf" else float(hits[-1])
    except Exception:
        return None


@st.cache_data(show_spinner=False)
def _video_quality_metrics(original: bytes, original_name: str, candidate: bytes, candidate_name: str) -> Dict[str, float] | None:
    ssim = _ffmpeg_video_metric(original, original_name, candidate, candidate_name, "ssim")
    psnr = _ffmpeg_video_metric(original, original_name, candidate, candidate_name, "psnr")
    if ssim is None and psnr is None:
        return None
    out: Dict[str, float] = {}
    if ssim is not None: out["Video SSIM"] = ssim
    if psnr is not None: out["Video PSNR"] = psnr
    return out


def _compression_analysis_rows(data: bytes, cands: List[dict]) -> List[dict]:
    ranked = sorted(cands, key=lambda x: len(x["data"]))
    out=[]
    input_mime=(st.session_state.get("input_mime") or mimetypes.guess_type(st.session_state.get("input_name", ""))[0] or "application/octet-stream")
    original_name=st.session_state.get("input_name", "input.bin")
    input_is_image=input_mime.startswith("image/")
    input_is_audio=input_mime.startswith("audio/")
    input_is_video=input_mime.startswith("video/")
    for rank,c in enumerate(ranked,1):
        preview_data, preview_name, preview_mime = _candidate_preview(c)
        q=None
        if not c.get("lossy"):
            if input_is_image:
                q=_image_quality_metrics(data,preview_data)
            elif input_is_audio:
                q={"Audio SNR":float("inf"),"Audio correlation":1.0,"Audio NRMSE":0.0}
            elif input_is_video:
                q={"Video SSIM":1.0,"Video PSNR":float("inf")}
        elif input_is_image and preview_mime.startswith("image/"):
            q=_image_quality_metrics(data,preview_data)
        elif input_is_audio and preview_mime.startswith("audio/"):
            q=_audio_quality_metrics(data,original_name,preview_data,preview_name)
        elif input_is_video and preview_mime.startswith("video/"):
            q=_video_quality_metrics(data,original_name,preview_data,preview_name)
        row={
            "Rank":rank, "Method":c["Method"], "Format":Path(c["name"]).suffix.lower() or "—",
            "Codec":c.get("codec") or c["Method"].split("_",1)[0],
            "Quality level":c.get("tier") or ("Exact" if not c.get("lossy") else "—"),
            "Stored bytes":len(c["data"]), "Compression ratio":len(data)/max(1,len(c["data"])),
            "Size reduction (%)":(1-len(c["data"])/max(1,len(data)))*100.0,
            "Estimated full DNA nt (HEDges)":hb.estimate_strands(c["data"],c.get("name","stored.bin"),c.get("mime"))*125,
            "Lossy":"Yes" if c.get("lossy") else "No",
        }
        if q:
            row.update(q)
            if q.get("SSIM") is not None: row["Quality (%)"]=100.0*float(q["SSIM"])
            elif q.get("Audio correlation") is not None: row["Quality (%)"]=100.0*float(q["Audio correlation"])
            elif q.get("Video SSIM") is not None: row["Quality (%)"]=100.0*float(q["Video SSIM"])
        elif not c.get("lossy"):
            row["Quality (%)"]=100.0
        out.append(row)
    return out


def _display_compression_table(rows: List[dict], selected: str) -> None:
    shown=[]
    for r in rows:
        psnr=r.get("PSNR") if r.get("PSNR") is not None else r.get("Video PSNR")
        shown.append({
            "Selected":"✓" if r["Method"]==selected else "", "Rank":r["Rank"], "Method":_format_method_name(r["Method"]),
            "Codec":r.get("Codec","—"), "Level":r.get("Quality level","—"), "Format":r["Format"],
            "Stored data":f"{r['Stored bytes']:,} B", "Compression ratio":f"{r['Compression ratio']:.3f}×",
            "Size reduction":f"{r['Size reduction (%)']:.2f}%", "Quality":f"{r['Quality (%)']:.2f}%" if r.get("Quality (%)") is not None else "—",
            "PSNR":"∞" if psnr==float('inf') else (f"{psnr:.2f} dB" if psnr is not None else "—"),
            "SSIM":f"{r['SSIM']:.4f}" if r.get("SSIM") is not None else (f"{r['Video SSIM']:.4f}" if r.get("Video SSIM") is not None else "—"),
            "Audio SNR":"∞" if r.get("Audio SNR")==float("inf") else (f"{r['Audio SNR']:.2f} dB" if r.get("Audio SNR") is not None else "—"),
            "Est. HEDges nt":f"{r['Estimated full DNA nt (HEDges)']:,}",
        })
    st.dataframe(pd.DataFrame(shown),width="stretch",hide_index=True)


def _compression_quality_fields(rows):
    if any(r.get('PSNR') is not None for r in rows):
        return [('PSNR','PSNR (dB) ↑'),('SSIM','SSIM ↑'),('MAE','MAE ↓')]
    if any(r.get('Video PSNR') is not None or r.get('Video SSIM') is not None for r in rows):
        return [('Video PSNR','PSNR (dB) ↑'),('Video SSIM','SSIM ↑')]
    if any(r.get('Audio correlation') is not None for r in rows):
        return [('Audio SNR','SNR (dB) ↑'),('Audio correlation','Correlation ↑'),('Audio NRMSE','NRMSE ↓')]
    return []


def _compression_quality_cards(row,rows):
    fields=_compression_quality_fields(rows)
    if not fields:return
    st.markdown('##### Compression quality')
    for col,(field,label) in zip(st.columns(len(fields)),fields):
        value=(row or {}).get(field)
        text='—' if value is None else '∞' if value==float('inf') else f'{value:.4f}' if 'SSIM' in field or 'correlation' in field or 'NRMSE' in field else f'{value:.2f}'
        col.metric(label,text,help='Compared with the original decoded content. ↑ higher is better; ↓ lower is better. ∞ means zero measured difference.')


def _compression_metric_chart(rows,field,label,selected,key):
    records=[{'Method':_format_method_name(r['Method']),'Value':float(r[field]),
              'Selected':r['Method']==selected,'Stored KB':r['Stored bytes']/1024}
             for r in rows if r.get(field) is not None and math.isfinite(r[field])]
    if not records:
        st.info('All measured candidates are exact (∞ dB).' if any(r.get(field)==float('inf') for r in rows) else 'This metric is unavailable for these candidates.')
        return
    spec={'height':max(220,28*len(records)), 'mark':{'type':'bar','cornerRadiusEnd':4},
          'encoding':{'y':{'field':'Method','type':'nominal','sort':None,'title':None,'axis':{'labelLimit':230}},
                      'x':{'field':'Value','type':'quantitative','title':label},
                      'color':{'condition':{'test':'datum.Selected','value':'#16a34a'},'value':'#60a5fa'},
                      'tooltip':[{'field':'Method'},{'field':'Value','format':'.4f'},{'field':'Stored KB','format':'.2f'},{'field':'Selected'}]}}
    st.vega_lite_chart(pd.DataFrame(records),spec,width='stretch',key=key)


def _render_compression_charts(rows: List[dict],selected='') -> None:
    if not rows:
        return
    st.markdown("#### Compression comparison")
    fields=_compression_quality_fields(rows)
    if fields:
        quality_tab,tradeoff_tab,size_tab=st.tabs(['Quality comparison','Size vs quality','File size'])
        with quality_tab:
            left,right=st.columns(2)
            for col,(field,label) in zip((left,right),fields[:2]):
                with col:
                    st.markdown('##### '+label)
                    _compression_metric_chart(rows,field,label,selected,'compression_'+field.replace(' ','_'))
            if any(r.get(fields[0][0])==float('inf') for r in rows):st.caption('∞ dB = identical pixels or samples; infinite values are omitted from the dB chart. Green = selected output.')
    else:
        tradeoff_tab=st.container();size_tab=None
    chart = pd.DataFrame([{
        "Method":_format_method_name(r["Method"]), "Codec":r.get("Codec","—"),
        "Level":r.get("Quality level","—"), "Stored bytes":r["Stored bytes"],
        "Quality (%)":r.get("Quality (%)"), "Selected":r['Method']==selected,
    } for r in rows])
    if fields and chart["Quality (%)"].notna().any():
        spec={
            "height":330,
            "mark":{"type":"circle","size":150,"opacity":0.85},
            "encoding":{
                "x":{"field":"Stored bytes","type":"quantitative","title":"Stored bytes (smaller is better)","scale":{"zero":False}},
                "y":{"field":"Quality (%)","type":"quantitative","title":"SSIM × 100" if 'SSIM' in fields[1][0] else "Correlation × 100","scale":{"zero":False}},
                "color":{"field":"Codec","type":"nominal"},
                "shape":{"field":"Level","type":"nominal"},
                "stroke":{"condition":{"test":"datum.Selected","value":"#111827"},"value":"transparent"},
                "strokeWidth":{"value":3},
                "tooltip":[{"field":"Method"},{"field":"Level"},{"field":"Stored bytes"},{"field":"Quality (%)","format":".2f"}],
            },
        }
        with tradeoff_tab:st.vega_lite_chart(chart.dropna(subset=["Quality (%)"]),spec,width="stretch",key="compression_tradeoff_chart")
    with size_tab if size_tab is not None else tradeoff_tab:
        _chart_lollipop([{"Method":_format_method_name(r['Method']),"Stored bytes":r['Stored bytes']} for r in rows],"Method","Stored bytes",key="compression_size_chart")
    if fields:
        with st.expander('Metric definitions'):
            st.write('PSNR: pixel fidelity in dB. SSIM: structural similarity (1 = identical). MAE: mean absolute RGB error on the 0–255 scale. These are measured after decoding the compressed file; they do not measure DNA recovery.')
            st.write('Images use local-window SSIM, up to 7 × 7, averaged over RGB. No resizing is applied. SSIM is unavailable for images smaller than 3 pixels on either axis. Audio uses decoded-waveform SNR, correlation and NRMSE; video uses FFmpeg PSNR/SSIM when available.')


# -----------------------------------------------------------------------------
# Legacy strand design
# -----------------------------------------------------------------------------

def _base4_si(sid: int, width: int = 8) -> str:
    chars = "ACGT"; n = int(sid); out = ["A"] * width
    for i in range(width-1, -1, -1):
        out[i] = chars[n % 4]; n //= 4
    if n:
        raise ValueError("Strand ID exceeds 8-nt base-4 index capacity")
    return "".join(out)


def _safe_filler(n: int, sid: int) -> str:
    motif = "ACGT" if sid % 2 == 0 else "TGCA"
    return "".join(motif[(i + sid) % 4] for i in range(n))


def build_legacy_strands(dna: str) -> List[dict]:
    dna = legacy_codec.clean_dna_text(dna)
    rows = []
    for sid, start in enumerate(range(0, len(dna), BODY_NT)):
        payload = dna[start:start+BODY_NT]
        filler = _safe_filler(BODY_NT-len(payload), sid)
        si = _base4_si(sid, LEGACY_SI_NT)
        full = FBR + si + payload + filler + RBR
        rows.append({
            "Strand ID": sid, "FBR":FBR, "SI":si, "Payload":payload,
            "Payload length":len(payload), "Filler":filler, "RBR":RBR,
            "Full strand":full, "Full strand length":len(full),
            "GC content":legacy_codec.gc_content(full),
            "Longest homopolymer":legacy_codec.longest_homopolymer(full),
        })
    return rows


def show_segmented(row: dict) -> None:
    st.markdown(
        '<div class="seq-box">'
        f'<span class="segment seg-fbr">FBR {html.escape(str(row.get("FBR","")))}</span>'
        f'<span class="segment seg-si">SI {html.escape(str(row.get("SI","")))}</span>'
        f'<span class="segment seg-body">PAYLOAD {html.escape(str(row.get("Payload","")))}</span>'
        + (f'<span class="segment seg-fill">FILLER {html.escape(str(row.get("Filler","")))}</span>' if row.get("Filler") else "") +
        f'<span class="segment seg-rbr">RBR {html.escape(str(row.get("RBR","")))}</span>'
        '</div>', unsafe_allow_html=True
    )


def _legacy_payload_from_consensus(seq: str, template: dict) -> str | None:
    seq = fastq_ecc.clean_dna(seq)
    # Consensus should restore the fixed-length full strand before legacy decode.
    if len(seq) != TARGET_STRAND_NT:
        return None
    start = len(FBR) + LEGACY_SI_NT
    n = int(template.get("Payload length", BODY_NT))
    return seq[start:start+n]


# -----------------------------------------------------------------------------
# Panel 1 — Input
# -----------------------------------------------------------------------------

def render_panel_1_upload() -> None:
    with st.container(border=True):
        step_header(1, "Input")
        left, right = st.columns(2, gap="large")
        with left:
            st.markdown("#### 📁 Source file")
            up = st.file_uploader("Upload file", type=None, key="input_upload")
            if up is not None:
                data = up.getvalue(); digest = hashlib.sha256(data).hexdigest()
                if st.session_state.get("input_sha") != digest:
                    _reset_after_input()
                    st.session_state.update({
                        "input_bytes":data, "input_name":up.name,
                        "input_mime":up.type or mimetypes.guess_type(up.name)[0] or "application/octet-stream",
                        "input_sha":digest,
                    })
            data = st.session_state.get("input_bytes", b"")
            name = st.session_state.get("input_name", "")
            if data:
                fmt, ext, mime = _detect_type(data, name)
                property_table([
                    {"Property":"File name","Value":name},
                    {"Property":"Detected type","Value":fmt},
                    {"Property":"File size","Value":f"{len(data):,} B"},
                    {"Property":"Binary length","Value":f"{len(data)*8:,} bits"},
                    {"Property":"SHA-256","Value":hashlib.sha256(data).hexdigest()},
                ])
        with right:
            data = st.session_state.get("input_bytes", b"")
            if data:
                _preview_bytes(data, st.session_state.get("input_name","input.bin"), st.session_state.get("input_mime"), "Uploaded file", key="input_text_preview")
                st.markdown("#### 🔢 Binary representation")
                bits=bytes_to_bits(data)
                seq_box(bits,6000)
                collect("Download binary (.txt)",bits.encode("ascii"),f"{Path(st.session_state.get('input_name','input')).stem}_binary.txt","text/plain",width="stretch",key="download_input_binary")
            else:
                st.info("Upload a file to start the pipeline.")


# -----------------------------------------------------------------------------
# Panel 2 — Compression
# -----------------------------------------------------------------------------

def _compressed_candidate_choice(cands: List[dict], choice: str) -> dict:
    compressed = [c for c in cands if c.get("Method") != "No compression"]
    if not compressed:
        raise ValueError("No valid compression candidate is available for this file type")
    if choice == "Auto best (smallest)":
        return min(compressed, key=lambda c: (len(c["data"]), c["Method"]))
    for c in compressed:
        label = f"{_format_method_name(c['Method'])} | {len(c['data']):,} B"
        if label == choice:
            return c
    raise ValueError(f"Unknown compressed version: {choice}")


def _balanced_candidate_choice(cands: List[dict]) -> dict:
    compressed = [c for c in cands if c.get("Method") != "No compression"]
    if not compressed:
        raise ValueError("No valid compression candidate is available for this file type")
    balanced = [c for c in compressed if c.get("tier") in {"Medium", "Exact", ""}]
    return min(balanced or compressed, key=lambda c:(len(c["data"]),c["Method"]))


def _apply_selected_compression(cand: dict) -> None:
    fingerprint = hashlib.sha256(cand["data"]).hexdigest() + "|" + str(cand["Method"])
    if st.session_state.get("compression_selected_fingerprint") == fingerprint:
        return
    _reset_after_compression()
    st.session_state.update({
        "stored_bytes": cand["data"],
        "stored_name": cand["name"],
        "stored_mime": cand["mime"],
        "compression_method": cand["Method"],
        "compression_lossy": cand["lossy"],
        "compression_selected_fingerprint": fingerprint,
    })


def _expected_compressor_methods(data: bytes, filename: str) -> List[str]:
    """Return the App2 compressor/version menu before the benchmark is run.

    This is UI-only discovery. The actual candidate must still be produced and
    validated by compression_candidates(); unavailable media codecs fail closed.
    """
    mime = (mimetypes.guess_type(filename)[0] or "application/octet-stream").lower()
    ext = Path(filename).suffix.lower()
    if mime.startswith("image/") or ext in {".png", ".jpg", ".jpeg", ".webp", ".bmp", ".tif", ".tiff"}:
        return ["png_optimized", "webp_lossless", "jpeg_light_q90", "jpeg_medium_q70", "jpeg_strong_q10", "webp_light_q90", "webp_medium_q70", "webp_strong_q10", "zstd_level19", "xz_preset9", "zlib_level9"]
    if mime.startswith("audio/") or ext in {".wav", ".mp3", ".flac", ".aac", ".m4a", ".ogg", ".opus"}:
        return ["flac_lossless", "opus_light_128k", "opus_medium_64k", "opus_strong_32k", "zstd_level19", "xz_preset9", "zlib_level9"]
    if mime.startswith("video/") or ext in {".mp4", ".mov", ".avi", ".mkv", ".webm"}:
        return ["h264_light_crf20", "h264_medium_crf28", "h264_strong_crf36", "h265_light_crf24", "h265_medium_crf30", "h265_strong_crf36", "zstd_level19", "xz_preset9", "zlib_level9"]
    if mime.startswith("text/") or ext in {".txt", ".csv", ".json", ".xml", ".md", ".py", ".m", ".html", ".css", ".js"}:
        return ["brotli_quality11", "zstd_level19", "xz_preset9"]
    if ext in {".pdf", ".doc", ".docx", ".ppt", ".pptx", ".xls", ".xlsx", ".epub"}:
        return ["zstd_level19", "xz_preset9", "zlib_level9"]
    return ["zstd_level19", "xz_preset9", "zlib_level9"]


def _candidate_by_method(cands: List[dict], method: str) -> dict:
    for c in cands:
        if c.get("Method") == method:
            return c
    available = [str(c.get("Method")) for c in cands if c.get("Method") != "No compression"]
    raise ValueError(
        f"Selected compressor/version '{method}' was not available on this system. "
        f"Available: {', '.join(available) if available else 'none'}"
    )


def render_panel_2_compression() -> None:
    with st.container(border=True):
        step_header(2, "Compression")
        data = st.session_state.get("input_bytes", b"")
        if not data:
            st.info("Complete Input first."); return
        name = st.session_state.get("input_name", "payload.bin")

        left, right = st.columns([0.9,1.1], gap="large")
        with left:
            st.markdown("#### 🗂️ Configuration Setup")
            storage_mode = st.radio("Storage method", ["No compression", "Compression"], horizontal=True, key="storage_mode")

            # V4.3: selection is visible BEFORE running compression.
            selection_mode = "Balanced (recommended)"
            requested_method = None
            if storage_mode == "Compression":
                selection_mode = st.radio(
                    "Compressor selection",
                    ["Balanced (recommended)"],
                    horizontal=True,
                    key="compression_selection_mode",
                )
                expected_methods = _expected_compressor_methods(data, name)
                if selection_mode == "Choose compressor / version":
                    previous_method = st.session_state.get("compression_requested_method")
                    if previous_method not in expected_methods:
                        previous_method = expected_methods[0]
                    requested_method = st.selectbox(
                        "Compressor / version",
                        expected_methods,
                        index=expected_methods.index(previous_method),
                        format_func=_format_method_name,
                        key="compression_requested_method",
                    )

            if st.button("Run Compression", type="primary", width="stretch", key="run_compression"):
                try:
                    cands = compression_candidates(data, name)
                    analysis_rows = _compression_analysis_rows(data, cands)
                    st.session_state["compression_rows"] = cands
                    st.session_state["compression_analysis_rows"] = analysis_rows
                    if storage_mode == "No compression":
                        cand = next(c for c in cands if c["Method"] == "No compression")
                        st.session_state["compression_version_choice"] = "No compression"
                    elif selection_mode == "Choose compressor / version":
                        cand = _candidate_by_method(cands, str(requested_method))
                        st.session_state["compression_version_choice"] = cand["Method"]
                    elif selection_mode == "Smallest output":
                        cand = _compressed_candidate_choice(cands, "Auto best (smallest)")
                        st.session_state["compression_version_choice"] = "Smallest output"
                    else:
                        cand = _balanced_candidate_choice(cands)
                        st.session_state["compression_version_choice"] = "Balanced (recommended)"
                    _apply_selected_compression(cand)
                except Exception as exc:
                    st.error(f"Compression failed: {exc}")

            cands = st.session_state.get("compression_rows") or []
            analysis_rows = st.session_state.get("compression_analysis_rows") or []
            if cands and storage_mode == "No compression":
                try:
                    cand = next(c for c in cands if c["Method"] == "No compression")
                    _apply_selected_compression(cand)
                except Exception as exc:
                    st.error(f"Storage selection failed: {exc}")

            # After benchmarking, expose every ACTUALLY available candidate and
            # allow switching without recompressing. This is intentionally not
            # hidden behind Auto-best.
            if cands and storage_mode == "Compression":
                st.markdown("##### Available compressed versions")
                compressed = [c for c in cands if c.get("Method") != "No compression"]
                actual_methods = [c["Method"] for c in compressed]
                current_method = st.session_state.get("compression_method")
                if current_method not in actual_methods and actual_methods:
                    current_method = actual_methods[0]
                if actual_methods:
                    selected_method = st.selectbox(
                        "Use compressed version",
                        actual_methods,
                        index=actual_methods.index(current_method),
                        format_func=lambda m: f"{_format_method_name(m)} · {_candidate_by_method(cands,m).get('tier') or 'Exact'} · {len(_candidate_by_method(cands,m)['data']):,} B",
                        key="compression_actual_version_selector",
                    )
                    try:
                        _apply_selected_compression(_candidate_by_method(cands, selected_method))
                    except Exception as exc:
                        st.error(f"Compressed version selection failed: {exc}")

                st.markdown("##### Candidate benchmark")
                _display_compression_table(
                    [r for r in analysis_rows if r["Method"] != "No compression"],
                    st.session_state.get("compression_method", ""),
                )

        with right:
            stored = st.session_state.get("stored_bytes")
            if stored is None:
                st.info("Run Compression to create the stored representation.")
            else:
                st.markdown("#### Selected output")
                analysis_rows = st.session_state.get("compression_analysis_rows") or []
                selected_row = next((r for r in analysis_rows if r["Method"]==st.session_state.get("compression_method")), None)
                p1,p2,p3=st.columns(3)
                p1.metric("Method",_format_method_name(st.session_state.get("compression_method","—")))
                p2.metric("Stored data",f"{len(stored):,} B")
                p3.metric("Compression ratio",f"{(selected_row or {}).get('Compression ratio',1.0):.2f}×")
                _compression_quality_cards(selected_row,analysis_rows)
                property_table([
                    {"Property":"Stored file","Value":st.session_state.get("stored_name")},
                    {"Property":"Size reduction","Value":f"{(selected_row or {}).get('Size reduction (%)',0):.2f}%"},
                    {"Property":"Estimated full DNA length (HEDges)","Value":f"{(selected_row or {}).get('Estimated full DNA nt (HEDges)',math.ceil(len(stored)*8/103)*125):,} nt"},
                ])
                selected_cand=next((c for c in (st.session_state.get("compression_rows") or []) if c.get("Method")==st.session_state.get("compression_method")),None)
                preview_data,preview_name,preview_mime=_candidate_preview(selected_cand or {"data":stored,"name":st.session_state.get("stored_name"),"mime":st.session_state.get("stored_mime")})
                _preview_bytes(preview_data, preview_name, preview_mime, "Selected preview", key="selected_stored_text_preview")
                collect("Download stored representation", stored, file_name=st.session_state.get("stored_name","stored.bin"), mime=st.session_state.get("stored_mime","application/octet-stream"), width="stretch")

        analysis_rows = st.session_state.get("compression_analysis_rows") or []
        if analysis_rows:
            st.markdown("#### Compression Metrics")
            before, after = st.columns(2, gap="large")
            with before:
                st.markdown("##### Before compression")
                _preview_bytes(data, name, st.session_state.get("input_mime"), "Original preview", key="compression_before_text_preview")
                fmt,_,_= _detect_type(data,name)
                property_table([{"Property":"Data type","Value":fmt},{"Property":"Original size","Value":f"{len(data):,} B"}])
            with after:
                st.markdown("##### After compression")
                stored=st.session_state.get("stored_bytes")
                if stored is not None:
                    selected_cand=next((c for c in (st.session_state.get("compression_rows") or []) if c.get("Method")==st.session_state.get("compression_method")),None)
                    preview_data,preview_name,preview_mime=_candidate_preview(selected_cand or {"data":stored,"name":st.session_state.get("stored_name"),"mime":st.session_state.get("stored_mime")})
                    _preview_bytes(preview_data,preview_name,preview_mime,"Selected preview",key="compression_after_text_preview")
                    property_table([{"Property":"Representation","Value":_format_method_name(st.session_state.get("compression_method","—"))},{"Property":"Stored data","Value":f"{len(stored):,} B"}])
            _render_compression_charts(analysis_rows,st.session_state.get('compression_method',''))


# -----------------------------------------------------------------------------
# Panel 3 — DNA Encoding
# -----------------------------------------------------------------------------

def render_panel_3_encoding() -> None:
    with st.container(border=True):
        step_header(3, "Encoding")
        stored = st.session_state.get("stored_bytes")
        if stored is None:
            st.info("Complete Compression first."); return

        left, right = st.columns(2, gap="large")
        with left:
            st.markdown("#### 🧬 Design")
            design = "TVH-V5.2 HEDges"
            st.write("**HEDGES · TVH-V5.2**")
            if design != "TVH-V5.2 HEDges":
                st.session_state["sid_protection_mode"] = hb.SID_MODE_HAMMING
            sid_mode = st.radio(
                "Strand ID (SI) protection",
                [hb.SID_MODE_HAMMING, hb.SID_MODE_NO_HAMMING],
                horizontal=True,
                key="sid_protection_mode",
                disabled=(design != "TVH-V5.2 HEDges"),
            )
            if design != "TVH-V5.2 HEDges":
                sid_mode = hb.SID_MODE_HAMMING
            if design == "TVH-V5.2 HEDges":
                est = hb.estimate_strands(
                    stored,
                    st.session_state.get("stored_name", "stored.bin"),
                    st.session_state.get("stored_mime", "application/octet-stream"),
                )
                st.metric("Estimated HEDges strands", est)
                sid_stats = hb.sid_codebook_stats(max(1, min(est, hb.strand_capacity(sid_mode))), sid_mode)
                property_table([
                    {"Property":"Useful payload / strand","Value":"103 bits"},
                    {"Property":"BODY","Value":"76 nt"},
                    {"Property":"Full strand","Value":"125 nt"},
                    {"Property":"SID mode","Value":sid_mode},
                    {"Property":"SID Hamming constraint","Value":"dH ≥ 3" if sid_mode==hb.SID_MODE_HAMMING else "None"},
                    {"Property":"Observed min pairwise HD","Value":sid_stats.get("min_pairwise_hamming") if sid_stats.get("min_pairwise_hamming") is not None else "—"},
                    {"Property":"SID GC range","Value":f"{sid_stats.get('min_gc_percent',0):.1f}–{sid_stats.get('max_gc_percent',0):.1f}%"},
                    {"Property":"SID max homopolymer","Value":sid_stats.get("max_homopolymer")},
                    {"Property":"File framing","Value":"Embedded name, domain, byte length and CRC32"},
                    {"Property":"Tail padding","Value":"Removed from embedded byte length"},
                ])
                if est > hb.strand_capacity(sid_mode):
                    st.error(f"The selected representation requires {est} HEDges strands; this build supports {hb.strand_capacity(sid_mode)} in this SI mode. Compress the file further.")
            if st.button("Run DNA Encoding", type="primary", width="stretch", key="run_encoding"):
                _reset_after_encoding()
                try:
                    if design == "TVH-V5.2 HEDges":
                        strands, meta = hb.encode_bytes(
                            stored,
                            sid_mode=sid_mode,
                            filename=st.session_state.get("stored_name", "stored.bin"),
                            mime=st.session_state.get("stored_mime", "application/octet-stream"),
                        )
                        st.session_state.update({"hedges_strands":strands,"hedges_meta":meta,"legacy_dna":None,"legacy_meta":None,"legacy_bits":bytes_to_bits(stored)})
                    else:
                        dna, bits, meta = legacy_codec.encode_bytes(stored, design)
                        st.session_state.update({"legacy_dna":dna,"legacy_bits":bits,"legacy_meta":meta,"hedges_strands":None,"hedges_meta":None})
                    st.session_state.update({"dna_design":design,"encoded_ready":True})
                except Exception as exc:
                    st.error(f"Encoding failed: {exc}")

            st.markdown("##### 📄 Encoded data properties")
            property_table([
                {"Property":"Stored data","Value":f"{len(stored):,} B"},
                {"Property":"Input bits","Value":f"{len(stored)*8:,} bits"},
                {"Property":"Selected DNA design","Value":DESIGN_DISPLAY.get(design,design)},
            ])

        with right:
            if not st.session_state.get("encoded_ready"):
                st.info("Run DNA Encoding to inspect the output.")
            else:
                design = st.session_state.get("dna_design")
                st.markdown("#### 🧬 DNA output")
                if design == "TVH-V5.2 HEDges":
                    strands: List[hb.HedgesStrand] = st.session_state.get("hedges_strands", [])
                    meta = st.session_state.get("hedges_meta", {})
                    if strands:
                        c1,c2,c3,c4 = st.columns(4)
                        c1.metric("Strands", len(strands)); c2.metric("BODY", "76 nt"); c3.metric("Full strand", "125 nt"); c4.metric("Rescue", meta.get("rescue_strands",0))
                        property_table([
                            {"Property":"DNA design","Value":"TVH-V5.2 HEDges"},
                            {"Property":"Payload / strand","Value":"103 bits"},
                            {"Property":"Internal checkpoints","Value":"4 bits (2 + 2)"},
                            {"Property":"Terminal verifier","Value":"8 nt"},
                            {"Property":"SID mode","Value":meta.get("sid_mode", hb.SID_MODE_HAMMING)},
                            {"Property":"SID Hamming constraint","Value":"dH ≥ 3" if meta.get("sid_mode",hb.SID_MODE_HAMMING)==hb.SID_MODE_HAMMING else "None"},
                            {"Property":"Observed min pairwise HD","Value":(meta.get("sid_codebook") or {}).get("min_pairwise_hamming", "—")},
                            {"Property":"SID GC range","Value":f"{float((meta.get('sid_codebook') or {}).get('min_gc_percent',0)):.1f}–{float((meta.get('sid_codebook') or {}).get('max_gc_percent',0)):.1f}%"},
                            {"Property":"SID max homopolymer","Value":(meta.get("sid_codebook") or {}).get("max_homopolymer", "—")},
                            {"Property":"Padding bits","Value":meta.get("padding_bits")},
                            {"Property":"Embedded file name","Value":meta.get("stored_filename") or "—"},
                            {"Property":"Framing overhead","Value":f"{int(meta.get('envelope_overhead_bytes',0))} B"},
                            {"Property":"Nominal BODY density","Value":f"{103/76:.4f} bit/nt"},
                        ])
                        bio_rows=[{
                            "Strand ID":x.strand_id,"Seed":x.seed_id,"Status":x.encode_status,"Certified":x.certified,"Rescue":x.rescue,
                            "BODY GC":f"{x.body_gc:.4f}","BODY longest HP":x.body_hp,"Full GC":f"{x.full_gc:.4f}","Full longest HP":x.full_hp,
                        } for x in strands]
                        st.dataframe(pd.DataFrame(bio_rows),width="stretch",hide_index=True,height=min(320,38*len(bio_rows)+38))
                        inspect=int(st.number_input("Inspect encoded strand",min_value=0,max_value=max(0,len(strands)-1),value=0,step=1,key="encoding_inspect_hedges"))
                        show_segmented(strands[inspect].row())
                        d1,d2=st.columns(2)
                        collect("Download HEDges BODY DNA",("\n".join(x.body for x in strands)+"\n").encode(),"hedges_body_dna.txt","text/plain",width="stretch")
                        collect("Download 125-nt strands",("\n".join(x.full_strand for x in strands)+"\n").encode(),"hedges_full_strands.txt","text/plain",width="stretch")
                else:
                    dna = st.session_state.get("legacy_dna", "")
                    c1,c2,c3,c4 = st.columns(4)
                    c1.metric("DNA length", f"{len(dna):,} nt")
                    c2.metric("GC content", f"{legacy_codec.gc_content(dna):.4f}")
                    c3.metric("Longest homopolymer", legacy_codec.longest_homopolymer(dna))
                    c4.metric("Information density", f"{len(stored)*8/max(1,len(dna)):.4f} bit/nt")
                    seq_box(dna)
                    collect("Download encoded DNA",dna.encode(),"encoded_dna.txt","text/plain",width="stretch")

        # The accepted App1 randomness analysis is retained and expanded to per-strand evaluation.
        if st.session_state.get("encoded_ready"):
            design=st.session_state.get("dna_design")
            if design=="TVH-V5.2 HEDges":
                strands=st.session_state.get("hedges_strands") or []
                overall="".join(x.body for x in strands)
                strand_sequences=[(int(x.strand_id),x.body,103) for x in strands]
                randomness_ui.render_randomness_analysis(
                    overall,strand_sequences,state_prefix="encoding",overall_bit_count=len(stored)*8,sequence_label="HEDges BODY (76 nt)"
                )
            else:
                dna=st.session_state.get("legacy_dna","")
                virtual_rows=build_legacy_strands(dna)
                strand_sequences=[(int(r["Strand ID"]),str(r.get("Payload","")),None) for r in virtual_rows if r.get("Payload")]
                randomness_ui.render_randomness_analysis(
                    dna,strand_sequences,state_prefix="encoding",overall_bit_count=len(stored)*8,sequence_label="encoded payload segment (primers/SI/filler excluded)"
                )


# -----------------------------------------------------------------------------
# Panel 4 — Strand Design + sequencing / read-level ECC
# -----------------------------------------------------------------------------

def _ensure_strands() -> List[dict]:
    design = st.session_state.get("dna_design")
    if not st.session_state.get("encoded_ready"):
        return []
    if design == "TVH-V5.2 HEDges":
        sid_mode = _current_hedges_sid_mode()
        rows = [{**s.row(), "SID mode": sid_mode} for s in (st.session_state.get("hedges_strands") or [])]
    else:
        rows = build_legacy_strands(st.session_state.get("legacy_dna", ""))
    st.session_state["strand_rows"] = rows
    return rows



def _direct_mutate(seq: str, rng: random.Random, sub_rate: float, ins_rate: float, del_rate: float) -> tuple[str, dict, List[dict]]:
    """Apply exact S/I/D events and retain a trace only for visualization/reporting."""
    seq = fastq_ecc.clean_dna(seq)
    out: List[str] = []
    counts = {"sub": 0, "ins": 0, "del": 0}
    events: List[dict] = []
    for clean_index, base in enumerate(seq):
        if rng.random() < del_rate:
            counts["del"] += 1
            events.append({"kind":"D","clean_index":clean_index,"original":base,"observed":None})
            continue
        observed = base
        if rng.random() < sub_rate:
            observed = rng.choice([b for b in "ACGT" if b != base])
            counts["sub"] += 1
            events.append({"kind":"S","clean_index":clean_index,"original":base,"observed":observed})
        out.append(observed)
        if rng.random() < ins_rate:
            inserted = rng.choice("ACGT")
            out.append(inserted)
            counts["ins"] += 1
            events.append({"kind":"I","clean_index":clean_index,"original":None,"observed":inserted})
    if seq and rng.random() < ins_rate:
        inserted = rng.choice("ACGT")
        out.append(inserted)
        counts["ins"] += 1
        events.append({"kind":"I","clean_index":len(seq)-1,"original":None,"observed":inserted,"terminal":True})
    return "".join(out), counts, events


def build_direct_error_rows(rows: List[dict], sub_rate: float, ins_rate: float, del_rate: float, seed: int) -> tuple[List[dict], dict]:
    """Mutate BODY/payload while FBR/SI/RBR remain structural and unchanged."""
    noisy_rows: List[dict] = []
    total = {"sub": 0, "ins": 0, "del": 0, "strands": 0}
    for row in rows:
        sid = int(row["Strand ID"])
        clean_body = fastq_ecc.clean_dna(str(row.get("Payload", "")) + str(row.get("Filler", "")))
        rng = random.Random(int(seed) + sid * 1_000_003)
        observed_body, counts, events = _direct_mutate(clean_body, rng, float(sub_rate), float(ins_rate), float(del_rate))
        full = str(row.get("FBR", FBR)) + str(row.get("SI", "")) + observed_body + str(row.get("RBR", RBR))
        noisy_rows.append({
            **row,
            "Clean BODY": clean_body,
            "Received BODY": observed_body,
            "Clean BODY length": len(clean_body),
            "Received BODY length": len(observed_body),
            "Received full strand": full,
            "Substitutions": counts["sub"],
            "Insertions": counts["ins"],
            "Deletions": counts["del"],
            "_events": events,
        })
        total["sub"] += counts["sub"]; total["ins"] += counts["ins"]; total["del"] += counts["del"]; total["strands"] += 1
    total["errors"] = total["sub"] + total["ins"] + total["del"]
    return noisy_rows, total


def _error_legend() -> None:
    st.markdown(
        '<span class="err-legend err-sub">Substitution A→G</span>'
        '<span class="err-legend err-ins">Insertion +C</span>'
        '<span class="err-legend err-del">Deletion −T</span>',
        unsafe_allow_html=True,
    )


def _direct_trace_html(clean: str, events: List[dict]) -> str:
    subs={int(e["clean_index"]):e for e in events if e.get("kind")=="S"}
    dels={int(e["clean_index"]):e for e in events if e.get("kind")=="D"}
    ins={}
    for e in events:
        if e.get("kind")=="I": ins.setdefault(int(e["clean_index"]),[]).append(e)
    out=[]
    for i,b in enumerate(clean):
        if i in dels:
            out.append(f'<span class="err-del" title="Deletion at BODY {i}">−{html.escape(b)}</span>')
        elif i in subs:
            ob=str(subs[i].get("observed") or "")
            out.append(f'<span class="err-sub" title="Substitution at BODY {i}">{html.escape(b)}→{html.escape(ob)}</span>')
        else:
            out.append(html.escape(b))
        for e in ins.get(i,[]):
            out.append(f'<span class="err-ins" title="Insertion after BODY {i}">+{html.escape(str(e.get("observed") or ""))}</span>')
    return ''.join(out)


def _render_direct_trace(row: dict) -> None:
    _error_legend()
    body=_direct_trace_html(str(row.get("Clean BODY","")), list(row.get("_events") or []))
    st.markdown(
        '<div class="seq-box">'
        f'<span class="segment seg-fbr">FBR {html.escape(str(row.get("FBR","")))}</span>'
        f'<span class="segment seg-si">SI {html.escape(str(row.get("SI","")))}</span>'
        f'<span class="segment seg-body">BODY {body}</span>'
        f'<span class="segment seg-rbr">RBR {html.escape(str(row.get("RBR","")))}</span>'
        '</div>', unsafe_allow_html=True,
    )


def _alignment_diff_html(reference: str, observed: str) -> str:
    """Display-only alignment; never used by consensus or decoder."""
    reference=fastq_ecc.clean_dna(reference); observed=fastq_ecc.clean_dna(observed)
    aa,bb,_=fastq_ecc._nw_align(reference,observed)
    out=[]
    for a,b in zip(aa,bb):
        if a==b:
            out.append(html.escape(a))
        elif a=='-':
            out.append(f'<span class="err-ins">+{html.escape(b)}</span>')
        elif b=='-':
            out.append(f'<span class="err-del">−{html.escape(a)}</span>')
        else:
            out.append(f'<span class="err-sub">{html.escape(a)}→{html.escape(b)}</span>')
    return ''.join(out)


def _decode_base4_si(si: str) -> int:
    n=0
    for c in fastq_ecc.clean_dna(si):
        n=n*4+"ACGT".index(c)
    return n


def _current_hedges_sid_mode() -> str:
    return hb.normalize_sid_mode((st.session_state.get("hedges_meta") or {}).get("sid_mode", hb.SID_MODE_HAMMING))


def _fastq_record_target(record: dict, design: str) -> tuple[int|None, dict, str]:
    """Resolve SID from sequence-intrinsic SI first; header SID is fallback only."""
    seq=fastq_ecc.clean_dna(record.get("sequence","")); q=list(record.get("phred") or [])
    # Full prepared strand: structural FBR/SI/RBR are unchanged by this simulator.
    if seq.startswith(FBR) and seq.endswith(RBR) and len(seq) >= len(FBR)+LEGACY_SI_NT+len(RBR):
        si=seq[len(FBR):len(FBR)+LEGACY_SI_NT]
        try:
            sid=hb.decode_si(si, _current_hedges_sid_mode())[0] if design=="TVH-V5.2 HEDges" else _decode_base4_si(si)
        except Exception:
            sid=None
        if sid is not None:
            if design=="TVH-V5.2 HEDges":
                a=len(FBR)+LEGACY_SI_NT; b=len(seq)-len(RBR)
                target=seq[a:b]; tq=q[a:b]
                return int(sid), {**record,"sequence":target,"phred":tq}, "Embedded SI"
            return int(sid), record, "Embedded SI"
    # External BODY-only FASTQ may explicitly carry a strand ID in its own header.
    sid=fastq_ecc.sid_from_header(record.get("header",""))
    if sid is not None:
        return int(sid), record, "FASTQ header"
    return None, record, "Unassigned"


def independent_fastq_consensus(data: bytes|str, design: str, method: str, sid_mode: str = hb.SID_MODE_AUTO) -> tuple[Dict[int,str],List[dict],List[dict],dict]:
    """FASTQ-only grouping/voting; no clean strand, payload, or record-order hint."""
    if design == "TVH-V5.2 HEDges":
        mode = hb.normalize_sid_mode(sid_mode)
        if mode == hb.SID_MODE_AUTO:
            decoder = lambda si: hb.decode_si_auto(si)[0]
        else:
            decoder = lambda si: hb.decode_si(si, mode)[0]
        target = "body"
    else:
        mode = None
        decoder = _decode_base4_si
        target = "full"
    cmap, rows, positions, stats = fastq_ecc.consensus_from_intrinsic_ids(
        data, fbr=FBR, rbr=RBR, si_nt=LEGACY_SI_NT,
        si_decoder=decoder, method=method, target=target,
    )
    if design == "TVH-V5.2 HEDges":
        stats = dict(stats)
        stats["sid_mode"] = mode
    return cmap, rows, positions, stats


def _parse_fasta_or_dna_text(data: bytes | str) -> tuple[Dict[int, str], str]:
    """Parse FASTA or one-sequence-per-line DNA. Headers may carry strand IDs."""
    if isinstance(data, (bytes, bytearray)):
        text = bytes(data).decode("utf-8-sig", errors="ignore")
    else:
        text = str(data)
    raw_lines = [x.strip() for x in text.splitlines() if x.strip()]
    if not raw_lines:
        raise ValueError("No DNA sequences found")

    seqs: Dict[int, str] = {}
    if any(x.startswith(">") for x in raw_lines):
        header = None
        parts: List[str] = []
        records: List[tuple[str, str]] = []
        for line in raw_lines:
            if line.startswith(">"):
                if header is not None:
                    records.append((header, fastq_ecc.clean_dna("".join(parts))))
                header = line[1:].strip()
                parts = []
            else:
                parts.append(line)
        if header is not None:
            records.append((header, fastq_ecc.clean_dna("".join(parts))))
        for i, (h, seq) in enumerate(records):
            if not seq:
                continue
            sid = fastq_ecc.sid_from_header(h)
            sid = i if sid is None else int(sid)
            if sid in seqs:
                raise ValueError(f"Duplicate strand ID {sid}")
            seqs[sid] = seq
        source_format = "FASTA"
    else:
        dna_lines = []
        for line in raw_lines:
            seq = fastq_ecc.clean_dna(line)
            if seq:
                dna_lines.append(seq)
        if not dna_lines:
            raise ValueError("No DNA sequences found")
        seqs = {i: seq for i, seq in enumerate(dna_lines)}
        source_format = "DNA TXT"
    return seqs, source_format


def _hedges_bodies_from_sequence_map(seqs: Dict[int, str], sid_mode: str = hb.SID_MODE_AUTO) -> tuple[Dict[int, str], List[dict], str]:
    """Accept HEDges BODY consensus sequences or full FBR|SI|BODY|RBR strands."""
    if not seqs:
        raise ValueError("No sequences available")
    full_like = all(seq.startswith(FBR) and seq.endswith(RBR) for seq in seqs.values())
    bodies: Dict[int, str] = {}
    parsed: List[dict] = []
    if full_like:
        # Recover SID from the visible SI; line/header order is not used as a hidden decode hint.
        for source_sid, full in sorted(seqs.items()):
            central = full[len(FBR):len(full)-len(RBR)]
            if len(central) <= hb.SI_NT:
                raise ValueError(f"Sequence {source_sid} is too short for SI + BODY")
            obs_si = central[:hb.SI_NT]
            mode = hb.normalize_sid_mode(sid_mode)
            if mode == hb.SID_MODE_AUTO:
                sid, canonical, si_d, detected_sid_mode = hb.decode_si_auto(obs_si)
            else:
                sid, canonical, si_d = hb.decode_si(obs_si, mode)
                detected_sid_mode = mode
            if sid in bodies:
                raise ValueError(f"Duplicate HEDges strand ID {sid}")
            body = central[hb.SI_NT:]
            bodies[sid] = body
            parsed.append({"Input ID": source_sid, "Strand ID": sid, "Observed SI": obs_si, "SI distance": si_d, "SID mode": detected_sid_mode, "BODY nt": len(body)})
        kind = "Full strands"
    else:
        for sid, seq in sorted(seqs.items()):
            bodies[int(sid)] = fastq_ecc.clean_dna(seq)
            parsed.append({"Input ID": sid, "Strand ID": int(sid), "BODY nt": len(bodies[int(sid)])})
        kind = "Consensus / BODY sequences"
    if sorted(bodies) != list(range(len(bodies))):
        raise ValueError("HEDges strand IDs must be contiguous from 0")
    return bodies, parsed, kind


def _looks_like_fastq(data: bytes, filename: str = "") -> bool:
    raw = bytes(data or b"")
    if raw[:2] == b"\x1f\x8b":
        try:
            raw = gzip.decompress(raw)
        except Exception:
            return False
    ext = Path(filename).suffix.lower()
    if ext in {".fastq", ".fq"}:
        return True
    text = raw.decode("utf-8-sig", errors="ignore")
    lines = [x for x in text.splitlines() if x.strip()]
    return len(lines) >= 4 and lines[0].startswith("@") and lines[2].startswith("+")


def _full_strand_map_from_direct_rows(noisy_rows: List[dict]) -> Dict[int, str]:
    return {int(r["Strand ID"]): fastq_ecc.clean_dna(r.get("Received full strand", "")) for r in noisy_rows}

def _levenshtein(a: str, b: str) -> int:
    a=fastq_ecc.clean_dna(a); b=fastq_ecc.clean_dna(b)
    prev=list(range(len(b)+1))
    for i,ca in enumerate(a,1):
        cur=[i]
        for j,cb in enumerate(b,1):
            cur.append(min(cur[-1]+1,prev[j]+1,prev[j-1]+(ca!=cb)))
        prev=cur
    return prev[-1]


def _consensus_reference_map(rows: List[dict]) -> Dict[int,str]:
    design=st.session_state.get("dna_design")
    if design=="TVH-V5.2 HEDges":
        return {int(r["Strand ID"]):fastq_ecc.clean_dna(r.get("Payload","")) for r in rows}
    return {int(r["Strand ID"]):fastq_ecc.clean_dna(r.get("Full strand","")) for r in rows}


def _entropy_from_counts(counts: List[int]) -> float:
    total=sum(counts)
    if total<=0: return 0.0
    out=0.0
    for c in counts:
        if c:
            p=c/total; out-=p*math.log2(p)
    return out


def _consensus_evaluation(rows: List[dict], cmap: Dict[int,str], position_rows: List[dict], method: str) -> tuple[List[dict],List[dict],dict]:
    """Post-hoc evaluation only. The clean sequence is never passed into voting."""
    refs=_consensus_reference_map(rows)
    bypos: Dict[int,List[dict]]={}
    for r in position_rows:
        bypos.setdefault(int(r.get("Strand ID",-1)),[]).append(r)
    strand_rows=[]; evaluated_positions=[]
    total_edits=0; total_den=0; probs=[]; entropies=[]; total_ties=0; total_low=0
    for sid in sorted(refs):
        ref=refs[sid]; seq=fastq_ecc.clean_dna(cmap.get(sid,"")); p=bypos.get(sid,[])
        edit=_levenshtein(ref,seq) if seq else len(ref)
        denom=max(1,len(ref),len(seq)); acc=1-edit/denom
        conf=[float(x.get("Confidence",0) or 0) for x in p]
        ents=[_entropy_from_counts([int(x.get(b,0) or 0) for b in "ACGT"]) for x in p]
        ties=sum(bool(x.get("Tie")) for x in p); low=sum(float(x.get("Confidence",0) or 0)<.70 for x in p)
        probs.extend(conf); entropies.extend(ents); total_ties+=ties; total_low+=low; total_edits+=edit; total_den+=denom
        strand_rows.append({
            "Strand ID":sid,"Coverage":next((int(x.get("Coverage",0) or 0) for x in (st.session_state.get("consensus_rows") or []) if int(x.get("Strand ID",-1))==sid),0),
            "Reference length":len(ref),"Consensus length":len(seq),"Consensus accuracy":acc,
            "Mean consensus probability":sum(conf)/max(1,len(conf)),"Minimum consensus probability":min(conf) if conf else 0.0,
            "Mean entropy (bits)":sum(ents)/max(1,len(ents)),"Tie positions":ties,"Low-confidence positions":low,
            "Consensus edit distance":edit,"Consensus sequence":seq,"Voting method":method,
        })
        for x in p:
            pos=int(x.get("Position",0) or 0)
            refnt=ref[pos-1] if 1<=pos<=len(ref) else "—"
            evaluated_positions.append({**x,"Voting entropy (bits)":_entropy_from_counts([int(x.get(b,0) or 0) for b in "ACGT"]),
                                        "Reference nt":refnt,"Consensus correct":"Yes" if refnt in "ACGT" and str(x.get("Consensus"))==refnt else "No",
                                        "Voting method":method})
    stats={
        "voting_method":method,"strand_count":len(refs),"covered_strands":len(cmap),
        "consensus_accuracy":1-total_edits/max(1,total_den),"consensus_edit_distance":total_edits,
        "mean_consensus_probability":sum(probs)/max(1,len(probs)),"minimum_consensus_probability":min(probs) if probs else 0.0,
        "mean_entropy_bits":sum(entropies)/max(1,len(entropies)),"tie_positions":total_ties,"low_confidence_positions":total_low,
    }
    return strand_rows,evaluated_positions,stats


def _render_consensus_analysis(rows: List[dict]) -> None:
    cmap=st.session_state.get("consensus_map") or {}
    if not cmap: return
    srows=st.session_state.get("consensus_eval_rows") or []
    prows=st.session_state.get("consensus_eval_positions") or []
    stats=st.session_state.get("consensus_eval_stats") or {}
    if not stats: return
    st.markdown("#### 🗳️ Consensus Reconstruction Analysis")
    v1,v2,v3,v4=st.columns(4)
    v1.metric("Consensus accuracy",f"{100*float(stats.get('consensus_accuracy',0)):.4f}%")
    v2.metric("Mean confidence",f"{100*float(stats.get('mean_consensus_probability',0)):.2f}%")
    v3.metric("Low-confidence positions",f"{int(stats.get('low_confidence_positions',0)):,}")
    v4.metric("Final DNA edit distance",f"{int(stats.get('consensus_edit_distance',0)):,}")
    property_table([
        {"Property":"Consensus method","Value":stats.get("voting_method","—")},
        {"Property":"Covered strands","Value":f"{stats.get('covered_strands',0)}/{stats.get('strand_count',0)}"},
        {"Property":"Minimum confidence","Value":f"{100*float(stats.get('minimum_consensus_probability',0)):.2f}%"},
        {"Property":"Mean voting entropy","Value":f"{float(stats.get('mean_entropy_bits',0)):.4f} bits"},
        {"Property":"Tie positions","Value":int(stats.get("tie_positions",0))},
    ])
    st.markdown("##### Consensus result by strand")
    shown=[]
    for r in srows:
        shown.append({
            "Strand ID":r["Strand ID"],"Coverage":r["Coverage"],"Reference length":r["Reference length"],"Consensus length":r["Consensus length"],
            "Accuracy":f"{100*float(r['Consensus accuracy']):.4f}%","Mean confidence":f"{100*float(r['Mean consensus probability']):.2f}%",
            "Minimum confidence":f"{100*float(r['Minimum consensus probability']):.2f}%","Mean entropy":f"{float(r['Mean entropy (bits)']):.4f}",
            "Ties":r["Tie positions"],"Low-confidence":r["Low-confidence positions"],"Edit distance":r["Consensus edit distance"],
        })
    st.dataframe(pd.DataFrame(shown),width="stretch",hide_index=True,height=min(350,38*len(shown)+38))
    if prows:
        maxsid=max(int(r.get("Strand ID",0)) for r in prows)
        sid=int(st.number_input("Inspect position-level evidence for strand",min_value=0,max_value=maxsid,value=0,step=1,key="inspect_consensus_evidence"))
        p=[r for r in prows if int(r.get("Strand ID",-1))==sid]
        table=[]
        for r in p:
            table.append({
                "Position":r.get("Position"),"A/C/G/T votes":f"{r.get('A',0)}/{r.get('C',0)}/{r.get('G',0)}/{r.get('T',0)}",
                "Gap votes":r.get("Gap votes",0),"Consensus":r.get("Consensus"),"Confidence":f"{100*float(r.get('Confidence',0)):.2f}%",
                "Entropy":f"{float(r.get('Voting entropy (bits)',0)):.4f}","Tie":"Yes" if r.get("Tie") else "No",
                "Reference nt":r.get("Reference nt"),"Consensus correct":r.get("Consensus correct"),
            })
        st.dataframe(pd.DataFrame(table),width="stretch",hide_index=True,height=320)
    d1,d2=st.columns(2)
    collect("Download consensus strand report",_df_csv(srows),"consensus_strand_report.csv","text/csv",width="stretch")
    collect("Download position-level evidence",_df_csv(prows),"consensus_position_evidence.csv","text/csv",width="stretch")


def render_panel_4_experiment() -> None:
    with st.container(border=True):
        step_header(4, "Strand Preparation")
        if not st.session_state.get("encoded_ready"):
            st.info("Complete DNA Encoding first."); return
        rows = st.session_state.get("strand_rows") or []
        left, right = st.columns(2, gap="large")

        with left:
            st.markdown("#### 🧵 Strand Preparation")
            prep_props=[
                {"Property":"Target total strand length","Value":"125 nt"},
                {"Property":"FBR","Value":FBR},
                {"Property":"SI","Value":"8 nt"},
                {"Property":"Payload / BODY","Value":"76 nt capacity"},
                {"Property":"RBR","Value":RBR},
            ]
            if st.session_state.get("dna_design")=="TVH-V5.2 HEDges":
                prep_props.insert(3,{"Property":"SID mode","Value":_current_hedges_sid_mode()})
            property_table(prep_props)
            if st.button("Run Strand Preparation", type="primary", width="stretch", key="prepare_strands"):
                try:
                    rows = _ensure_strands()
                    for k in ["direct_error_rows","direct_error_stats","fastq_bytes","fastq_sim_stats","fastq_origin","consensus_map","consensus_rows","consensus_position_rows","consensus_stats","consensus_eval_rows","consensus_eval_positions","consensus_eval_stats","decode_result"]:
                        st.session_state.pop(k,None)
                except Exception as exc:
                    st.error(f"Strand preparation failed: {exc}")
            rows = st.session_state.get("strand_rows") or []
            if rows:
                lengths=[len(fastq_ecc.clean_dna(r.get("Full strand",""))) for r in rows]
                gcs=[float(r.get("GC content",0) or 0) for r in rows]
                hps=[int(r.get("Longest homopolymer",0) or 0) for r in rows]
                c1,c2,c3,c4=st.columns(4)
                c1.metric("Strands",len(rows)); c2.metric("Total strand nt",f"{sum(lengths):,}"); c3.metric("Average GC",f"{np.mean(gcs):.4f}"); c4.metric("Max strand HP",max(hps) if hps else 0)
                table=[]
                for r in rows[:500]:
                    table.append({
                        "Strand ID":r.get("Strand ID"),"SI":r.get("SI"),"Payload length":r.get("Payload length"),
                        "Filler length":len(str(r.get("Filler",""))),"Total length":len(fastq_ecc.clean_dna(r.get("Full strand",""))),
                        "GC content":f"{float(r.get('GC content',0)):.4f}","Longest homopolymer":r.get("Longest homopolymer"),
                        "HEDges status":r.get("HEDges status","—"),
                    })
                st.markdown("##### 🧾 Strand table")
                st.dataframe(pd.DataFrame(table),width="stretch",hide_index=True,height=min(310,38*len(table)+38))
                sid=int(st.number_input("Inspect clean strand",min_value=0,max_value=max(0,len(rows)-1),value=0,step=1,key="inspect_sid"))
                show_segmented(rows[sid])
                d1,d2=st.columns(2)
                collect("Download prepared strands CSV",_df_csv(rows),"prepared_strands.csv","text/csv",width="stretch")
                collect("Download full strands TXT",("\n".join(r["Full strand"] for r in rows)+"\n").encode(),"full_strands.txt","text/plain",width="stretch")

        with right:
            st.markdown("#### 🧪 Sequencing / Error Model")
            rows=st.session_state.get("strand_rows") or []
            if not rows:
                st.info("Run Strand Preparation first.")
            else:
                error_model=st.radio(
                    "Model",
                    ["Direct sequence errors","FASTQ reads"],
                    horizontal=True,
                    key="sequencing_error_model",
                )

                if error_model=="Direct sequence errors":
                    e1,e2,e3=st.columns(3)
                    sub=e1.number_input("Substitution rate",0.0,0.10,0.003,0.001,format="%.4f",key="direct_sub")
                    ins=e2.number_input("Insertion rate",0.0,0.05,0.0002,0.0001,format="%.4f",key="direct_ins")
                    dele=e3.number_input("Deletion rate",0.0,0.05,0.0002,0.0001,format="%.4f",key="direct_del")
                    seed=st.number_input("Random seed",0,value=7,step=1,key="direct_seed")
                    if st.button("Add errors to DNA sequences",type="primary",width="stretch",key="run_direct_errors"):
                        noisy_rows,stats=build_direct_error_rows(rows,float(sub),float(ins),float(dele),int(seed))
                        st.session_state["direct_error_rows"]=noisy_rows
                        st.session_state["direct_error_stats"]=stats
                        st.session_state.pop("decode_result",None)
                    noisy_rows=st.session_state.get("direct_error_rows") or []
                    stats=st.session_state.get("direct_error_stats") or {}
                    if noisy_rows:
                        m1,m2,m3,m4=st.columns(4)
                        m1.metric("Substitutions",stats.get("sub",0)); m2.metric("Insertions",stats.get("ins",0)); m3.metric("Deletions",stats.get("del",0)); m4.metric("Total errors",stats.get("errors",0))
                        shown=[{
                            "Strand ID":r["Strand ID"],
                            "Clean BODY nt":r["Clean BODY length"],
                            "Received BODY nt":r["Received BODY length"],
                            "Substitutions":r["Substitutions"],
                            "Insertions":r["Insertions"],
                            "Deletions":r["Deletions"],
                        } for r in noisy_rows]
                        st.dataframe(pd.DataFrame(shown),width="stretch",hide_index=True,height=min(300,38*len(shown)+38))
                        inspect=int(st.number_input("Inspect noisy strand",min_value=0,max_value=max(0,len(noisy_rows)-1),value=0,step=1,key="inspect_direct_sid"))
                        _render_direct_trace(noisy_rows[inspect])
                        d1,d2=st.columns(2)
                        collect("Download noisy strands TXT",("\n".join(r["Received full strand"] for r in noisy_rows)+"\n").encode(),"direct_error_full_strands.txt","text/plain",width="stretch")
                        collect("Download error report CSV",_df_csv(shown),"direct_error_report.csv","text/csv",width="stretch")

                else:
                    source=st.radio("FASTQ source",["Simulate from prepared strands","Upload FASTQ"],horizontal=True,key="fastq_source")
                    design=st.session_state.get("dna_design")
                    if source=="Simulate from prepared strands":
                        q1,q2,q3=st.columns(3)
                        coverage=q1.number_input("Coverage / noisy read sets",1,200,10,1,key="fq_cov")
                        error_model_display=q2.selectbox("Error probability model",["Error rate","Q-score"],key="fq_error_model")
                        seed=q3.number_input("Random seed",0,value=11,step=1,key="fq_seed")
                        error_model="Independent Bernoulli per nucleotide" if error_model_display=="Error rate" else "Phred Q-score driven"

                        qscore_profile="Fixed Q"; q_mean=30.0; q_std=4.0; q_start=35.0; q_end=20.0; q_fixed=30.0
                        if error_model_display=="Q-score":
                            qscore_profile=st.selectbox("Q-score profile",["Normal Q distribution","Linear quality decay","Fixed Q"],key="fq_q_profile")
                            if qscore_profile=="Normal Q distribution":
                                a,b=st.columns(2)
                                q_mean=a.number_input("Mean Q",2.0,40.0,30.0,1.0,key="fq_q_mean")
                                q_std=b.number_input("Q standard deviation",0.0,15.0,4.0,0.5,key="fq_q_std")
                            elif qscore_profile=="Linear quality decay":
                                a,b=st.columns(2)
                                q_start=a.number_input("Start Q",2.0,40.0,35.0,1.0,key="fq_q_start")
                                q_end=b.number_input("End Q",2.0,40.0,20.0,1.0,key="fq_q_end")
                            else:
                                q_fixed=st.number_input("Fixed Q",2.0,40.0,30.0,1.0,key="fq_q_fixed")

                        e1,e2,e3,e4=st.columns(4)
                        sub=e1.number_input("Substitution rate",0.0,0.50,0.0200,0.001,format="%.4f",key="fq_sub",disabled=(error_model_display=="Q-score"))
                        ins=e2.number_input("Insertion rate",0.0,0.20,0.0000,0.0005,format="%.4f",key="fq_ins")
                        dele=e3.number_input("Deletion rate",0.0,0.20,0.0000,0.0005,format="%.4f",key="fq_del")
                        dropout=e4.number_input("Strand dropout rate",0.0,1.0,0.0000,0.001,format="%.4f",key="fq_dropout")
                        if st.button("Simulate FASTQ",type="primary",width="stretch",key="simulate_fastq"):
                            fq,stats=fastq_ecc.simulate_fastq(
                                rows, coverage=int(coverage), sub_rate=float(sub), ins_rate=float(ins), del_rate=float(dele),
                                seed=int(seed), error_probability_model=error_model, qscore_profile=qscore_profile,
                                q_mean=float(q_mean), q_std=float(q_std), q_start=float(q_start), q_end=float(q_end), q_fixed=float(q_fixed),
                                strand_dropout_rate=float(dropout),
                            )
                            st.session_state["fastq_bytes"]=fq; st.session_state["fastq_sim_stats"]=stats; st.session_state["fastq_origin"]="simulated"
                            for k in ["consensus_map","consensus_rows","consensus_position_rows","consensus_stats","consensus_eval_rows","consensus_eval_positions","consensus_eval_stats","decode_result"]:
                                st.session_state.pop(k,None)
                    else:
                        up=st.file_uploader("Upload FASTQ / FASTQ.gz",type=["fastq","fq","txt","gz"],key="fastq_upload")
                        if up is not None:
                            st.session_state["fastq_bytes"]=up.getvalue(); st.session_state["fastq_sim_stats"]=None; st.session_state["fastq_origin"]="uploaded"
                            for k in ["consensus_map","consensus_rows","consensus_position_rows","consensus_stats","consensus_eval_rows","consensus_eval_positions","consensus_eval_stats","decode_result"]:
                                st.session_state.pop(k,None)

                    fq=st.session_state.get("fastq_bytes")
                    parsed=[]
                    if fq:
                        try:
                            parsed=fastq_ecc.parse_fastq(fq)
                            lens=[len(r["sequence"]) for r in parsed]; qvals=[q for r in parsed for q in r["phred"]]
                            f1,f2,f3,f4=st.columns(4)
                            f1.metric("Total reads",f"{len(parsed):,}"); f2.metric("Mean read length",f"{np.mean(lens):.1f} nt" if lens else "—")
                            f3.metric("Mean Q",f"{np.mean(qvals):.2f}" if qvals else "—"); f4.metric("Read length range",f"{min(lens)}–{max(lens)} nt" if lens else "—")
                        except Exception as exc:
                            st.error(f"FASTQ parse failed: {exc}")
                        stats=st.session_state.get("fastq_sim_stats")
                        if stats:
                            ss=stats.get("summary",{})
                            s1,s2,s3,s4=st.columns(4)
                            s1.metric("Substitutions",ss.get("sub",0)); s2.metric("Insertions",ss.get("ins",0)); s3.metric("Deletions",ss.get("del",0)); s4.metric("Dropped strands",ss.get("dropped",0))
                            if stats.get("per_strand"):
                                with st.expander("Injected errors by strand",expanded=False):
                                    st.dataframe(pd.DataFrame(stats["per_strand"]),width="stretch",hide_index=True,height=260)
                            property_table([
                                {"Property":"Error probability model","Value":"Q-score" if stats.get("error_probability_model")=="Phred Q-score driven" else "Error rate"},
                                {"Property":"Q-score profile","Value":stats.get("qscore_profile","—")},
                                {"Property":"Mean payload Q","Value":f"{float(stats.get('mean_payload_q',0)):.2f}"},
                            ])
                        if parsed:
                            with st.expander("FASTQ preview",expanded=False):
                                preview=[]
                                for r in parsed[:3]: preview += ["@"+r["header"],r["sequence"],"+",r["quality"]]
                                st.code("\n".join(preview),language="text")
                            if st.session_state.get("fastq_origin")=="simulated":
                                inspect_read=int(st.number_input("Inspect simulated read",min_value=1,max_value=len(parsed),value=1,step=1,key="inspect_fastq_read"))-1
                                rr=parsed[inspect_read]
                                sid,_,assignment_source=_fastq_record_target(rr,design)
                                clean_ref=next((fastq_ecc.clean_dna(r.get("Full strand","")) for r in rows if int(r.get("Strand ID",-1))==sid),"") if sid is not None else ""
                                if clean_ref:
                                    _error_legend()
                                    st.markdown(f'<div class="seq-box">{_alignment_diff_html(clean_ref,rr["sequence"])}</div>',unsafe_allow_html=True)
                        collect("Download FASTQ",fq,"sequencing_reads.fastq","text/plain",width="stretch")

                        st.divider()
                        st.markdown("#### 🗳️ Consensus Reconstruction")
                        voting=st.radio("Consensus method",["Majority voting","Q-score likelihood voting"],horizontal=True,key="fq_voting")
                        if st.button("Run consensus reconstruction",type="primary",width="stretch",key="run_consensus"):
                            try:
                                sid_mode_for_consensus = (
                                    _current_hedges_sid_mode() if design=="TVH-V5.2 HEDges" and st.session_state.get("fastq_origin")=="simulated"
                                    else hb.SID_MODE_AUTO
                                )
                                cmap,srows,prows,cstats=independent_fastq_consensus(fq,design,voting,sid_mode=sid_mode_for_consensus)
                                erows,eprows,estats=_consensus_evaluation(rows,cmap,prows,voting)
                                st.session_state.update({
                                    "consensus_map":cmap,"consensus_rows":srows,"consensus_position_rows":prows,"consensus_stats":cstats,
                                    "consensus_eval_rows":erows,"consensus_eval_positions":eprows,"consensus_eval_stats":estats,"consensus_sid_mode":sid_mode_for_consensus if design=="TVH-V5.2 HEDges" else None,
                                })
                                st.session_state.pop("decode_result",None)
                            except Exception as exc:
                                st.error(f"Consensus failed: {exc}")
                        cmap=st.session_state.get("consensus_map") or {}
                        if cmap:
                            cs=st.session_state.get("consensus_stats",{})
                            c1,c2,c3=st.columns(3)
                            c1.metric("Usable reads",cs.get("usable_reads",0)); c2.metric("Covered strands",cs.get("covered_strands",0)); c3.metric("Voting",st.session_state.get("fq_voting","—"))
                            property_table([{"Property":"Read assignment","Value":cs.get("assignment","—")},{"Property":"Unassigned reads","Value":cs.get("unassigned_reads",0)}])
                            collect("Download consensus sequences",("\n".join(f">strand{sid}\n{seq}" for sid,seq in sorted(cmap.items()))+"\n").encode(),"consensus_dna.fasta","text/plain",width="stretch")

        if st.session_state.get("strand_rows") and st.session_state.get("consensus_map"):
            _render_consensus_analysis(st.session_state.get("strand_rows") or [])


# -----------------------------------------------------------------------------
# Panel 5 — Decode
# -----------------------------------------------------------------------------

def _with_intrinsic_restoration(result: dict) -> dict:
    """Attach a deterministic lossless inner-payload view without decode hints.

    `result["data"]` remains the exact byte stream reconstructed from DNA. If
    those bytes are themselves a recognized lossless wrapper (zlib/gzip/xz/
    bzip2, or a one-member ZIP), `restored_data` exposes the decompressed payload.
    No original file, SHA, compression-method session state, or filename is used.
    """
    data = result.get("data")
    if data is None:
        return result
    restored = restore_lossless_payload(data)
    if restored is None:
        return result
    result = dict(result)
    result.update({
        "stored_data": data,
        "stored_format": result.get("format"),
        "stored_extension": result.get("extension"),
        "stored_mime": result.get("mime"),
        "restored_data": restored["data"],
        "restored_name": restored["name"],
        "restored_format": restored["format"],
        "restored_extension": restored["extension"],
        "restored_mime": restored["mime"],
        "restoration_method": restored["method"],
        "restoration_validation": restored["validation"],
        "wrapper_format": restored["wrapper_format"],
    })
    return result


def _decode_legacy_from_rows(rows: List[dict], consensus: Dict[int,str] | None = None):
    design = st.session_state.get("dna_design")
    if consensus is None:
        payload_dna = "".join(str(r.get("Payload","")) for r in rows)
    else:
        parts=[]
        for r in rows:
            sid=int(r["Strand ID"])
            if sid not in consensus:
                raise ValueError(f"Missing consensus for strand {sid}")
            p=_legacy_payload_from_consensus(consensus[sid],r)
            if p is None:
                raise ValueError(f"Consensus strand {sid} is {len(consensus[sid])} nt; legacy decoding requires the consensus stage to restore the 125-nt strand length.")
            parts.append(p)
        payload_dna="".join(parts)
    data,bits,meta=legacy_codec.decode_dna(payload_dna,design,st.session_state.get("legacy_meta") or {})
    result = {"status":"OK","data":data,"bits":bits,"meta":meta,"format":_detect_type(data,st.session_state.get("stored_name","recovered.bin"))[0],"extension":Path(st.session_state.get("stored_name","recovered.bin")).suffix or ".bin","mime":st.session_state.get("stored_mime","application/octet-stream"),"strand_results":[]}
    return _with_intrinsic_restoration(result)


def _decode_hedges_bodies(bodies: Dict[int,str], max_edits: int, sid_mode: str = hb.SID_MODE_AUTO):
    results, full_bits, rec_status, chosen, candidates = hb.decode_bodies(bodies,max_edits=max_edits,sid_mode=sid_mode)
    out={"status":rec_status,"data":chosen.data if chosen else None,"bits":chosen.bits if chosen else full_bits,"full_bits":full_bits,
         "format":chosen.format if chosen else None,"extension":chosen.extension if chosen else ".bin","mime":chosen.mime if chosen else "application/octet-stream",
         "filename":chosen.filename if chosen else None,"domain":chosen.domain if chosen else None,
         "padding_bits":chosen.padding_bits if chosen else None,"reconstruction_hypothesis":chosen.hypothesis if chosen else None,
         "validation":chosen.validation if chosen else None,"reconstruction_candidates":[c.to_dict() for c in candidates],
         "strand_results":[r.__dict__ for r in results]}
    return _with_intrinsic_restoration(out)


def render_panel_5_decoding() -> None:
    with st.container(border=True):
        step_header(5, "Decoding")
        rows = st.session_state.get("strand_rows") or []
        design=st.session_state.get("dna_design") or "TVH-V5.2 HEDges"
        left,right=st.columns(2,gap="large")

        with left:
            st.markdown("#### 📥 Decoder input")
            sources=[]
            if rows:
                sources.append("Current clean strands")
            if st.session_state.get("direct_error_rows"):
                sources.append("Current direct-error strands")
            if st.session_state.get("consensus_map"):
                sources.append("Current FASTQ consensus")
            sources.append("Upload decode file")
            source=st.radio("Decode source",sources,horizontal=True,key="decode_source")

            upload_bytes=None
            upload_name=""
            upload_kind=None
            upload_sequences: Dict[int,str]={}
            import_voting="Q-score likelihood voting"

            if source=="Upload decode file":
                up=st.file_uploader(
                    "Upload full strands, consensus DNA, FASTQ, or FASTQ.gz",
                    type=["txt","fa","fasta","fna","fastq","fq","gz"],
                    key="decode_upload",
                )
                if up is not None:
                    upload_bytes=up.getvalue(); upload_name=up.name
                    if _looks_like_fastq(upload_bytes,upload_name):
                        upload_kind="FASTQ"
                        try:
                            recs=fastq_ecc.parse_fastq(upload_bytes)
                            lens=[len(r["sequence"]) for r in recs]
                            m1,m2=st.columns(2)
                            m1.metric("Detected input","FASTQ")
                            m2.metric("Reads",f"{len(recs):,}")
                            if lens:
                                property_table([{"Property":"Read length range","Value":f"{min(lens)}–{max(lens)} nt"}])
                        except Exception as exc:
                            st.error(f"FASTQ parse failed: {exc}")
                        property_table([{"Property":"Read assignment","Value":"Auto: embedded SI → explicit FASTQ header fallback"}])
                        import_voting=st.radio("Consensus method",["Majority voting","Q-score likelihood voting"],horizontal=True,key="decode_fastq_voting")
                    else:
                        try:
                            upload_sequences,fmt=_parse_fasta_or_dna_text(upload_bytes)
                            if design=="TVH-V5.2 HEDges":
                                _,preview_rows,kind=_hedges_bodies_from_sequence_map(upload_sequences, sid_mode=hb.SID_MODE_AUTO)
                                upload_kind=kind
                                m1,m2=st.columns(2)
                                m1.metric("Detected input",kind)
                                m2.metric("Sequences",len(upload_sequences))
                                if preview_rows:
                                    st.dataframe(pd.DataFrame(preview_rows),width="stretch",hide_index=True,height=min(260,38*len(preview_rows)+38))
                            else:
                                upload_kind=fmt
                                m1,m2=st.columns(2)
                                m1.metric("Detected input",fmt)
                                m2.metric("Sequences",len(upload_sequences))
                        except Exception as exc:
                            st.error(f"DNA input parse failed: {exc}")

            max_edits=2
            if design=="TVH-V5.2 HEDges":
                with st.expander("Advanced HEDges decoder",expanded=False):
                    max_edits=st.slider("Maximum edit radius per BODY (K)",0,4,2,1,key="hedges_k")

            if st.button("Decode and reconstruct",type="primary",width="stretch",key="run_decode"):
                try:
                    if design=="TVH-V5.2 HEDges":
                        parsed_rows=[]
                        current_sid_mode = _current_hedges_sid_mode()
                        decode_sid_mode = hb.SID_MODE_AUTO if source=="Upload decode file" else (
                            st.session_state.get("consensus_sid_mode", current_sid_mode) if source=="Current FASTQ consensus" else current_sid_mode
                        )
                        if source=="Current clean strands":
                            bodies={s.strand_id:s.body for s in (st.session_state.get("hedges_strands") or [])}
                        elif source=="Current direct-error strands":
                            bodies={int(r["Strand ID"]):fastq_ecc.clean_dna(r.get("Received BODY","")) for r in (st.session_state.get("direct_error_rows") or [])}
                        elif source=="Current FASTQ consensus":
                            bodies={int(k):fastq_ecc.clean_dna(v) for k,v in (st.session_state.get("consensus_map") or {}).items()}
                            if sorted(bodies)!=list(range(len(bodies))):
                                raise ValueError("Consensus strand IDs are not contiguous from 0")
                        else:
                            if not upload_bytes:
                                raise ValueError("Upload a decode file")
                            if upload_kind=="FASTQ":
                                cmap,srows,prows,cstats=independent_fastq_consensus(upload_bytes,"TVH-V5.2 HEDges",import_voting,sid_mode=hb.SID_MODE_AUTO)
                                bodies,parsed_rows,kind=_hedges_bodies_from_sequence_map({int(k):fastq_ecc.clean_dna(v) for k,v in cmap.items()}, sid_mode=hb.SID_MODE_AUTO)
                            else:
                                if not upload_sequences:
                                    upload_sequences,_=_parse_fasta_or_dna_text(upload_bytes)
                                bodies,parsed_rows,kind=_hedges_bodies_from_sequence_map(upload_sequences, sid_mode=hb.SID_MODE_AUTO)
                        result=_decode_hedges_bodies(bodies,int(max_edits),sid_mode=decode_sid_mode); result["source"]=source; result["parsed_rows"]=parsed_rows; result["sid_mode"]=decode_sid_mode
                        if source=="Upload decode file":
                            result["input_format"]=upload_kind
                            if upload_kind=="FASTQ":
                                result["consensus_method"]=import_voting
                    else:
                        if not rows:
                            raise ValueError("SM / Rule-base decoding requires the current strand design context")
                        if source=="Current clean strands":
                            result=_decode_legacy_from_rows(rows,None)
                        elif source=="Current direct-error strands":
                            cmap=_full_strand_map_from_direct_rows(st.session_state.get("direct_error_rows") or [])
                            result=_decode_legacy_from_rows(rows,cmap)
                        elif source=="Current FASTQ consensus":
                            result=_decode_legacy_from_rows(rows,st.session_state.get("consensus_map") or {})
                        else:
                            if not upload_bytes:
                                raise ValueError("Upload a decode file")
                            if upload_kind=="FASTQ":
                                cmap,_,_,_=independent_fastq_consensus(upload_bytes,design,import_voting)
                            else:
                                cmap=upload_sequences or _parse_fasta_or_dna_text(upload_bytes)[0]
                            result=_decode_legacy_from_rows(rows,cmap)
                        result["source"]=source
                        if source=="Upload decode file": result["input_format"]=upload_kind
                    st.session_state["decode_result"]=result
                except Exception as exc:
                    st.session_state["decode_result"]={"status":"FAIL","error":str(exc),"data":None,"strand_results":[],"source":source}

        with right:
            st.markdown("#### 📤 Decoded output")
            result=st.session_state.get("decode_result")
            if not result:
                st.info("Run decoding to see the recovered file.")
            else:
                status=result.get("status")
                if status in {"OK","OK_WEAK"}:
                    st.success(f"Reconstruction {status}")
                elif status=="AMBIGUOUS": st.warning("Reconstruction AMBIGUOUS")
                else: st.error(f"Reconstruction FAIL{': '+result.get('error','') if result.get('error') else ''}")
                if result.get("input_format"):
                    property_table([{"Property":"Input format","Value":result.get("input_format")},{"Property":"Decode source","Value":result.get("source","—")}])
                if design=="TVH-V5.2 HEDges" and result.get("strand_results"):
                    sr=result["strand_results"]
                    c1,c2,c3=st.columns(3)
                    c1.metric("HEDges OK",sum(r["status"]=="OK" for r in sr)); c2.metric("FAIL",sum(r["status"]=="FAIL" for r in sr)); c3.metric("AMBIGUOUS",sum(r["status"]=="AMBIGUOUS" for r in sr))
                    st.dataframe(pd.DataFrame(sr),width="stretch",hide_index=True)
                data=result.get("data")
                if data is not None:
                    ext=result.get("extension") or ".bin"
                    outname=result.get("filename") or ("recovered"+ext)
                    display_data=result.get("restored_data") if result.get("restored_data") is not None else data
                    display_name=result.get("restored_name") or outname
                    display_mime=result.get("restored_mime") or result.get("mime") or "application/octet-stream"
                    display_format=result.get("restored_format") or result.get("format") or "Unknown"
                    c1,c2,c3=st.columns(3)
                    c1.metric("Recovered bytes",f"{len(display_data):,}")
                    c2.metric("Detected format",display_format)
                    c3.metric("Padding removed",f"{result.get('padding_bits',0)} bits" if result.get("padding_bits") is not None else "—")
                    props=[]
                    if result.get("reconstruction_hypothesis"):
                        props.extend([
                            {"Property":"Framing hypothesis","Value":result.get("reconstruction_hypothesis")},
                            {"Property":"Validation","Value":result.get("validation") or "—"},
                        ])
                    if result.get("restored_data") is not None:
                        props.extend([
                            {"Property":"Recovered wrapper","Value":result.get("wrapper_format") or result.get("format") or "—"},
                            {"Property":"Lossless restoration","Value":result.get("restoration_method") or "—"},
                        ])
                    if props:
                        property_table(props)
                    _preview_bytes(display_data,display_name,display_mime,"Reconstructed file",key="decode_reconstructed_text_preview")
                    collect("Download reconstructed file",display_data,display_name,display_mime,width="stretch",key="download_reconstructed_file")
                    if result.get("restored_data") is not None:
                        with st.expander("Recovered stored representation",expanded=False):
                            property_table([
                                {"Property":"Stored format","Value":result.get("stored_format") or result.get("format") or "—"},
                                {"Property":"Stored bytes","Value":f"{len(data):,} B"},
                            ])
                            collect("Download recovered stored representation",data,outname,result.get("stored_mime") or result.get("mime","application/octet-stream"),width="stretch",key="download_recovered_wrapper")
                if result.get("reconstruction_candidates"):
                    with st.expander("Autonomous reconstruction hypotheses",expanded=False):
                        st.dataframe(pd.DataFrame(result["reconstruction_candidates"]),width="stretch",hide_index=True)


# -----------------------------------------------------------------------------
# Panel 6 — Summarization
# -----------------------------------------------------------------------------

def _bit_compare(a: str,b: str) -> dict:
    n=max(len(a),len(b)); err=sum((i>=len(a) or i>=len(b) or a[i]!=b[i]) for i in range(n))
    return {"bits":n,"errors":err,"ber":err/max(1,n),"accuracy":1-err/max(1,n)}


def render_panel_6_analysis() -> None:
    with st.container(border=True):
        step_header(6,"Summarization")
        result=st.session_state.get("decode_result")
        if not result:
            st.info("Complete Decoding first."); return
        original=st.session_state.get("input_bytes",b"")
        stored=st.session_state.get("stored_bytes",b"")
        recovered=result.get("restored_data") if result.get("restored_data") is not None else result.get("data")
        design=st.session_state.get("dna_design")

        st.markdown("#### File comparison")
        c1,c2,c3=st.columns(3,gap="large")
        with c1:
            st.markdown("##### Original")
            _preview_bytes(original,st.session_state.get("input_name","input.bin"),st.session_state.get("input_mime"),"Original",key="summary_original_text_preview")
            property_table([{"Property":"Size","Value":f"{len(original):,} B"}])
        with c2:
            st.markdown("##### Stored / encoded representation")
            _preview_bytes(stored,st.session_state.get("stored_name","stored.bin"),st.session_state.get("stored_mime"),"Stored",key="summary_stored_text_preview")
            property_table([{"Property":"Method","Value":st.session_state.get("compression_method","—")},{"Property":"Size","Value":f"{len(stored):,} B"}])
        with c3:
            st.markdown("##### Decoded")
            if recovered is not None:
                decoded_name=result.get("restored_name") or ("recovered"+(result.get("extension") or ".bin"))
                decoded_mime=result.get("restored_mime") or result.get("mime")
                decoded_format=result.get("restored_format") or result.get("format") or "Unknown"
                _preview_bytes(recovered,decoded_name,decoded_mime,"Decoded",key="summary_decoded_text_preview")
                property_table([{"Property":"Size","Value":f"{len(recovered):,} B"},{"Property":"Detected format","Value":decoded_format}])
            else:
                st.error("No uniquely reconstructed output is available.")

        st.markdown("#### Compression analysis")
        comp_rows=st.session_state.get("compression_analysis_rows") or []
        selected=next((r for r in comp_rows if r.get("Method")==st.session_state.get("compression_method")),None)
        compression_report=[
            {"Metric":"Selected method","Value":st.session_state.get("compression_method","—")},
            {"Metric":"Original size","Value":f"{len(original):,} B"},
            {"Metric":"Stored size","Value":f"{len(stored):,} B"},
            {"Metric":"Compression ratio","Value":f"{(selected or {}).get('Compression ratio',len(original)/max(1,len(stored))):.3f}×"},
            {"Metric":"Size reduction","Value":f"{(selected or {}).get('Size reduction (%)',(1-len(stored)/max(1,len(original)))*100):.2f}%"},
        ]
        if selected and selected.get("SSIM") is not None:
            ps=selected.get("PSNR")
            compression_report += [
                {"Metric":"PSNR","Value":"∞" if ps==float("inf") else f"{ps:.2f} dB"},
                {"Metric":"SSIM","Value":f"{selected.get('SSIM'):.4f}"},
                {"Metric":"MAE","Value":f"{selected.get('MAE'):.3f}"},
            ]
        property_table(compression_report)

        st.markdown("#### Encode-decode analysis")
        strand_rows=st.session_state.get("strand_rows") or []
        if design=="TVH-V5.2 HEDges":
            hs=st.session_state.get("hedges_strands") or []
            dna_len=sum(len(x.body) for x in hs); full_len=sum(len(x.full_strand) for x in hs)
            gc=float(np.mean([x.body_gc for x in hs])) if hs else 0; hp=max([x.body_hp for x in hs] or [0])
        else:
            dna=st.session_state.get("legacy_dna",""); dna_len=len(dna); full_len=sum(len(r.get("Full strand","")) for r in strand_rows); gc=legacy_codec.gc_content(dna) if dna else 0; hp=legacy_codec.longest_homopolymer(dna) if dna else 0
        encoding_report=[
            {"Metric":"DNA design","Value":DESIGN_DISPLAY.get(design,design)},
            *([{"Metric":"SID mode","Value":_current_hedges_sid_mode()}] if design=="TVH-V5.2 HEDges" else []),
            {"Metric":"Encoded DNA length","Value":f"{dna_len:,} nt"},
            {"Metric":"Full strand DNA length","Value":f"{full_len:,} nt"},
            {"Metric":"GC content","Value":f"{gc:.4f}"},
            {"Metric":"Longest homopolymer","Value":hp},
            {"Metric":"Strand count","Value":len(strand_rows)},
        ]
        rand_rows=st.session_state.get("encoding_rand_rows") or []
        rdf=pd.DataFrame(rand_rows) if rand_rows else pd.DataFrame()
        if not rdf.empty and "Base entropy" in rdf:
            ok=rdf[rdf.get("Analysis status","")=="OK"] if "Analysis status" in rdf else rdf
            if not ok.empty:
                encoding_report += [
                    {"Metric":"Mean strand base entropy","Value":f"{pd.to_numeric(ok['Base entropy'],errors='coerce').mean():.4f}"},
                    {"Metric":"Mean strand GC deviation","Value":f"{pd.to_numeric(ok['GC deviation'],errors='coerce').mean():.4f}"},
                    {"Metric":"Mean strand |γ|","Value":f"{pd.to_numeric(ok['|γ|'],errors='coerce').mean():.4f}"},
                    {"Metric":"Mean strand Alpha RMSE","Value":f"{pd.to_numeric(ok['Alpha RMSE'],errors='coerce').mean():.4f}"},
                ]
        property_table(encoding_report)

        if st.session_state.get("direct_error_stats"):
            st.markdown("#### Direct sequence error report")
            ds=st.session_state.get("direct_error_stats") or {}
            property_table([
                {"Metric":"Substitutions","Value":ds.get("sub",0)},
                {"Metric":"Insertions","Value":ds.get("ins",0)},
                {"Metric":"Deletions","Value":ds.get("del",0)},
                {"Metric":"Total errors","Value":ds.get("errors",0)},
                {"Metric":"Affected strand set","Value":f"{ds.get('strands',0)} strands"},
            ])

        if st.session_state.get("consensus_map"):
            st.markdown("#### FASTQ / ECC report")
            es=st.session_state.get("consensus_eval_stats") or {}
            sim=st.session_state.get("fastq_sim_stats") or {}; ss=sim.get("summary",{})
            property_table([
                {"Metric":"FASTQ source","Value":st.session_state.get("fastq_origin","—")},
                {"Metric":"Error probability model","Value":("Q-score" if sim.get("error_probability_model")=="Phred Q-score driven" else "Error rate" if sim else "—")},
                {"Metric":"Q-score profile","Value":sim.get("qscore_profile","—") if sim else "—"},
                {"Metric":"Consensus method","Value":es.get("voting_method",st.session_state.get("fq_voting","—"))},
                {"Metric":"FASTQ reads","Value":ss.get("records",(st.session_state.get("consensus_stats") or {}).get("records","—"))},
                {"Metric":"Injected substitutions","Value":ss.get("sub","—")},
                {"Metric":"Injected insertions","Value":ss.get("ins","—")},
                {"Metric":"Injected deletions","Value":ss.get("del","—")},
                {"Metric":"Dropped strands","Value":ss.get("dropped","—")},
                {"Metric":"Consensus accuracy","Value":f"{100*float(es.get('consensus_accuracy',0)):.4f}%"},
                {"Metric":"Mean consensus confidence","Value":f"{100*float(es.get('mean_consensus_probability',0)):.2f}%"},
                {"Metric":"Consensus DNA edit distance","Value":es.get("consensus_edit_distance","—")},
                {"Metric":"Low-confidence positions","Value":es.get("low_confidence_positions","—")},
                {"Metric":"Tie positions","Value":es.get("tie_positions","—")},
            ])

        st.markdown("#### Recovery Quality Report")
        if recovered is None:
            property_table([{"Metric":"Decode status","Value":result.get("status","FAIL")},{"Metric":"Recovery class","Value":"Failed / ambiguous"}])
            return
        stored_sha=hashlib.sha256(stored).hexdigest(); recovered_sha=hashlib.sha256(recovered).hexdigest(); exact=stored_sha==recovered_sha
        cmp=_bit_compare(bytes_to_bits(stored),bytes_to_bits(recovered))
        rec_report=[
            {"Metric":"Decode source","Value":result.get("source","—")},
            {"Metric":"Decode status","Value":result.get("status","—")},
            {"Metric":"Output status","Value":"Openable / reconstructed"},
            {"Metric":"Payload bit accuracy","Value":f"{100*cmp['accuracy']:.6f}%"},
            {"Metric":"BER","Value":f"{cmp['ber']:.6e}"},
            {"Metric":"Checksum","Value":"Pass" if exact else "Fail"},
            {"Metric":"Recovery class","Value":"Exact" if exact else "Recovered but not exact"},
        ]
        if design=="TVH-V5.2 HEDges":
            rec_report += [
                {"Metric":"Autonomous framing hypothesis","Value":result.get("reconstruction_hypothesis","—")},
                {"Metric":"Padding removed","Value":f"{result.get('padding_bits',0)} bits" if result.get("padding_bits") is not None else "—"},
            ]
        property_table(rec_report)
        if exact: st.success("100% stored-payload recovery — SHA-256 match")
        else: st.error("Recovered output does not exactly match the selected stored representation.")

        if design=="TVH-V5.2 HEDges" and result.get("full_bits") and "?" not in result.get("full_bits","") and strand_rows:
            useful=len(stored)*8; n=len(strand_rows)
            d1,d2,d3=st.columns(3)
            d1.metric("Useful BODY density",f"{useful/(n*76):.4f} bit/nt")
            d2.metric("SI + BODY density",f"{useful/(n*84):.4f} bit/nt")
            d3.metric("Full strand density",f"{useful/(n*125):.4f} bit/nt")


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def render_app() -> None:
    st.set_page_config(page_title=APP_TITLE,page_icon="🧬",layout="wide")
    apply_style()
    st.markdown('<div class="hero-card"><div class="hero-title">🧬 DNA Data Storage System</div></div>',unsafe_allow_html=True)
    render_stepper()
    render_panel_1_upload()
    render_panel_2_compression()
    render_panel_3_encoding()
    render_panel_4_experiment()
    render_panel_5_decoding()
    render_panel_6_analysis()


if __name__=="__main__":
    render_app()
