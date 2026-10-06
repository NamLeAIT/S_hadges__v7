from __future__ import annotations

import hashlib
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from dna_randomness_metrics import analyze_dna_randomness, DNASequenceError


def _f(v: Any, d: int = 4) -> str:
    try: return f"{float(v):.{d}f}"
    except Exception: return "—"


def _reference_line(fig: go.Figure, value: float, label: str) -> None:
    fig.add_hline(y=value, line_dash="dash", annotation_text=label, annotation_position="top left")


def _hash_payload(overall: str, strands: Iterable[Tuple[int, str, Optional[int]]]) -> str:
    h = hashlib.sha256(overall.encode())
    for sid, seq, bits in strands:
        h.update(f"|{sid}|{bits}|{seq}".encode())
    return h.hexdigest()


def _one_row(sid: int, seq: str, payload_bits: Optional[int]) -> Dict[str, Any]:
    seq = "".join(c for c in str(seq).upper() if c in "ACGT")
    if len(seq) < 4:
        return {"Strand ID": sid, "Length (nt)": len(seq), "Analysis status": "Too short"}
    try:
        report = analyze_dna_randomness(seq, original_bit_count=payload_bits, max_lag=min(20, len(seq)-1)).to_dict()
        q = report["sequence_quality"]; t = report["translational_active_particle"]; i = report["inverse_ising"]; l = report["logic_3input_1output"]
        return {
            "Strand ID": sid,
            "Length (nt)": int(q["length_nt"]),
            "Information density (bits/nt)": q.get("information_density_bits_per_nt"),
            "Base entropy": q["base_entropy_normalized"],
            "GC fraction": q["gc_fraction"],
            "GC deviation": q["gc_deviation_from_0_5"],
            "Longest homopolymer": int(q["longest_homopolymer"]),
            "HP segments ≥2": int(q["homopolymer_segments_ge_2"]),
            "V": t["velocity_V_L_per_nt"],
            "D": t["diffusion_D_L2_per_nt"],
            "|γ|": i["abs_polarization_gamma"],
            "|λ|": i["abs_ising_lambda"],
            "|h|": i["abs_ising_h"],
            "Neighbour corr.": i["nearest_neighbour_correlation"],
            "Alpha RMSE": l["alpha_rmse_from_0_5"],
            "Beta RMSE": l["beta_rmse_from_0_125"],
            "Analysis status": "OK",
        }
    except Exception as exc:
        return {"Strand ID": sid, "Length (nt)": len(seq), "Analysis status": f"Unavailable: {exc}"}


def _plot_strand_metric(df: pd.DataFrame, metric: str, key: str, ideal: float | None = None) -> None:
    if metric not in df.columns:
        return
    d = df[["Strand ID", metric]].copy()
    d[metric] = pd.to_numeric(d[metric], errors="coerce")
    d = d.dropna()
    if d.empty:
        return
    fig = go.Figure(go.Scatter(x=d["Strand ID"], y=d[metric], mode="lines+markers", name=metric))
    if ideal is not None:
        _reference_line(fig, ideal, f"Ideal = {ideal}")
    fig.update_layout(title=f"{metric} by strand", xaxis_title="Strand ID", yaxis_title=metric, height=350, margin=dict(l=40,r=20,t=55,b=40))
    st.plotly_chart(fig, width="stretch", key=key)


def _render_overall(report: Dict[str, Any], prefix: str) -> None:
    q=report["sequence_quality"]; t=report["translational_active_particle"]; i=report["inverse_ising"]; l=report["logic_3input_1output"]
    m1,m2,m3,m4,m5=st.columns(5)
    m1.metric("Base entropy", _f(q.get("base_entropy_normalized")))
    m2.metric("GC deviation", _f(q.get("gc_deviation_from_0_5")))
    m3.metric("|Polarization γ|", _f(i.get("abs_polarization_gamma")))
    m4.metric("Alpha RMSE", _f(l.get("alpha_rmse_from_0_5")))
    m5.metric("Beta RMSE", _f(l.get("beta_rmse_from_0_125")))
    st.caption("Entropy closer to 1 is better; GC deviation, |γ|, |λ| and α/β RMSE closer to 0 indicate less structural bias.")

    tab1,tab2,tab3,tab4,tab5=st.tabs(["Base composition","Translational MSD","Inverse Ising","Logic α / β","Summary"])
    with tab1:
        vals=[q.get("a_fraction",0),q.get("c_fraction",0),q.get("g_fraction",0),q.get("t_fraction",0)]
        fig=go.Figure(go.Bar(x=list("ACGT"), y=vals, text=[f"{x:.3f}" for x in vals], textposition="auto"))
        _reference_line(fig,0.25,"Ideal = 0.25")
        fig.update_layout(title="Nucleotide composition",xaxis_title="Nucleotide",yaxis_title="Fraction",yaxis_range=[0,max(.35,max(vals)*1.2)])
        st.plotly_chart(fig,width="stretch",key=f"{prefix}_base")
    with tab2:
        fig=go.Figure()
        fig.add_trace(go.Scatter(x=t.get("lags_nt",[]),y=t.get("msd",[]),mode="lines+markers",name="Observed MSD"))
        fig.add_trace(go.Scatter(x=t.get("lags_nt",[]),y=t.get("msd_fitted",[]),mode="lines",name="Fitted MSD"))
        fig.update_layout(title="Translational active-particle MSD",xaxis_title="Lag (nt)",yaxis_title="Mean squared displacement")
        st.plotly_chart(fig,width="stretch",key=f"{prefix}_msd")
        a,b,c=st.columns(3); a.metric("Velocity V",_f(t.get("velocity_V_L_per_nt"),6)); b.metric("Diffusion D",_f(t.get("diffusion_D_L2_per_nt"),6)); c.metric("Fit R²",_f(t.get("fit_r_squared"),4))
    with tab3:
        labels=["|γ|","|λ|","|h|","|Neighbour corr.|"]
        vals=[abs(float(i.get("abs_polarization_gamma",0))),abs(float(i.get("abs_ising_lambda",0))),abs(float(i.get("abs_ising_h",0))),abs(float(i.get("nearest_neighbour_correlation",0)))]
        fig=go.Figure(go.Bar(x=labels,y=vals,text=[f"{x:.4f}" for x in vals],textposition="auto"))
        fig.update_layout(title="Local binary structure",yaxis_title="Absolute value")
        st.plotly_chart(fig,width="stretch",key=f"{prefix}_ising")
        st.caption("λ and h are transparent pseudo-likelihood estimates; they are not claimed to be numerically identical to the unavailable manuscript SI-S4 closed-form expression.")
    with tab4:
        contexts=[f"{x:03b}" for x in range(8)]
        alpha=[l.get("alpha",{}).get(c) for c in contexts]
        beta=[l.get("beta",{}).get(c,0) for c in contexts]
        fig=go.Figure(go.Bar(x=contexts,y=alpha)); _reference_line(fig,.5,"Ideal = 0.5")
        fig.update_layout(title="Conditional output balance α",xaxis_title="3-bit context",yaxis_title="P(output=0 | context)",yaxis_range=[0,1])
        st.plotly_chart(fig,width="stretch",key=f"{prefix}_alpha")
        fig=go.Figure(go.Bar(x=contexts,y=beta)); _reference_line(fig,.125,"Ideal = 0.125")
        fig.update_layout(title="Context occurrence balance β",xaxis_title="3-bit context",yaxis_title="Context fraction")
        st.plotly_chart(fig,width="stretch",key=f"{prefix}_beta")
    with tab5:
        rows=[
            ("DNA length",int(q.get("length_nt",0)),"nt"),("Base entropy normalized",q.get("base_entropy_normalized"),""),
            ("GC fraction",q.get("gc_fraction"),""),("GC deviation from 0.5",q.get("gc_deviation_from_0_5"),""),
            ("Longest homopolymer",int(q.get("longest_homopolymer",0)),"nt"),("Homopolymer segments ≥2",int(q.get("homopolymer_segments_ge_2",0)),"segments"),
            ("Translational velocity V",t.get("velocity_V_L_per_nt"),"L/nt"),("Translational diffusion D",t.get("diffusion_D_L2_per_nt"),"L²/nt"),
            ("|Polarization γ|",i.get("abs_polarization_gamma"),""),("|Inverse-Ising λ|",i.get("abs_ising_lambda"),""),
            ("|Inverse-Ising h|",i.get("abs_ising_h"),""),("Alpha RMSE from 0.5",l.get("alpha_rmse_from_0_5"),""),("Beta RMSE from 0.125",l.get("beta_rmse_from_0_125"),""),
        ]
        if "information_density_bits_per_nt" in q:
            rows.insert(1,("Information density",q["information_density_bits_per_nt"],"bits/nt"))
        st.dataframe(pd.DataFrame(rows,columns=["Metric","Value","Unit"]),width="stretch",hide_index=True)


def render_randomness_analysis(
    overall_dna: str,
    strand_sequences: List[Tuple[int, str, Optional[int]]],
    *,
    state_prefix: str,
    overall_bit_count: Optional[int] = None,
    sequence_label: str = "Payload/BODY",
) -> None:
    """Full App1-style randomness analysis plus a per-strand extension."""
    if not overall_dna:
        return
    current_hash=_hash_payload(overall_dna,strand_sequences)
    report_key=f"{state_prefix}_rand_overall"; rows_key=f"{state_prefix}_rand_rows"; hash_key=f"{state_prefix}_rand_hash"
    with st.container(border=True):
        st.markdown("#### Randomization Analysis")
        st.caption(f"Per-strand randomness is evaluated on the encoded {sequence_label}; primers/SI are excluded so fixed structural regions do not bias the encoding comparison.")
        if st.button("Run Analysis",type="primary",width="stretch",key=f"{state_prefix}_run_rand"):
            try:
                overall=analyze_dna_randomness(overall_dna,original_bit_count=overall_bit_count,max_lag=min(20,len(overall_dna)-1)).to_dict()
                rows=[_one_row(sid,seq,bits) for sid,seq,bits in strand_sequences]
                st.session_state[report_key]=overall; st.session_state[rows_key]=rows; st.session_state[hash_key]=current_hash
            except Exception as exc:
                st.warning(f"Randomization analysis unavailable: {exc}")
        if st.session_state.get(hash_key) not in (None,current_hash) and st.session_state.get(report_key):
            st.info("Encoded DNA changed. Run Analysis again to refresh randomness results.")
        overall=st.session_state.get(report_key); rows=st.session_state.get(rows_key) or []
        if not overall or not rows:
            return

        tab_a,tab_b=st.tabs(["Per-strand randomness","Overall encoded DNA"])
        with tab_a:
            df=pd.DataFrame(rows)
            ok=df[df.get("Analysis status","")=="OK"].copy() if "Analysis status" in df else df.copy()
            if not ok.empty:
                m1,m2,m3,m4,m5=st.columns(5)
                m1.metric("Mean entropy",_f(pd.to_numeric(ok["Base entropy"],errors="coerce").mean()))
                m2.metric("Mean GC deviation",_f(pd.to_numeric(ok["GC deviation"],errors="coerce").mean()))
                m3.metric("Max homopolymer",int(pd.to_numeric(ok["Longest homopolymer"],errors="coerce").max()))
                m4.metric("Mean |γ|",_f(pd.to_numeric(ok["|γ|"],errors="coerce").mean()))
                m5.metric("Mean α RMSE",_f(pd.to_numeric(ok["Alpha RMSE"],errors="coerce").mean()))
            st.dataframe(df,width="stretch",hide_index=True,height=min(450,38*len(df)+38))
            if not ok.empty:
                r1,r2=st.columns(2)
                with r1:
                    _plot_strand_metric(ok,"Base entropy",f"{state_prefix}_strand_entropy",ideal=1.0)
                    _plot_strand_metric(ok,"Longest homopolymer",f"{state_prefix}_strand_hp")
                    _plot_strand_metric(ok,"V",f"{state_prefix}_strand_v")
                    _plot_strand_metric(ok,"|γ|",f"{state_prefix}_strand_gamma",ideal=0.0)
                    _plot_strand_metric(ok,"Alpha RMSE",f"{state_prefix}_strand_alpha",ideal=0.0)
                with r2:
                    _plot_strand_metric(ok,"GC fraction",f"{state_prefix}_strand_gc",ideal=0.5)
                    _plot_strand_metric(ok,"GC deviation",f"{state_prefix}_strand_gcdev",ideal=0.0)
                    _plot_strand_metric(ok,"D",f"{state_prefix}_strand_d")
                    _plot_strand_metric(ok,"|λ|",f"{state_prefix}_strand_lambda",ideal=0.0)
                    _plot_strand_metric(ok,"Beta RMSE",f"{state_prefix}_strand_beta",ideal=0.0)
            csv=df.to_csv(index=False).encode("utf-8")
            from exports import collect
            collect("Download per-strand randomness CSV",csv,"strand_randomness_metrics.csv","text/csv",width="stretch",key=f"{state_prefix}_rand_csv")
        with tab_b:
            _render_overall(overall,state_prefix)
