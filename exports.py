"""One visible home for secondary downloads; session-local bytes only."""
import streamlit as st


def collect(label,data,file_name=None,mime=None,**kwargs):
    name=file_name or 'output.bin'
    if hasattr(data,'getvalue'):data=data.getvalue()
    if isinstance(data,str):data=data.encode()
    st.session_state.setdefault('all_exports',{})[name]=data
    return False
