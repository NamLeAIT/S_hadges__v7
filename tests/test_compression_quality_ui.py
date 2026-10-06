import io
import numpy as np
import pytest
from PIL import Image
from pathlib import Path
from streamlit.testing.v1 import AppTest
import original_ui as ui


def png(array):
    buf=io.BytesIO();Image.fromarray(array).save(buf,format='PNG');return buf.getvalue()


def test_quality_matches_local_ssim_reference_and_psnr():
    structural_similarity=pytest.importorskip('skimage.metrics').structural_similarity
    rng=np.random.default_rng(123)
    ref=rng.integers(0,256,(32,32,3),dtype=np.uint8)
    got=np.clip(ref.astype(int)+rng.integers(-20,21,ref.shape),0,255).astype(np.uint8)
    metrics=ui._image_quality_metrics(png(ref),png(got))
    assert metrics['SSIM']==pytest.approx(structural_similarity(ref,got,channel_axis=2,data_range=255),abs=1e-12)
    mse=np.mean((ref.astype(float)-got.astype(float))**2)
    assert metrics['PSNR']==pytest.approx(10*np.log10(255**2/mse))
    exact=ui._image_quality_metrics(png(ref),png(ref))
    assert exact=={'PSNR':float('inf'),'SSIM':1.0,'MAE':0.0}
    assert ui._image_quality_metrics(png(ref),png(got[:20])) is None


def test_image_compression_displays_quality_cards_and_charts():
    # Deterministic patterned image creates measurable lossy candidates.
    y,x=np.indices((40,40))
    array=np.stack([(x*9)%256,(y*9)%256,((x+y)*7)%256],axis=-1).astype(np.uint8)
    app=AppTest.from_file(str(Path(__file__).resolve().parents[1]/'app.py'),default_timeout=120)
    app.session_state['input_bytes']=png(array)
    app.session_state['input_name']='quality.png';app.session_state['input_mime']='image/png'
    app.run();assert not app.exception
    app.radio(key='storage_mode').set_value('Compression').run()
    next(b for b in app.button if b.label=='Run Compression').click().run();assert not app.exception
    assert {'PSNR (dB) ↑','SSIM ↑','MAE ↓'}<={m.label for m in app.metric}
    assert {'Quality comparison','Size vs quality','File size'}<={t.label for t in app.tabs}
    assert len(app.get('vega_lite_chart'))+len(app.get('arrow_vega_lite_chart'))>=4
    assert any('∞'==m.value for m in app.metric) or any(m.label=='PSNR (dB) ↑' and m.value not in ('—','∞') for m in app.metric)
