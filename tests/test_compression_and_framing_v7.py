import io

import pytest

import compression_pipeline as cp
import framing
import hedges_bridge as hb
import legacy_codec as lc


def test_general_compressor_is_self_contained_and_exact():
    raw=(b"DNA storage compression benchmark\n"*200)
    _best,candidates=cp.run_compression_benchmark("notes.txt",raw)
    methods={c.method for c in candidates}
    assert {"zstd_level19","xz_preset9"} <= methods
    for candidate in candidates:
        restored=framing.restore_lossless_payload(candidate.data)
        assert restored is not None
        assert restored["data"]==raw
        assert restored["name"]=="notes.txt"


def test_image_lossy_codecs_have_three_named_levels():
    Image=pytest.importorskip("PIL.Image")
    image=Image.new("RGB",(48,48),(40,120,200))
    buf=io.BytesIO();image.save(buf,"BMP")
    _best,candidates=cp.run_compression_benchmark("picture.bmp",buf.getvalue())
    levels={}
    for candidate in candidates:
        if candidate.lossy:
            levels.setdefault(candidate.codec_family,set()).add(candidate.quality_tier)
    assert levels["JPEG"]=={"Light","Medium","Strong"}
    assert levels["WebP"]=={"Light","Medium","Strong"}


def test_hedges_embeds_name_domain_length_crc_and_padding():
    raw=b"arbitrary bytes\x00\xff"*7
    strands,meta=hb.encode_bytes(raw,filename="result.bin",mime="application/octet-stream")
    decoded,bits,status,chosen,candidates=hb.decode_bodies(
        {strand.strand_id:strand.body for strand in strands},
        max_edits=0,
        sid_mode=meta["sid_mode"],
    )
    assert status=="OK"
    assert chosen is not None and chosen.data==raw
    assert chosen.filename=="result.bin"
    assert chosen.domain=="binary"
    assert chosen.padding_bits==meta["padding_bits"]
    assert "CRC32" in chosen.validation


def test_public_rinf_is_shuffled_implementation():
    assert lc.DISPLAY["RINF_SHUFFLED_B16"]=="R∞"
    data=b"R infinity public mapping"
    dna,_,meta=lc.encode_bytes(data,"RINF_SHUFFLED_B16")
    assert lc.decode_dna(dna,"RINF_SHUFFLED_B16",meta)[0]==data
    assert dna!=lc.encode_bytes(data,"RINF_B16")[0]


def test_ui_exposes_only_hedges():
    import original_ui as ui
    assert "RINF_B16" not in ui.DESIGNS
    assert ui.DESIGNS==["TVH-V5.2 HEDges"]
