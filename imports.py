"""Explicit input package reader. Evaluation references stay outside decoding."""
import io
import zipfile
import pandas as pd


def package(raw):
    allowed={'R1.fastq','R2.fastq','decoder_settings.json','evaluation/reference_strands.csv','evaluation/expected_hash.txt','evaluation/channel_config.json'}
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        members=z.infolist()
        if len({m.filename for m in members})!=len(members):raise ValueError('Duplicate ZIP members')
        if any(m.filename not in allowed for m in members):raise ValueError('Use sequencing_output.zip, not All_results.zip')
        if sum(m.file_size for m in members)>100_000_000:raise ValueError('Package exceeds 100 MB')
        files={m.filename:z.read(m) for m in members}
    if 'R1.fastq' not in files:raise ValueError('Package needs R1.fastq')
    if 'evaluation/expected_hash.txt' in files:
        value=files['evaluation/expected_hash.txt'].decode('ascii').strip()
        if len(value)!=64 or set(value.lower())-set('0123456789abcdef'):raise ValueError('Invalid expected SHA-256')
    return files


def references(raw):
    frame=pd.read_csv(io.BytesIO(raw),keep_default_na=False)
    if not {'Strand ID','Payload'}.issubset(frame.columns):raise ValueError('Reference needs Strand ID and Payload columns')
    if len(frame)>65536:raise ValueError('Too many reference strands')
    ids=pd.to_numeric(frame['Strand ID'],errors='raise')
    if any(ids<0) or any(ids%1) or ids.duplicated().any():raise ValueError('Reference strand IDs must be unique nonnegative integers')
    frame['Strand ID']=ids.astype(int)
    for col in ('Payload','Filler'):
        if col not in frame:frame[col]=''
        if any(not isinstance(s,str) or len(s)>600 or set(s)-set('ACGT') for s in frame[col]):raise ValueError('Invalid reference DNA')
    return frame.to_dict('records')
