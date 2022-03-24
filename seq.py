import streamlit as st
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import numpy as np

st.title("Labinftec Blastn")

text = st.text_area("Insira sua sequência").format("fasta")
if st.button("Analisar"):
    result_handle = NCBIWWW.qblast("blastn", "nr", text)
    save_file = open("blast.xml", "w")
    save_file.write(result_handle.read())
    save_file.close()
    result_handle.close()

    result_handle = open("blast.xml")
    blast_records = NCBIXML.parse(result_handle)
    blast_record = next(blast_records)
    E_VALUE_THRESH = 0.04

    df = pd.DataFrame(columns=['Sequência', 'Tamanho', 'E-value'])
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < E_VALUE_THRESH:
                df = df.append({'Sequência': alignment.title, 'Tamanho': alignment.length, 'E-value': hsp.expect}, ignore_index=True)
# st.table(df) para tabela
    st.write(df)
