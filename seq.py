import streamlit as st
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
import numpy as np

text = st.text_area("Insira sua sequência")

if text != "Insira sua sequência":
    st.write("Analisando a sequência")
    st.balloons()

    blast = st.radio("Escolha o tipo de blast",("blastp","blastn","tblastn"))
    input_fasta = Seq(text)
    result_handle = NCBIWWW.qblast("blastn", "nr", input_fasta)

    st.write("O blast está sendo realizado, aguarde")
    st.balloons()
# O repositório é do GitHub e por isso, os arquivos são necessários para rodar
    with open("blast.xml", "w") as out_handle:
        out_handle.write(result_handle.read())
    result_handle.close()

    st.write("Carregando resultados")
    st.balloons()

    from Bio.Blast import NCBIXML
    result_handle = open("blast.xml")
    blast_records = NCBIXML.parse(result_handle)
    E_VALUE_THRESH = 0.01
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < E_VALUE_THRESH:
                    st.write("****Alignment****")
                    st.write("Sequence:", alignment.title)
                    st.write("Length:", alignment.length)
                    st.write("e value:", hsp.expect)
                    st.write(hsp.query[0:75] + "...")
                    st.write(hsp.match[0:75] + "...")
                    st.write(hsp.sbjct[0:75] + "...")
else:
    st.write("Insira a sequência acima")
