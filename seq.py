import streamlit as st
import pandas as pd
from Bio import SeqIO
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
    with open("my_blast.xml", "w") as out_handle:
        out_handle.write(result_handle.read())
    result_handle.close()

    st.write("Salvando o resultado do blast, caso você queira abrir depois")
    st.balloons()
    result_handle = open("my_blast.xml")
    from Bio.Blast import NCBIXML
    blast_record = NCBIXML.read(result_handle)
    E_VALUE_THRESH = 0.01
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < E_VALUE_THRESH:
                st.write(alignment.title)
                st.write(hsp.expect)
                st.write(hsp.score)
                st.write(hsp.query)
                st.write(hsp.match)
                st.write(hsp.sbjct)
                time.sleep(2)
    result_handle.close()
    st.write("O arquivo foi salvo e pode ser baixado em my_blast.xml")
