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
    result_handle =
    NCBIWWW.qblast(blast, "nr", input_fasta)
    save_file = open("blast.xml", "w")
    save_file.write(result_handle.read())
    save_file.close()
    result_handle.close()

    st.write("Blast finalizado, por favor aguarde enquanto o resultado é lido")

    record = SeqIO.read("blast.xml", "blast-xml")
    result = []
    for alignment in record.alignments:
        for hsp in alignment.hsps:
            result.append([alignment.hit_id, alignment.length, hsp.expect, "..." + hsp.query[-35:] + "[" + str(hsp.match) + "]" + hsp.sbjct[:35] + "...", alignment.title, alignment.accession])
    result = pd.DataFrame(result, columns = ["Identificador", "Tamanho (aminoácidos)", "E-value", "Alinhamento", "Título", "Acesso"])
    result = result.sort_values("E-value").reset_index(drop=True)

    st.table(result)
    st.write("Fim do programa")
