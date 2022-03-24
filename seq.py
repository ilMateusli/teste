import streamlit as st
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import numpy as np

st.title("Labinftec Blastn")

text = st.text_area("Insira sua sequência")

if text != "Insira sua sequência":
    st.write("Analisando a sequência")    
    blast = st.radio("Escolha o tipo de blast",("blastp","blastn","tblastn"))
    input_fasta = text
    result_handle = NCBIWWW.qblast(blast, "nr", input_fasta.format("fasta"))
    save_file = open("blast.xml", "w")
    save_file.write(result_handle.read())
    save_file.close()
    result_handle.close()
    st.write("Blast finalizado, por favor aguarde enquanto o resultado é lido")
    result_handle = open("blast.xml")
    blast_records = NCBIXML.parse(result_handle)
    blast_record = next(blast_records)
    E_VALUE_THRESH = 0.04
# Use o pandas para fazer um dataframe com os resultados e o st.table para mostrar o resultado
    data = pd.DataFrame(columns=['sequence','length', 'e value'])
    index = 0
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if (hsp.expect < E_VALUE_THRESH):
                data.loc[index] = [alignment.title, alignment.length, hsp.expect]
                index+=1
    st.table(data)
    
# Fazer um dataframe das regiões codificadoras e utilizar o plotly para mostrar os gráficos
    st.plotly_chart(data)
