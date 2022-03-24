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
    input_fasta = "ATCGATCGATCGTGATCGATCGTAGGTACGAGTG"
    result_handle = NCBIWWW.qblast("blastn", "nr", input_fasta.format("fasta"))
    save_file = open("blast.xml", "w")
    save_file.write(result_handle.read())
    save_file.close()
    result_handle.close()
    st.write("Blast finalizado, por favor aguarde enquanto o resultado é lido")
    st.write(result_handle)
