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

    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < E_VALUE_THRESH:
                st.write("****Alignment****")
                st.write("sequence:", alignment.title)
                st.write("length:", alignment.length)
                st.write("e value:", hsp.expect)
