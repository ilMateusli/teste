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
    input_fasta = """>seq1
{}
""".format(text)

    with open("fasta.fasta","w") as fasta:
        fasta.write(input_fasta)

    for seq_record in SeqIO.parse("fasta.fasta","fasta"):
        print(seq_record.id)
        print(repr(seq_record.seq))
        print(len(seq_record))

        result_handle = NCBIWWW.qblast(blast, "nr", seq_record.seq)
        with open("my_blast.xml", "w") as out_handle:
            out_handle.write(result_handle.read())
        result_handle.close()

        result_handle = open("my_blast.xml")
        blast_records = NC
        BIWWW.parse(result_handle)

        E_VALUE_THRESH = 0.04
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect < E_VALUE_THRESH:
                        st.write("****Alignment****")
                        st.write("sequence:", alignment.title)
                        st.write("length:", alignment.length)
                        st.write("e value:", hsp.expect)
                        st.write(hsp.query[0:75] + "...")
                        st.write(hsp.match[0:75] + "...")
                        st.write(hsp.sbjct[0:75] + "...")
                        st.balloons()
                        st.balloons()
        st.success("Sequência Analisada")
