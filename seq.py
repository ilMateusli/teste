import streamlit as st
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import numpy as np

st.title("Labinftec Blastn")

st.write("""
    # Analisador de sequências
    Este analisador utiliza o sistema de blast para analisar uma sequência de DNA e retornar o resultado do blast.

    **Para utilizar o analisador, basta inserir uma sequência e clicar no botão.**
    """)

st.write("Insira uma sequência de DNA para ser analisada:")
input_fasta = st.text_area("", "CAACCGCTATTCCTCTTTTTGCAGGGGTTTTTCAAAATTATCAAGTTCCTCTTTTATCAGTATATGTTCAAGCTGCAAATTTACATTTATCAGTTTTGAGAGATGTTTCAGTGTTTGGACAAAGAAGACCTTTTAATATAGGGATAAATAATCAACAACGGCCTAGCCTCCCAGGTTTATCTGTTCTTGACGGGACAGAATTTGCTTATGGGACCTCCTCAAATTTGCCATCCGCTGTATACAGAAAAAGCGGAACGGTAGATTCGCTGGAT")

if st.button("Analisar sequência"):
    result_handle = NCBIWWW.qblast("blastn", "nr", input_fasta.format("fasta"))
    save_file = open("blast.xml", "w")
    save_file.write(result_handle.read())
    save_file.close()
    result_handle.close()
    st.write("Blast finalizado, por favor aguarde enquanto o resultado é lido")
    result_handle = open("blast.xml")
    blast_records = NCBIXML.parse(result_handle)
    blast_record = next(blast_records)
    E_VALUE_THRESH = 0.04
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < E_VALUE_THRESH:
                data = {'sequence': alignment.title, 'length': alignment.length, 'e value': hsp.expect, 'sbjct': hsp.sbjct[0:75] + "..."}
                df = pd.DataFrame(data, columns = ['sequence', 'length', 'e value', 'sbjct'])
                st.write("****Alignment****")
                st.dataframe(df)

st.write("""
    ### Referências:
    * https://www.ncbi.nlm.nih.gov/books/NBK279671/
    * https://streamlit.io/docs/
""")
