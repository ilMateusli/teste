import streamlit as st
import pandas as pd
import numpy as np
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline

# Criando caixa de texto para inserir sequência
st.title('Blastn e Blastx')
st.write('Insira sua sequência:')
seq = st.text_area('')

# Criando botão para rodar o Blastn
if st.button('Blastn'):
    st.write('Rodando Blastn...')
    # Criando arquivo fasta
    rec = SeqRecord(Seq(seq, IUPAC.unambiguous_dna), id='seq')
    SeqIO.write(rec, 'seq.fasta', 'fasta')
    # Rodando o blastn
    result_handle = NCBIWWW.qblast("blastn", "nt", 'seq.fasta')
    # Salvando o resultado
    save_file = open("my_blast.xml", "w")
    save_file.write(result_handle.read())
    save_file.close()
    result_handle.close()
    # Lendo o resultado
    result_handle = open("my_blast.xml")
    blast_records = NCBIXML.parse(result_handle)
    # Imprimindo o resultado
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                st.write('****Alignment****')
                st.write('sequence:', alignment.title)
                st.write('length:', alignment.length)
                st.write('e value:', hsp.expect)
                st.write(hsp.query[0:75] + '...')
                st.write(hsp.match[0:75] + '...')
                st.write(hsp.
                         subject[0:75] + '...')
