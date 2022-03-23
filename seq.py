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

# Para o resultado do blast aparecer na tela, é necessário rodar o código e clicar no botão "Blastn"
# Criando caixa de texto para inserir sequência
st.title('Blastn e Blastx')
st.write('Insira sua sequência:')
seq = st.text_area('')

# Criando botão para rodar o Blastx
if st.button('Blastx'):
    st.write('Rodando Blastx...')
    result_handle = NCBIWWW.qblast("blastx", "nr", seq)
    blast_records = NCBIXML.read(result_handle)
    st.write(blast_records)

# Criando botão para rodar o Blastn
if st.button('Blastn'):
    st.write('Rodando Blastn...')
    result_handle = NCBIWWW.qblast("blastn", "nt", seq)
    blast_records = NCBIXML.read(result_handle)
    st.write(blast_records)
