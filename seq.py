# Escreva um código que faça o Blastn e Blastx de sequências (que serão inseridas em uma caixa de texto) usando o Streamlit e Biopython:

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
seq = st.text_input('')

# Criando botão para rodar o Blastn
if st.button('Blastn'):
    st.write('Rodando Blastn...')
    # Transformando sequência em um objeto SeqRecord
    seq_obj = Seq(seq, IUPAC.unambiguous_dna)
    seq_record = SeqRecord(seq_obj)
    # Salvando o objeto em um arquivo FASTA
    SeqIO.write(seq_record, 'seq.fasta', 'fasta')
    # Rodando o Blastn
    blastn_cline = NcbiblastnCommandline(query='seq.fasta', db='nt', evalue=0.001,
                                         outfmt=5, out='blastn.xml')
    stdout, stderr = blastn_cline()
    # Abrindo o arquivo XML
    blast_result_handle = open('blastn.xml')
    blast_records = NCBIXML.parse(blast_result_handle)
    # Lendo o arquivo XML
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                st.write(alignment.title)
                st.write(hsp.expect)
                st.write(hsp.query[0:75] + '...')
                st.write(hsp.match[0:75] + '...')
                st.write(hsp.sbjct[0:75] + '...')

# Criando botão para rodar o Blastx
if st.button('Blastx'):
    st.write('Rodando Blastx...')
    # Transformando sequência em um objeto SeqRecord
    seq_obj = Seq(seq, IUPAC.unambiguous_dna)
    seq_record = SeqRecord(seq_obj)
    # Salvando o objeto em um arquivo FASTA
    SeqIO.write(seq_record, 'seq.fasta', 'fasta')
    # Rodando o Blastx
    blastx_cline = NcbiblastnCommandline(query='seq.fasta', db='nr', evalue=0.001,
                                         outfmt=5, out='blastx.xml')
    stdout, stderr = blastx_cline()
    # Abrindo o arquivo XML
    blast_result_handle = open('blastx.xml')
    blast_records = NCBIXML.parse(blast_result_handle)
    # Lendo o arquivo XML
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                st.write(alignment.title)
                st.write(hsp.expect)
                st.write(hsp.query[0:75] + '...')
                st.write(hsp.match[0:75] + '...')
                st.write(hsp.sbjct[0:75] + '...')
                
# Criando botão para limpar a tela
if st.button('Limpar'):
    st.write('Limpeza concluída')