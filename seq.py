import streamlit as st
import pandas as pd
from Bio import SeqIO
from Bio.Blast import NCBIWWW
import numpy as np

st.title('Pipeline Detecção Viral')

# Importando arquivos de texto.
# Obs. - Todos os arquivos devem estar salvo no mesmo diretório,
# uma vez que o programa apenas realiza a análise. 


st.subheader('Trabalho Prático II - Cinética Química e Bioquímica - CA064')
st.subheader('Integrantes do grupo: Carlos, Julio e Ylia')

# Adicionando informação no sidebar (lateral).
# st.sidebar.header('Options')
st.sidebar.markdown('Todos os arquivos de texto, presentes na pasta do projeto.')
st.sidebar.info('Os arquivos FASTA estão no diretório do projeto.')


# Criando função simples para testes, obtendo sequencias de outro arquivo.
def create_fasta(sequence, name_fasta):
    records = []
    for i, line in enumerate(sequence.splitlines()):
        records.append(SeqIO.SeqRecord(SeqIO.Seq(line.upper()), id=name_fasta))
    fasta_string = ">" + \
        str(records[0].id) + "\n" + \
        str(records[0].seq) + "-------\n"
    return fasta_string


# Criando função simples para testes, obtendo sequencias do blast.
def blast(sequence):
    handle = NCBIWWW.qblast('blastn', 'nt', sequence)
    return handle.read()

# Lembrar: Fazer a cobenção dos arquivos .txt no Streamlit
# Lembrar: Retirar opções disponiveis no menu lateral do programa. (Não usar)
# Lembrar: (Apagar) colocar fazer um select de qual arquivo de fasta eu escolho para rodar o Programa.


# Upload do arquivo Fasta.

@st.cache(persist=True)
def load_data_fasta():
    # Lista arquivos no diretório.
    arq_fasta = list(filter(lambda x: x.endswith('.fas'), open('..\\pipeline\\fasta', 'r')))
    # Concatenando todas as sequencias do Fasta em apenas um arquivo.

    all_sequence = []
    for i, seq in enumerate(arq_fasta):
        for single_record in SeqIO.parse(arq_fasta[i], "fasta"):
            all_sequence.append(">" +
                                str(single_record.id) + "\n" +
                                str(single_record.seq))
    # Retorna uma lista com as sequencias do Fasta.
    return all
    # Lembrar: Criar uma função para utilizar Entrez do Biopython, e receber os links dos arquivos para download
    # Lembrar: (Supor-Tall)retirar todas as linhas de cabeçalho apos a criação das listas.


# Realiza a requisição do Blast.
@st.cache
def load_data_blast(sequence):
    handle = NCBIWWW.qblast('tblastn', 'refseq_rna', sequence)
    return handle.read()


# Retornar resultado do blast em tabela para o streamlit
def blast_to_csv(blast_record):

    blast_records = NCBIWWW.parse(blast_record)
    # A cada record é criado um dict com as informações necessárias para a comparação.
    for blast_record in blast_records:
        alignment_list = []
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                alignment_dict = {
                    'alignment_id': alignment.title,
                    'alignment_length': alignment.length,
                    'query_start': hsp.query_start,
                    'query_end': hsp.query_end,
                    'hit_start': hsp.sbjct_start,
                    'hit_end': hsp.sbjct_end,
                    'expect': hsp.expect,
                    'score': hsp.score,
                    'bits': hsp.bits,
                    'identities': hsp.identities,
                    'positive': hsp.positives,
                    'gaps': hsp.gaps,
                    'query_seq': hsp.query,
                    'match_seq': hsp.match,
                    'sbjct_seq': hsp.sbjct
                }
                alignment_list.append(alignment_dict)
        return alignment_list
