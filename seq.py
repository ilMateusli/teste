st.title("Viral Detector")
st.write("Streamlit + BioPython")

with open("blast.html", "r") as file:
   html = file.read()

# Inserção do texto
text = st.text_area("Insira as sequências a serem analisadas", "")

# Convertendo o texto para um arquivo fasta
if text != "":
    file = open("sequencias.fasta", "w")
    file.write(text)
    file.close()

# Converter o arquivo fasta em uma string
sequences = ""
for seq_record in SeqIO.parse("sequencias.fasta", "fasta"):
    sequences += str(seq_record.seq)
sequences = sequences.replace("\n", "")

# Realizando o blast
if sequences != "":
    file = open("sequencias.fasta", "w")
    file.write(sequences)
    file.
    close()
    subprocess.call(["blastn", "-query", "sequencias.fasta", "-db", "refseq_rna", "-out", "sequencias.txt"])


# Gerando o html com o resultado do blast
with open("sequencias.txt", "r") as file:
   html = file.read()

st.write(html, unsafe_allow_html=True)
