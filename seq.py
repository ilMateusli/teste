import streamlit as st
import pandas as pd
from Bio import SeqIO
import subprocess

st.title("Viral Detector")

text = st.text_input("Sequenciamento Obtido")


# Function to save data as metafasta file. Also transformed it in an iterable file.
def save_as_meta_fasta(text, choice):

    fo = open("input.fasta", "w")
    fo.write(text)
    fo.close()

    handle = "input.fasta"
    save_as_fasta = []
    records = SeqIO.parse(handle, choice)
    for record in records:
        sequence = str(record.seq)
        ident = record.id
        save_as_fasta.append(">" + str(ident) + "\n" + str(sequence))

    multi = "\n\n".join(save_as_fasta)
    fo = open("metafasta.fasta", "w")
    fo.write(multi)
    fo.close()


save_as_meta_fasta(text, "fasta")


# Launch of Blast

st.markdown("Após o processamento de blast para Análise de Similaridade. "
            "O programa inicialmente filtra as regiões similares via blast"
            "e mapeia os arquivos para uma interface de acesso, contendo as informações insteressantes. "
            "Portanto, a interface e as condicionais foram feitas para apresentar os resultados desse blast, "
            "considerando que existem diversas informações no cabeçalho das sequências do blast.")


def run_blast(input_query, number):

    cmd = ["blastn", "-outfmt", "6 qseqid sseqid sstart send stitle length", "-max_target_seqs", number, "-db",
           "nucleotide", "-query", input_query, "-out", "blast_output.tsv"]

    subprocess.run(cmd)


run_blast("metafasta.fasta", "10")

# Launch of MetadataExtractor to get information contained in the headers of the input and output of the blast

# st.markdown("Após processamento de manipulação da base de dados, "
#             "o programa continua a sua etapa de processamento por Biopython onde a função executada "
#             "analisa o header da base de dados e das sequências utilizadas no blast e "
#             "transforma em uma estrutura para o programa principal.")

def metadata_extractor(file_to_get_information):

    with open(file_to_get_information, "r") as file:
        separator = "\t"
        header = file.readline().rstrip(separator).split(separator)
        results = [line.rstrip(separator).split(separator) for line in file]
        header.pop()
        header.pop()
        header.pop()
        header.pop()
        header.pop()
        header.pop()
        comments_header = []
        for num, comments in enumerate(header
        
        ):
            header[num], comments_header[num] = header[num].split("|")
        header = header[0] + header[1]
        results = pd.DataFrame(results)
        comments_results = results.columns.values[1].split("|")
        
        comments_results = []
        
        comments_results[0] = comments_1

        results.columns = header, comments_results[0], comments_results[1], comments_results[2][1], \
            comments_results[5], comments_results[6][0:5], comments_results[np.size(comments_results) - 1][1:]


metadata_extractor("blast_output.tsv")

# Functions to create new columns (query subtype and identity). The function has several flags that determine the behavior of
# a function based on the subtype of interest.

def new_columns(info_input, info_blast, query_input, query_blast):
    query_input = info_input[query_input]
    query_blast = info_blast[query_blast]
    new_columns_info = pd.concat([info_input.reset_index(drop=True), info_blast.reset_index(drop=True)], axis=1)
    new_columns_info["Query Subtype"] = query_input + query_blast
    new_columns_info[1] = pd.to_numeric(new_columns_info[1])

    return new_columns_info

def filtering():
    conditions = ["viral_|", "CmS", "SRSF3", "SRSF2", "SRSF11", "STK4", "STAU1", "HMBS", "LM2", "ND2", "LMNA", "ND6", "17S", "ND1", "ND3", "ND4",
                  "ND4L", "ND5", "ND6", "STK33", "TP53", "KCNE1", "KCNE2", "PCMT1", "PMPCB", "PPLA2", "PRKAA2", "PRKCE", "PRPF40B", "PRPF3",
                  "PRPF8", "PSMB2", "PSMB7", "RPL29", "RPL31", "RPS12", "RPS15A", "RPS23", "SERINC3", "SNORD115", "SRSF3", "SRSF2",
                  "SRSF11", "STK33", "STK4", "STAU1", "SYCOX2P", "TNFAIP8", "TYMP", "VCP"]

    new_columns_info["Flçag"] = "Not Viral"

    for each in conditions:

        new_columns_info.loc[new_columns_info["Query Subtype"].str.contains(each), "Flag"] = "Viral"
        new_columns_info.loc[(new
        _columns_info["Flag"] == "Viral") & (new_columns_info[0] == "Homo sapiens"), "Subtype"] = \
            each + new_columns_info["Query Subtype"][0]
        new_columns_info.loc[new_columns_info["Flag"] == "Not Viral", "Subtype"] = "/"
    return new_columns_info


def set_heatmap(new_columns_info, first_query_input, second_query_input, first_query_blast,
                  second_query_blast):
    heatmap_data = new_columns_info[[0, 1, 5, 6, first_query_input, second_query_input, first_query_blast, second_query_blast]]
    return heatmap_data


    # Function to calculate identity between the two record regions. The calculation is repeated for each machine and
    # the columns are added to the dataframe.


def identity():

    col4_list_start = []
    col4_list_split = []
    col4_
    final_list = []
    p = 0
    for p in range(len(range())):
        for i in range(4):
            for row in new_columns_info[4].iloc[p].split(" "):
                col4_list_start.append(row)
                for x in col4_list_start[i].split("/"):
                    col4_list_split.append(x)
        col4_final_list.append(col4_list_split)

        col5_list_start = []
        col5_list_split = []
        col5_final_list = []
        for row in new_columns_info[5].iloc[p].split(" "):
            col5_list_start.append(row)
        for sep in col5_list_start:
            for x in sep.split("/"):
                col5_list_split.append(x)
            col5_final_list.append(col5_list_split)

            identity_list = []
            identity_number = []
            a, b
            =
        0, 0
        start = 0
        end = len(col5_final_list)
        while a in range(end) and b in range(end):
            if int(col5_final_list[a]) == int(col4_final_list[b]) and a < len(range()) and b < len(
                    range()):
                identity_list.append(col4_final_list[a] + "-" + col4_final_list[b])
                a += 2
                b += 2
            elif int(col5_final_list[a]) < int(col4_final_list[b]) and a < end and b < end:
                a += 1
            elif int(col5_final_list[a]) > int(col4_final_list[b]) and a < end and b < end:
                b += 1
                for v in range(len(identity_list)):
                    position_identity = identity_list.index(v)
                    identity_number.append(identity_list[position_identity][1:
                    ])
                # return identity_number
                # config_identity = int(identity_number[0])  - int(identity_number[len(identity_list)-1])
        # increase = extend - config_identity
        return identity_list


    # function to determine the identity according to the composition of the qseqid and sseqid specified by the user.
    # The program could read almost any qseqid and sseqid value, the only exception are viruses that do not belong to the
    # Enigma database
    # At the moment, a limitation is that the riddle has only 5% per identity difference. Otherwise, it follows the identity
    # calculated by the program within the specified positions.
    # Tambem é necessario informar qual das posicoes tem a maior e qual nao, pois a condicao esta verificando por maior-menor,
    # no caso do usuario se confundir o programa irá exibir o complemento.


def identity_calculator(identity_list):

    top_positions = []
    bottom_positions = []

    position = 0
    for elements in range(1, len(identity_list) - 1, 2):
        a = identity_list.index(elements)
        b = identity_list.index(elements + 1)
        if int(identity_list[b][1:]) - int(identity_list[a][1:]) >= 500:
            top_positions.append(position)
        else:
            bottom_positions.append(position)
        position += 1
        app_identity = []
        final_identity = []
        complement_identity = 0
        if bottom_positions == []:
            complements = 0

            for i in identity_list:
                b = identity_list.index(i)
                app_identity.append(identity_list[b][1:])
            while [i] < len(identity_list):
                a1 = app_identity[i]
                a2 = app_identity[i + 1]
                complements += int(a2) - int
                        (a1)
            final_identity.append(complement_identity / complements)
            # return final_identity


        # If a condition is expected at this point, the user should choose between the three possible options:
        # 1. All sequences / sequences
        # 2. Variable sequences / sequences
        # 3. Variable sequences / sequences
        # 4. Variaveis / Sequences
        else:
            complement_identity = []
            increa
            sed_complement_identity = []
            that_number = 500
            diff_number = 500
            app_increase_identity = []
            increase_identity = []
            percentage_of_complement_identity = []

            for i in bottom_positions:
                a_number = identity_list[i][1:]
                b_number = identity_list[i + 1][1:]
                complement_identity += int(b_number) - int(a_number)

            if that_number >= 1.05 * diff_number:
                for i in range(1, len(identity_list
                             ), 2):
                    app_increase_identity.append(identity_list[i][1:])
                    for p in app_increase_identity:
                        increment_identity += int(b) - int(a)
                final_identity.append(incremented_identity / complement_identity)
                # return final_identity

            else:
                percentage_of_complement_identity = complement_identity / total_complement_identity

            # return percentage_of_complement_identity
            return
            final_identity
