import argparse

def parse():
    parser = argparse.ArgumentParser(description='Alignement de reads sur génome')
    parser.add_argument('-g','--genome', required=True, help='Fichier FASTA du génome')
    parser.add_argument('-r','--reads', required=False, help='Fichier FASTQ.gz des reads')
    parser.add_argument('-o','--out', required=True, help='Fichier de sortie')
    parser.add_argument('-k', type=int, default=20, help='Taille des k-mers')
    parser.add_argument('-i', type=int, default=70, help="Pourcentage d'identité ou un alignement est accepté")

    return parser.parse_args()