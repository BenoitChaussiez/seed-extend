import argparse

def parse():
    parser = argparse.ArgumentParser(description='Alignement de reads sur génome')
    parser.add_argument('--genome', required=True, help='Fichier FASTA du génome')
    parser.add_argument('--reads', required=True, help='Fichier FASTQ.gz des reads')
    parser.add_argument('--out', required=True, help='Fichier de sortie')
    parser.add_argument('-k', type=int, default=20, help='Taille des k-mers')
    parser.add_argument('--score-min', type=int, default=50, help='Score minimum')