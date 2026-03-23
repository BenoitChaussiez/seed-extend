from parser import *
from Preprocess import *
from Suffixe_table import *
from Lecture import *
from Alignement import *
from DecoupeKmer import *
from Donnees_sim import *
import time
import psutil
import os

def main_test(genome, k, output_file):
    start = time.time()
    # Lire le génome
    genome = read_genome(genome)

    # Construire la table des suffixes
    table_suffixes = suffixe_table(genome)

    #Simuler des reads à partir du génome de gA
    reads = generate_reads(genome, read_length=150, n_reads=1000)

    # Aligner les reads sur les génomes
    nb_read_alignés = 0
    nb_read_total = len(reads)
    # Écrire les résultats dans un fichier de sortie
    with open("output.tsv", "w") as f:
        for read in reads:
            num_read = 1
            kmers = extraire_kmers_read_biopython(fichier_source, num_read, k, step=5)['kmers']
            for kmer in kmers:
                positions = est_present_dicho(table_suffixes, kmer, genome)
                if positions:
                    for pos in positions:
                        aligné, approx_pos, score = align_read(read, genome, pos)
                        if aligné:
                            line = f"{num_read}\t{approx_pos}\t{score}\t{read}\n"
                            nb_read_alignés += 1
                            break  # Si on trouve un alignement pour ce read, on peut passer au read suivant
            num_read += 1
        f.write(line)
            

        end = time.time()
        process = psutil.Process(os.getpid())
        mem = process.memory_info().rss  # en octets
        f.write("# Résumé\n")
        f.write(f"Nombre total de reads: {nb_read_total}\n")
        f.write(f"Nombre de reads alignés: {nb_read_alignés}\n")
        f.write(f"Probabilité d'alignement: {nb_read_alignés / nb_read_total:.2%}\n")
        f.write(f"Temps d'exécution: {end - start:.2f} secondes\n")
        f.write(f"Utilisation mémoire: {mem / (1024 * 1024):.2f} MB\n")

    print(f"Résultats écrits dans {output_file}")

def main_reel(genome, reads, k, output_file):
    start = time.time()
    # Lire le génome
    genome = read_genome(genome)

    # Construire la table des suffixes
    table_suffixes = suffixe_table(genome)

    #Simuler des reads à partir du génome de gA
    fichier_source = reads
    fichier_sortie = 'SRR10971381_1_qualite20.fastq.gz'
    reads_list = filtrer_fastq(fichier_source, fichier_sortie, qualite_min=20)
    # Aligner les reads sur les génomes
    nb_read_alignés = 0
    nb_read_total = len(reads_list)
    # Écrire les résultats dans un fichier de sortie
    with open("output.tsv", "w") as f:
        for read in reads_list:
            num_read = 1
            kmers = extraire_kmers_read_biopython(fichier_source, num_read, k, step=5)['kmers']
            for kmer in kmers:
                positions = est_present_dicho(table_suffixes, kmer, genome)
                if positions:
                    for pos in positions:
                        aligné, approx_pos, score = align_read(read, genome, pos)
                        if aligné:
                            line = f"{kmers['nom']}\t{approx_pos}\t{score}\t{read}\n"
                            nb_read_alignés += 1
                            break  # Si on trouve un alignement pour ce read, on peut passer au read suivant
            num_read += 1
        f.write(line)
            

        end = time.time()
        process = psutil.Process(os.getpid())
        mem = process.memory_info().rss  # en octets
        f.write("# Résumé\n")
        f.write(f"Nombre total de reads: {nb_read_total}\n")
        f.write(f"Nombre de reads alignés: {nb_read_alignés}\n")
        f.write(f"Probabilité d'alignement: {nb_read_alignés / nb_read_total:.2%}\n")
        f.write(f"Temps d'exécution: {end - start:.2f} secondes\n")
        f.write(f"Utilisation mémoire: {mem / (1024 * 1024):.2f} MB\n")

    print(f"Résultats écrits dans {output_file}")

if __name__ == "__main__":
    args = parse()
    main_test(args.genome, args.k, args.out)
    #main_reel(args.genome, args.reads, args.k, args.out)

