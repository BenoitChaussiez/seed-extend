from parser import *
from Preprocess import filtrer_fastq
from Suffixe_table import *
from Alignement import *
from DecoupeKmer import extraire_kmers_read_biopython, extraire_kmers_fichier , extraire_kmers_sequence
from Donnees_sim import *
import time
import psutil
import os

def main_test(genome, k, output_file, align_threshold=0.8):
    start = time.time()
    
    # Lire le génome
    genome_sequence = read_genome(genome)
    print(f"Génome lu: {len(genome_sequence)} bases")
    
    # Construire la table des suffixes
    table_suffixes = suffixe_table(genome_sequence)
    print(f"Table des suffixes construite: {len(table_suffixes)} entrées")
    
    # Simuler des reads sur les deux brins
    reads = generate_reads(genome_sequence, read_length=150, n_reads=1000, both_strands=True)
    print(f"Reads simulés: {len(reads)}")
    
    # Aligner les reads
    nb_read_alignés = 0
    nb_read_total = len(reads)
    nb_forward = 0
    nb_reverse = 0
    
    with open(output_file, "w") as f:
        for num_read, read in enumerate(reads, start=1):
            # Extraire les k-mers
            read_rc = reverse_complement(read)

            for read_seq in [read, read_rc]:
                kmers = extraire_kmers_sequence(read_seq, k, step=5)

                for kmer in kmers:
                    positions = est_present_dicho(table_suffixes, kmer, genome_sequence)
                    if positions:
                        for pos in positions:
                            # align_read retourne (aligné, position, score)
                            aligné, align_pos, score, strand = align_read(read, genome_sequence, pos)
                            if aligné:
                                # Déterminer le brin (optionnel, pour les statistiques)
                                # On peut vérifier si le read original ou son reverse complement est dans le génome
                                if strand == '+':
                                    nb_forward += 1
                                else:
                                    nb_reverse += 1
                                
                                line = f"{num_read}\t{align_pos}\t{score}\t{strand}\t{read}\n"
                                f.write(line)
                                nb_read_alignés += 1
                                found_alignment = True
                                break
                        if found_alignment:
                            break
        
        end = time.time()
        process = psutil.Process(os.getpid())
        mem = process.memory_info().rss
        
        f.write("# Résumé\n")
        f.write(f"Nombre total de reads: {nb_read_total}\n")
        f.write(f"Nombre de reads alignés: {nb_read_alignés}\n")
        if nb_read_total > 0:
            f.write(f"Taux d'alignement: {nb_read_alignés / nb_read_total:.2%}\n")
            f.write(f"Reads alignés sur brin direct (+): {nb_forward}\n")
            f.write(f"Reads alignés sur brin inverse (-): {nb_reverse}\n")
        f.write(f"Seuil d'alignement: {align_threshold:.0%}\n")
        f.write(f"Temps d'exécution: {end - start:.2f} secondes\n")
        f.write(f"Utilisation mémoire: {mem / (1024 * 1024):.2f} MB\n")
    
    print(f"\nRésultats écrits dans {output_file}")
    print(f"Alignement: {nb_read_alignés}/{nb_read_total} ({nb_read_alignés/nb_read_total:.2%})")
    print(f"  - Brin direct (+): {nb_forward}")
    print(f"  - Brin inverse (-): {nb_reverse}")


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


