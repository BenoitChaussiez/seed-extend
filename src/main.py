from parser import *
from Preprocess import filtrer_fastq
from Suffixe_table import *
from Alignement import *
from DecoupeKmer import extraire_kmers_read_biopython, extraire_kmers_fichier , extraire_kmers_sequence
from Donnees_sim import *
import time
import psutil
import os
from Bio import SeqIO
import gzip

def main_test(genome, k, id_threshold, output_file, align_threshold=0.8):
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
                            aligné, align_pos, score, strand = align_read(read, genome_sequence, pos, id_threshold)
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


def main_reel(genome, reads, k, id_threshold, output_file, max_reads=10000):  
    start = time.time()
    
    genome = read_genome(genome)
    print(f"Génome lu: {len(genome)} bases")

    table_suffixes = suffixe_table(genome)
    print(f"Table des suffixes construite: {len(table_suffixes)} entrées")

    fichier_source = reads
    reads_list = []
    
    print(f"\nLecture des {max_reads} premiers reads du fichier: {fichier_source}")
    
    if fichier_source.endswith('.gz'):
        import gzip
        opener = gzip.open
        mode = 'rt'
    else:
        opener = open
        mode = 'r'
    
    with opener(fichier_source, mode) as f:
        for i, record in enumerate(SeqIO.parse(f, 'fastq')):
            if i >= max_reads:
                break
            reads_list.append(str(record.seq))
    
    print(f"Reads chargés: {len(reads_list)}")
    
    nb_read_alignés = 0
    nb_forward = 0
    nb_reverse = 0
    
    with open(output_file, "w") as out_f:
        out_f.write("# Resultats d'alignement\n")
        out_f.write("# Format: read_number\tposition\tscore\tstrand\tread\n")
        out_f.write("# Pour chaque read, on teste le brin direct (+) et le brin inverse (-)\n\n")
        
        for num_read, read in enumerate(reads_list, start=1):
            if num_read % 1000 == 0:
                print(f"Traitement read {num_read}/{len(reads_list)}, alignés: {nb_read_alignés}...")
            
            found_alignment = False
            
            read_rc = reverse_complement(read)
            
            for read_seq, strand in [(read, '+'), (read_rc, '-')]:
                kmers = extraire_kmers_sequence(read_seq, k, step=5)
                
                for idx, kmer in enumerate(kmers):
                    positions = est_present_dicho(table_suffixes, kmer, genome)
                    
                    if positions:
                        for pos in positions[:5]:
                            aligné, approx_pos, score, align_strand = align_read(read_seq, genome, pos, id_threshold)
                            
                            if aligné:
                                out_f.write(f"{num_read}\t{approx_pos}\t{score}\t{strand}\t{read}\n")
                                nb_read_alignés += 1
                                if strand == '+':
                                    nb_forward += 1
                                else:
                                    nb_reverse += 1
                                found_alignment = True
                                break
                        if found_alignment:
                            break
                
                if found_alignment:
                    break
        
        end = time.time()
        process = psutil.Process(os.getpid())
        mem = process.memory_info().rss
        
        out_f.write("# RESUME\n")
        out_f.write(f"Total reads analyses: {len(reads_list)}\n")
        out_f.write(f"Reads alignes: {nb_read_alignés}\n")
        if len(reads_list) > 0:
            out_f.write(f"Taux d'alignement: {nb_read_alignés / len(reads_list):.2%}\n")
            out_f.write(f"Reads alignes sur brin direct (+): {nb_forward}\n")
            out_f.write(f"Reads alignes sur brin inverse (-): {nb_reverse}\n")
        out_f.write(f"Temps d'execution: {end - start:.2f} secondes\n")
        out_f.write(f"Utilisation memoire: {mem / (1024 * 1024):.2f} MB\n")
        out_f.write(f"Parametres: k={k}, step=5, max_reads={max_reads}\n")
    
    print(f"\n{'='*60}")
    print(f"Résultats sauvegardés dans: {output_file}")
    print(f"Reads alignés: {nb_read_alignés}/{len(reads_list)} ({nb_read_alignés/len(reads_list):.2%})")
    print(f"  - Brin direct (+): {nb_forward}")
    print(f"  - Brin inverse (-): {nb_reverse}")
    print(f"{'='*60}")

if __name__ == "__main__":
    args = parse()
    
    if hasattr(args, 'reads') and args.reads:
        main_reel(args.genome, args.reads, args.k, args.i, args.out)
    else:
        main_test(args.genome, args.k, args.i, args.out)


