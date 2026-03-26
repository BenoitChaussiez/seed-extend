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


def main_reel(genome, reads, k, output_file, max_reads=10000):  
    start = time.time()
    
    genome = read_genome(genome)
    print(f"Génome lu: {len(genome)} bases")

    table_suffixes = suffixe_table(genome)
    print(f"Table des suffixes construite: {len(table_suffixes)} entrées")
    
    print(f"\nExtrait de la table des suffixes (10 premiers):")
    for i in range(min(10, len(table_suffixes))):
        suffix_pos = table_suffixes[i]
        print(f"  Position {suffix_pos}: {genome[suffix_pos:suffix_pos+30]}...")

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
    
    for num_read, read in enumerate(reads_list, start=1):
        print(f"\n{'='*60}")
        print(f"Read {num_read}: {read[:50]}...")
        print(f"Longueur: {len(read)}")
        
        kmers = extraire_kmers_sequence(read, k, step=5)
        print(f"Nombre de kmers (k={k}, step=5): {len(kmers)}")
        
        if kmers:
            print(f"Premier kmer: {kmers[0]}")
            
            if kmers[0] in genome:
                print(f"   Le kmer est dans le génome!")
                pos = genome.find(kmers[0])
                print(f"    Position: {pos}")
            else:
                print(f"  le kmer n'est PAS dans le génome")
                
                read_rc = reverse_complement(read)
                kmers_rc = extraire_kmers_sequence(read_rc, k, step=5)
                if kmers_rc and kmers_rc[0] in genome:
                    print(f"   Mais le reverse complement du kmer est dans le génome!")
                    pos = genome.find(kmers_rc[0])
                    print(f"    Position: {pos}")
                else:
                    print(f"  Le kmer n'est pas trouvé, même en reverse complement")
        
        found = False
        for idx, kmer in enumerate(kmers[:3]):  
            positions = est_present_dicho(table_suffixes, kmer, genome)
            print(f"  Kmer {idx} ('{kmer}'): {len(positions)} positions trouvées par dichotomie")
            
            if positions:
                print(f"    Positions: {positions[:5]}")
                for pos in positions[:1]:
                    print(f"\n  Test alignement à la position {pos}:")
                    aligné, approx_pos, score, strand = align_read(read, genome, pos)
                    print(f"    aligné={aligné}, score={score}, strand={strand}")
                    
                    if aligné:
                        print(f"   ALIGNEMENT RÉUSSI!")
                        found = True
                        break
            if found:
                break
        
        if not found:
            print(f"\n  AUCUN ALIGNEMENT trouvé pour ce read")

    print(f"\n{'='*60}")
    print("Diagnostic terminé")

if __name__ == "__main__":
    args = parse()
    
    if hasattr(args, 'reads') and args.reads:
        main_reel(args.genome, args.reads, args.k, args.out)
    else:
        main_test(args.genome, args.k, args.out)


