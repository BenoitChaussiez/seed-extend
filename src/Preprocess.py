import gzip
from Bio import SeqIO
import time
import argparse

def parse_preprocess():
    parser = argparse.ArgumentParser(description='Preprocessing of FASTQ files to filter out low-quality reads')
    parser.add_argument('-r','--reads', required=True, help='fastq files')
    parser.add_argument('-o','--out', required=True, help='Fichier de sortie')
    parser.add_argument('-q', '--qualite', type=int, default=20, help='Minimum quality score (default: 20)')
    return parser.parse_args()

def filtrer_fastq(entree, sortie, qualite_min=20):
    """    Filtre un fichier FASTQ compressé selon le score de qualité moyen
    
    Paramètres:
        entree (str): Chemin vers le fichier FASTQ.gz d'entrée
        sortie (str): Chemin vers le fichier FASTQ.gz de sortie
        qualite_min (int): Score Phred minimum (défaut: 20)
    
    Retourne:
        int: Nombre de reads conservés
        """
    import time
    start = time.time()
    
    total = 0
    conserves = 0
    
    with gzip.open(entree, 'rt') as f_in:
        with gzip.open(sortie, 'wt') as f_out:
            for i, record in enumerate(SeqIO.parse(f_in, 'fastq')):
                total += 1
                
                if total % 100_000 == 0:
                    elapsed = time.time() - start
                    rate = total / elapsed
                    print(f"\rReads: {total:,} | Vitesse: {rate:,.0f}/s | Temps: {elapsed:.0f}s", end="")
                
                qualites = record.letter_annotations.get("phred_quality")
                if qualites and sum(qualites)/len(qualites) >= qualite_min:
                    conserves += 1
                    SeqIO.write(record, f_out, 'fastq')
    
    print(f"\nTerminé: {total} reads en {time.time()-start:.1f}s")
    return conserves


def main():
    args = parse_preprocess()
    fichier=args.reads
    fichier_sortie=args.out
    qualité=args.qualite
    filtrer_fastq(fichier,fichier_sortie,qualité)
if __name__ =="__main__":
    main()