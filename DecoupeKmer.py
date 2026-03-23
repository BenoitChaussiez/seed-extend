read = "GTTGCTTCCCGTAGGAGTTTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTACTGATCGAAGTCTTGGTGAGCCGTTACCTCACCAACAAACTAATCAGACGCATCCCCATCCAATACCGAAATTCTTTAATGT"

k = 25
step = 5

#for i in range(0, len(read) - k + 1, step):
    #kmer = read[i:i+k]
    #print(f"{i:3d}: {kmer}")
import gzip
from Bio import SeqIO
def extraire_kmers_read_biopython(fichier_fastq, num_read=1, k=25, step=5):
    """
    Extrait les k-mers d'un read FASTQ.gz.
    
    Args:
        fichier_fastq: Chemin du fichier FASTQ.gz
        num_read: Numéro du read (1-indexé)
        k: Taille des k-mers
        step: Pas d'échantillonnage
    
    Returns:
        Dict avec 'nom', 'sequence', 'kmers'
    """
    with gzip.open(fichier_fastq, 'rt') as f:
        for i, record in enumerate(SeqIO.parse(f, 'fastq')):
            if i == num_read - 1:  
                sequence = str(record.seq).upper()
                
                kmers = []
                for j in range(0, len(sequence) - k + 1, step):
                    kmer = sequence[j:j+k]
                    
                    kmers.append(kmer)
                
                return {
                    'nom': record.id,
                    'sequence': sequence,
                    'kmers': kmers
                }
fichier = 'SRR10971381_1_100_premiers.fastq.gz'

resultat = extraire_kmers_read_biopython(fichier, num_read=1, k=25, step=5)
print(resultat['kmers'])