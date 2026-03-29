import gzip
from Bio import SeqIO
def extraire_kmers_read_biopython(fichier_fastq, num_read=1, k=12, step=1):
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

def extraire_kmers_fichier(fichier, num_read=1, k=12, step=1, format='fastq'):
    """
    Extrait les k-mers d'un read depuis un fichier FASTQ/FASTA (compressé ou non).
    
    Args:
        fichier: Chemin du fichier (peut être .fastq, .fastq.gz, .fasta, .fasta.gz)
        num_read: Numéro du read (1-indexé)
        k: Taille des k-mers
        step: Pas d'échantillonnage
        format: Format du fichier ('fastq' ou 'fasta')
    
    Returns:
        Dict avec 'nom', 'sequence', 'kmers'
    """
    if fichier.endswith('.gz'):
        opener = gzip.open
        mode = 'rt'
    else:
        opener = open
        mode = 'r'
    
    with opener(fichier, mode) as f:
        for i, record in enumerate(SeqIO.parse(f, format)):
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
    
    return None

def extraire_kmers_sequence(sequence, k, step):
    """
    Extrait les k-mers d'une séquence en mémoire.
    
    Args:
        sequence: La séquence ADN (string)
        k: Taille des k-mers
        step: Pas d'échantillonnage
    
    Returns:
        List de k-mers

    Doctest:
    >>> ADN= 'ATGCGC'
    >>> extraire_kmers_sequence(ADN,3,3)
    ['ATG', 'CGC']
    """
    kmers = []
    for j in range(0, len(sequence) - k + 1, step):
        kmers.append(sequence[j:j+k])
    return kmers
