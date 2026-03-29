import random

def generate_reads(genome, read_length=150, n_reads=1000, both_strands=True):
    """
    Génère des reads aléatoires à partir du génome.
    
    Args:
        genome: Séquence du génome
        read_length: Longueur des reads
        n_reads: Nombre de reads à générer
        both_strands: Si True, génère des reads sur les deux brins (50% chaque)
    
    Returns:
        List de reads (strings)
    """
    reads = []
    genome_length = len(genome)
    
    def reverse_complement(seq):
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
        return ''.join(complement.get(base, 'N') for base in reversed(seq))
    
    for _ in range(n_reads):
        start_pos = random.randint(0, genome_length - read_length)
        read = genome[start_pos:start_pos + read_length]
        
        if both_strands and random.choice([True, False]):
            read = reverse_complement(read)
        
        reads.append(read)
    
    return reads