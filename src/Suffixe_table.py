def read_genome(file_path: str) -> str:
    """
    Lit un fichier FASTA et retourne la séquence génomique complète sous forme de chaîne de caractères.
    
    :param file_path: chemin vers le fichier FASTA
    :return: séquence génomique complète
    """
    genome = []
    
    with open(file_path, 'r') as fasta:
        for line in fasta:
            line = line.strip()
            if not line or line.startswith(">"):
                continue  
            genome.append(line)
    genome.append("$") 
    return ''.join(genome)

def suffixe_table(genome: str) -> list:
    """
    Construit la table des suffixes pour une séquence génomique donnée.
    
    :param genome: séquence génomique complète
    :return: table des suffixes (liste d'indices)
    >>> genome = "ATTGT$"
    >>> suffixe_table(genome)
    [5, 0, 3, 4, 2, 1]
    """
    suffixe_table = [i for i in range(len(genome)-1, -1, -1)]
    suffixe_table_triée = sorted(suffixe_table, key=lambda i: genome[i:])
    return suffixe_table_triée

def est_present_dicho(suffix_table, kmer, genome):
    """
    Recherche dichotomique pour trouver toutes les positions où un k-mer est présent dans le génome en utilisant la table des suffixes.
    
    :param suffix_table: table des suffixes (liste d'indices triés)
    :param kmer: le k-mer à rechercher
    :param genome: la séquence génomique
    :return: liste des positions de départ du k-mer dans le génome
    >>> genome = "ACGTACGT"
    >>> sa = suffixe_table(genome)
    >>> est_present_dicho(sa, "ACG", genome)
    [0, 4]
    >>> est_present_dicho(sa, "TTT", genome)
    []
    """
    left, right = 0, len(suffix_table) - 1
    found = False
    while left <= right:
        mid = (left + right) // 2
        suffix_start = suffix_table[mid]
        suffix = genome[suffix_start:]
        if suffix.startswith(kmer):
            found = True
            break
        elif suffix < kmer:
            left = mid + 1
        else:
            right = mid - 1
    
    if not found:
        return []
    
    left_bound = mid
    while left_bound > 0 and genome[suffix_table[left_bound - 1]:].startswith(kmer):
        left_bound -= 1
    
    right_bound = mid
    while right_bound < len(suffix_table) - 1 and genome[suffix_table[right_bound + 1]:].startswith(kmer):
        right_bound += 1
    
    positions = sorted([suffix_table[i] for i in range(left_bound, right_bound + 1)])
    return positions

# genome = read_genome('GCF_000862245.1_ViralProj15330_genomic.fna')
# table_suffixes = suffixe_table(genome)
# print(table_suffixes)
# pos = est_present_dicho(table_suffixes, "GAGTGGGTTGTTCCCACTCACTCCA", genome)
# print(pos)
# print(genome[0:30])  # Affiche une portion du génome pour vérifier la présence du k-mer

# def kmers(read: str, k: int) -> set:
#     """
#     Génère un ensemble de tous les k-mers présents depuis le read souhaité.
    
#     :param read: la séquence du read à partir de laquelle générer les k-mers
#     :param k: la longueur des k-mers
#     :return: un ensemble de k-mers
#     """
#     return set([read[i:i+k] for i in range(0, (len(read) - k + 1), 5)])

# print(kmers("TTAAAACTGGGAGTGGGTTGTTCCCACTCACTCCACCCAT", 25))
