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
    
    return ''.join(genome)
genome = read_genome("/Users/paullemonnier/Desktop/Master_MISO/Master 1/S2/MABS/Projet_seed-extend/GCF_000862245.1_ViralProj15330_genomic.fna")

suffixe_table = [i for i in range(len(genome)-1, -1, -1)]
print(suffixe_table)
def tri_suffixes(genome: str) -> list:
    """
    Trie les suffixes d'une séquence génomique et retourne une table de suffixes.
    
    :param genome: séquence génomique
    :return: table de suffixes triée
    """
    pass