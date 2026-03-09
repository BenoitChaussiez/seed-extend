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
    genome.append("$")  # Ajouter un caractère de fin de séquence
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


