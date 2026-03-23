def read_multi_fasta(file_path: str) -> dict:
    """
    Lit un fichier multi-FASTA et retourne un dictionnaire {id_sequence: sequence}.
    
    :param file_path: chemin vers le fichier FASTA
    :return: dictionnaire contenant les séquences
    """
    sequences = {}
    seq_id = None
    seq_lines = []
    seq_test=[]

    with open(file_path, 'r') as fasta:
        for line in fasta:
            line = line.strip()
            if not line:
                continue  

            if line.startswith(">"):
                if seq_id:
                    sequences[seq_id] = ''.join(seq_lines)
                seq_id = line[1:].strip()  
                seq_lines = []
            elif line.startswith("(") or line.startswith("."):
                seq_lines.append(line)
            else:
                seq_lines.append(line)
        

        if seq_id:
            sequences[seq_id] = ''.join(seq_lines)
                

    return sequences
