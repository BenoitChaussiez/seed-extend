import parasail

def reverse_complement(sequence):
    """Calcule le complément inverse d'une séquence ADN"""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(complement.get(base, 'N') for base in reversed(sequence))


def align_read(read, genome, pos):
    """
    Aligne un read sur le génome à partir d'une position approximative.
    Teste les deux brins (direct et complément inverse).
    
    Returns:
        tuple: (aligné, position, score)
    """
    window = len(read) + 50
    start = max(0, pos - window)
    end = min(len(genome), pos + window)
    region = genome[start:end]
    
    best_align = False
    best_pos = None
    best_score = -1
    
    # Tester les deux brins
    for read_seq in [read, reverse_complement(read)]:
        result = parasail.sg_trace(read_seq, region, 5, 1, parasail.dnafull)
        
        # Calculer l'identité
        aligned_query = result.traceback.query
        aligned_ref = result.traceback.ref
        
        matches = 0
        total = 0
        for q, r in zip(aligned_query, aligned_ref):
            if q != '-' and r != '-':
                total += 1
                if q == r:
                    matches += 1
        
        identity = (matches / total * 100) if total > 0 else 0
        
        if identity >= 80:
            # Calculer la position
            if hasattr(result, 'ref_begin'):
                approx_position = start + result.ref_begin
            else:
                ref_start = 0
                for q, r in zip(aligned_query, aligned_ref):
                    if r != '-':
                        break
                    ref_start += 1
                approx_position = start + ref_start
            
            if result.score > best_score:
                best_align = True
                best_pos = approx_position
                best_score = result.score
    
    if best_align:
        return True, best_pos, best_score
    else:
        return False, None, 0