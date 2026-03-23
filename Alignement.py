import parasail
def align_read(read, genome, pos):
    window = len(read) + 50
    start = max(0, pos - window)
    end = min(len(genome), pos + window)

    region = genome[start:end]

    result = parasail.sg_trace(read, region, 10, 1, parasail.dnafull)
    score = result.score

    aligned_query = result.traceback.query
    aligned_ref   = result.traceback.ref

    matches = sum(q == r for q, r in zip(aligned_query, aligned_ref))
    identity_pct = matches / len(aligned_query) * 100

    if identity_pct >= 90:
        approx_position = start + result.ref_begin
        return True, approx_position, score
    else:
        return False, None, score