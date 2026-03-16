import parasail
def align_read(read, genome, pos):
    window = len(read)
    region = genome[pos-window:pos+window]
    result = parasail.sg_trace(read, region, 10, 1, parasail.dnafull)
    # aligned sequences
    aligned_query = result.traceback.query
    aligned_ref   = result.traceback.ref

    # compute % identity
    matches = sum(q == r for q, r in zip(aligned_query, aligned_ref))
    identity_pct = matches / len(aligned_query) * 100
    if identity_pct < 90:
        print(f"Low identity: {identity_pct:.2f}%")
    else:
        print(f"Score: {result.score}")
        print(f"Percentage identity: {identity_pct:.2f}%")
        print(f"Alignment:\n{aligned_query}\n{aligned_ref}")

