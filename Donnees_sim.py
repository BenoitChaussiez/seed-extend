import random

def generate_reads(genome, read_length=150, n_reads=1000):
    reads = []
    for _ in range(n_reads):
        start = random.randint(0, len(genome) - read_length)
        read = genome[start:start+read_length]
        reads.append(read)
    return reads
