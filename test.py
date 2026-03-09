from Bio import SeqIO
import gzip

def lire_fasta(chemin):
    """Lit un fichier FASTA (compressé ou non)"""
    genome = ""
    if chemin.endswith('.gz'):
        ouverture = gzip.open(chemin, 'rt')
    else:
        ouverture = open(chemin, 'r')
    
    with ouverture as f:
        for record in SeqIO.parse(f, 'fasta'):
            genome += str(record.seq).upper()
    return genome

def lire_fastq(chemin, n_max=None):
    """Lit un fichier FASTQ compressé et retourne une liste de reads"""
    reads = []
    with gzip.open(chemin, 'rt') as f:
        for i, record in enumerate(SeqIO.parse(f, 'fastq')):
            if n_max and i >= n_max:
                break
            reads.append({
                'nom': record.id,
                'sequence': str(record.seq).upper(),
                'qualite': record.letter_annotations["phred_quality"] if hasattr(record, 'letter_annotations') else None
            })
    return reads
test = lire_fasta('GCF_000862245.1_RHINo_genomic.fna.gz')
print(test)
test= lire_fastq('SRR10971381_1.fastq(1).gz')
