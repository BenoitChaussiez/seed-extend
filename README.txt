Voici le README du projet 2 de MABS

Prerequis :
Numpy
biopython
sys

Description:
Ce programme permet de trouver des reads qui sont présent dans certains génomes de virus afin de determiner qu'elle est le virus qui est en cause pour une ensemble de patient

dossier /data:
GCF_000862245.1_RHINo_genomic.fna
GCF_000864105.1_GRIPPEA_genomic.fna
GCF_000864765.1_HIV_genomic.fna
GCF_003972325.1_CORONAVIRUS_genomic.fna
SRR10971381_1_qualite20.fastq.gz
SRR10971381_2_qualite20.fastq.gz
SRR10971381_1.fastq(1).gz
SRR10971381_2.fastq(1).gz

dossier /src:
Preprocess.py
parser.py
Suffixe_table.py
DecoupeKmer.py
main.py

Etape primordial Preprocess :
Pour gagner du temps et eviter d'inserer des reads NUL qui sont avec un quality score inférieur à 20 on a décidé de faire du preprocessing afin de retenir uniquement les read ayant un score supérieur à 20.
Pour cela on ouvre le fichier et on calcule le score du READ on fait la moyenne , si le score est supérieur à 20 on réecirt le read dans un nouveau fichier sinon il n'est pas gardé.

Parser pour le preprocess:
-r fichier du reads
-o fichier de sortie
-q qualité minimal
Exemple de ligne de commande :
python3 Preprocess.py -r SRR10971381_1.fastq.gz -o SRR10971381_1_qualite20.fastq.gz 

Etape finale pour l'alignement :


Parser:
-g fichier du genome
-r fichier du reads
-o fichier de sortie 
-k taille de kmer par défaut=25

Exemple de ligne de commande :
python3 main.py -g GCF_000862245.1_RHINo_genomic.fna -r SRR10971381_1_qualite20.fastq.gz -o Fichier_sortie.tsv 


Pour utiliser le programme il faut d'abord faire du preprocessing sur les fichiers notamment sur les fichiers de READ


Ligne de Commande de preprocess



