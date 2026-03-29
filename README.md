# Projet 2 : à la recherche de la petite bête : Alignement avec une stratégie Seed-and-Extend 

## Description

Ce projet implémente une méthode d’alignement de reads sur un génome basée sur une approche **seed-and-extend**.

Le pipeline repose sur :

* une **table des suffixes** pour la phase de *seeding*
* un alignement local via la bibliothèque **Parasail** pour la phase d’*extension*

L’objectif est de retrouver la position d’origine de reads dans un génome donné.

---

## Organisation du projet

### `Suffixe_table.py`

* `read_genome(file_path)`

  * Lit un fichier FASTA et retourne la séquence du génome.
* `suffixe_table(genome)`

  * Construit la table des suffixes (liste des indices triés).
* `est_present_dicho(suffix_table, kmer, genome)`

  * Recherche les positions d’un k-mer dans le génome via une recherche dichotomique.

---

### `DecoupeKmer.py`

* `extraire_kmers_sequence(sequence, k, step)`

  * Génère les k-mers d’une séquence.
* `extraire_kmers_read_biopython(...)`

  * Extrait les k-mers d’un read depuis un FASTQ.gz.
* `extraire_kmers_fichier(...)`

  * Version générique pour FASTA/FASTQ.

---

### `Alignement.py`

* `reverse_complement(sequence)`

  * Calcule le complément inverse d’une séquence ADN.
* `align_read(read, genome, pos, id_threshold)`

  * Effectue l’alignement d’un read sur une région du génome.
  * Utilise la bibliothèque **Parasail** (`sg_trace`).
  * Teste les deux brins (+ et −).
  * Retourne :

    * statut d’alignement
    * position approximative
    * score
    * brin

---

### `Donnees_sim.py`

* `generate_reads(genome, ...)`

  * Génère des reads simulés à partir du génome.
  * Possibilité de générer sur les deux brins.

---

### `Preprocess.py`

* `filtrer_fastq(entree, sortie, qualite_min)`

  * Filtre les reads d’un fichier FASTQ.gz selon leur qualité moyenne (score Phred)
  * Écrit les reads filtrés dans un nouveau fichier FASTQ.gz

* `parse()`

  * Parse les arguments en ligne de commande pour le prétraitement

#### Utilisation du script de prétraitement

```bash
python Preprocess.py -r reads.fastq.gz -o reads_filtered.fastq.gz -q 20
```

#### Paramètres

* `-r / --reads` : fichier FASTQ.gz d’entrée (obligatoire)
* `-o / --out` : fichier FASTQ.gz de sortie (obligatoire)
* `-q / --qualite` : seuil de qualité minimum (défaut : 20)


---

### `main.py`

Contient deux modes d’exécution :

#### `main_test`

* Génère des reads simulés à partir du génome
* Permet de valider le pipeline dans un cas contrôlé

#### `main_reel`

* Utilise des reads issus d’un fichier FASTQ
* Applique le pipeline complet sur des données réelles

---

## Méthode

### 1. Seeding

* Les reads sont découpés en k-mers
* Chaque k-mer est recherché dans le génome via la table des suffixes
* Les positions correspondantes servent de points de départ

### 2. Extend

* Une région du génome est extraite autour de chaque position candidate
* Un alignement local est effectué avec **Parasail**
* Les deux orientations du read sont testées

### 3. Filtrage

* Un alignement est conservé si son identité dépasse un seuil
* Le meilleur alignement (score maximal) est retenu

---

## Installation

### Prérequis

* Python ≥ 3.8
* Bibliothèques :

  * parasail
  * biopython
  * psutil

### Installation des dépendances

```bash
pip install parasail biopython psutil
```

---

## Tests

Des tests simples sont intégrés sous forme de **doctests** directement dans les fichiers Python.

Ces tests permettent de vérifier rapidement le bon fonctionnement des fonctions principales (table des suffixes, extraction de k-mers, alignement, etc.).

### Exécution des doctests

Depuis le terminal, lancer :

```bash
python -m doctest -v Suffixe_table.py
python -m doctest -v Alignement.py
python -m doctest -v DecoupeKmer.py
```

### Objectif

Les doctests couvrent :

* la construction et l’utilisation de la table des suffixes
* l’extraction de k-mers
* les fonctions utilitaires (complément inverse, génération de reads)

Ces tests sont volontairement simples et servent principalement à valider le comportement attendu des fonctions de base.

---

## Utilisation

### Commande générale

```bash
python main.py -g genome.fasta -o output.tsv -k 12
```

### Avec reads simulés (mode test)

```bash
python main.py -g genome.fasta -o output.tsv -k 12 -i 80
```

### Avec reads réels

```bash
python main.py -g genome.fasta -r reads.fastq.gz -o output.tsv -k 12 -i 80
```

---

## Paramètres

* `-g / --genome` : fichier FASTA du génome (obligatoire)
* `-r / --reads` : fichier FASTQ des reads (optionnel)
* `-o / --out` : fichier de sortie
* `-k` : taille des k-mers (défaut : 20)
* `-i` : pourcentage d'identité (défaut : 70)

---

## Format de sortie

## Format de sortie

Le fichier `.tsv` contient pour chaque read aligné :

```
read_id    position    score    strand    sequence
```

Un résumé est ajouté à la fin du fichier :

```
# RÉSUMÉ
Total reads analysés: <nombre>
Reads alignés: <nombre>
Taux d'alignement: <pourcentage>
Reads alignés sur brin direct (+): <nombre>
Reads alignés sur brin inverse (-): <nombre>
Temps d'exécution: <secondes>
Utilisation mémoire: <MB>
Paramètres: k=<valeur>, step=5, max_reads=<valeur>
```


---