from Bio import SeqIO
from Bio.Align import substitution_matrices
import random

# Cargar BLOSUM62 moderna
blosum62 = substitution_matrices.load("BLOSUM62")

def blosum_score(seq1, seq2):
    score = 0
    for a, b in zip(seq1, seq2):
        if a == '-' or b == '-':
            continue
        try:
            score += blosum62[(a, b)]
        except KeyError:
            try:
                score += blosum62[(b, a)]
            except KeyError:
                score -= 1
    return score

def align_sequences(seqs):
    max_len = max(len(s) for s in seqs)
    return [s.ljust(max_len, "-") for s in seqs]

def strong_mutation(seq, prob=0.3):
    mutated = ""
    for c in seq:
        if random.random() < prob:
            mutated += "-"  # inserta gap
        mutated += c
    return mutated

def average_score(aligned):
    total = 0
    count = 0
    for i in range(len(aligned)):
        for j in range(i + 1, len(aligned)):
            total += blosum_score(aligned[i], aligned[j])
            count += 1
    return total / count, total

def run_bfoa():
    # Semilla diferente por corrida para variación
    random.seed()

    # Cargar secuencias desde multifasta
    records = list(SeqIO.parse("multifasta.fasta", "fasta"))
    original = [str(r.seq) for r in records]

    num_bacterias = 10
    iteraciones = 50
    best_fitness = float('-inf')
    best_score = 0

    for _ in range(iteraciones):
        for _ in range(num_bacterias):
            mutated = [strong_mutation(s, prob=random.uniform(0.2, 0.4)) for s in original]
            aligned = align_sequences(mutated)
            fitness, score = average_score(aligned)

            if fitness > best_fitness:
                best_fitness = fitness
                best_score = score

    return round(best_fitness, 2), iteraciones, round(best_score, 2)
