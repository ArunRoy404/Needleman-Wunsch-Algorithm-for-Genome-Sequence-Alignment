import pandas as pd
import numpy as np
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas

codon_dict = {
    'UUU': 'phenylalanine', 'UUC': 'phenylalanine',
    'UUA': 'leucine',       'UUG': 'leucine',
    'CUU': 'leucine',       'CUC': 'leucine',
    'CUA': 'leucine',       'CUG': 'leucine',
    'AUU': 'isoleucine',    'AUC': 'isoleucine',
    'AUA': 'isoleucine',    'AUG': 'methionine',
    'GUU': 'valine',        'GUC': 'valine',
    'GUA': 'valine',        'GUG': 'valine',
    'UCU': 'serine',        'UCC': 'serine',
    'UCA': 'serine',        'UCG': 'serine',
    'CCU': 'proline',       'CCC': 'proline',
    'CCA': 'proline',       'CCG': 'proline',
    'ACU': 'threonine',     'ACC': 'threonine',
    'ACA': 'threonine',     'ACG': 'threonine',
    'GCU': 'alanine',       'GCC': 'alanine',
    'GCA': 'alanine',       'GCG': 'alanine',
    'UAU': 'tyrosine',      'UAC': 'tyrosine',
    'UAA': 'stop',          'UAG': 'stop',
    'CAU': 'histidine',     'CAC': 'histidine',
    'CAA': 'glutamine',     'CAG': 'glutamine',
    'AAU': 'asparagine',    'AAC': 'asparagine',
    'AAA': 'lysine',        'AAG': 'lysine',
    'GAU': 'aspartic acid', 'GAC': 'aspartic acid',
    'GAA': 'glutamic acid', 'GAG': 'glutamic acid',
    'UGU': 'cysteine',      'UGC': 'cysteine',
    'UGA': 'stop',          'UGG': 'tryptophan',
    'CGU': 'arginine',      'CGC': 'arginine',
    'CGA': 'arginine',      'CGG': 'arginine',
    'AGU': 'serine',        'AGC': 'serine',
    'AGA': 'arginine',      'AGG': 'arginine',
    'GGU': 'glycine',       'GGC': 'glycine',
    'GGA': 'glycine',       'GGG': 'glycine'
}


def needleman_wunsch(seq1, seq2, match_score=5, mismatch_penalty=-3, gap_penalty=-3):
    n, m = len(seq1), len(seq2)

    score_matrix = np.zeros((n + 1, m + 1), dtype=int)
    for i in range(1, n + 1):
        score_matrix[i][0] = score_matrix[i - 1][0] + gap_penalty
    for j in range(1, m + 1):
        score_matrix[0][j] = score_matrix[0][j - 1] + gap_penalty

    
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = score_matrix[i - 1][j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_penalty)
            delete = score_matrix[i - 1][j] + gap_penalty
            insert = score_matrix[i][j - 1] + gap_penalty
            score_matrix[i][j] = max(match, delete, insert)


    aligned_seq1 = []
    aligned_seq2 = []
    i, j = n, m
    while i > 0 and j > 0:
        score = score_matrix[i][j]
        diagonal = score_matrix[i - 1][j - 1]
        up = score_matrix[i - 1][j]
        left = score_matrix[i][j - 1]

        if score == diagonal + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_penalty):
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif score == up + gap_penalty:
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append('-')
            i -= 1
        else:
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j - 1])
            j -= 1

    while i > 0:
        aligned_seq1.append(seq1[i - 1])
        aligned_seq2.append('-')
        i -= 1

    while j > 0:
        aligned_seq1.append('-')
        aligned_seq2.append(seq2[j - 1])
        j -= 1

    return ''.join(reversed(aligned_seq1)), ''.join(reversed(aligned_seq2))


    
def annotate_mutations(aligned_seq1, aligned_seq2):
    mutations = []
    positions = []
    deletion_done = 0
    similarity = 0
    for i, (base1, base2) in enumerate(zip(aligned_seq1, aligned_seq2), start=1):
        if base1 != base2:
            if base1 == '-':
                mutations.append(f'Insertion    at position {i}: {base2}')
                positions.append(i-deletion_done)
            elif base2 == '-':
                mutations.append(f'Deletion     at position {i}: {base1}')
                deletion_done+=1
            else:
                mutations.append(f'Substitution at position {i}: {base1} -> {base2}')
                positions.append(i-deletion_done)
        else:
            similarity += 1
    return mutations, positions, similarity



def get_codons(sequence):
    codons =[]
    for i in range(0, len(sequence),3):
        
        if len(sequence[i:i+3]) == 3:
            codons.append(sequence[i:i+3])
    return codons


def save_to_pdf(aligned_seq1, aligned_seq2, mutations, mutated_codons, file_path):
    c = canvas.Canvas(file_path, pagesize=letter)
    text_object = c.beginText(40, 750)

    text_object.setFont("Courier-Bold", 10)
    text_object.setFillColor(colors.black)
    text_object.textLine("Aligned Sequences:")

    # Display aligned sequences
    for i in range(0, len(aligned_seq1), 80):
        segment1 = aligned_seq1[i:i + 80]
        segment2 = aligned_seq2[i:i + 80]
        for char1, char2 in zip(segment1, segment2):
            if char2 == '-':
                text_object.setFillColor(colors.red)
            else:
                text_object.setFillColor(colors.black)
            text_object.textOut(char1)
        text_object.textLine()

        for char1, char2 in zip(segment1, segment2):
            if char1 == char2:
                text_object.setFillColor(colors.black)
            elif char1 == '-':
                text_object.setFillColor(colors.blue)
            elif char1 != char2 and char2 != "-":
                text_object.setFillColor(colors.green)
            text_object.textOut(char2)
        text_object.textLine()
        text_object.textLine()
    
    text_object.textLine(f"Mutation Nucleotide: {len(mutations)}")    
    text_object.textLine("\n\n\nMutations:")
    for mutation in mutations:
        text_object.textLine(mutation)    
    text_object.textLine()
    text_object.textLine()

   
    text_object.textLine(f"\n\n\nMutated Amino Acids: {len(mutated_codons)}")
    for amino_acid in mutated_codons:
        mutation_line = f"{amino_acid} -> {codon_dict[amino_acid]}"
        text_object.textLine(mutation_line)        
        

    c.drawText(text_object)
    c.save()
    print(f"PDF saved as {file_path}.")



# def main():
#     original_file = "original_sequence.txt"
#     mutated_file = "mutated_sequence.txt"
#     original_sequence = read_file(original_file)
#     mutated_sequence = read_file(mutated_file)

#     original_sequence_codons = get_codons(original_sequence)
#     mutated_sequence_codons = get_codons(mutated_sequence)

#     aligned_seq1, aligned_seq2 = needleman_wunsch(original_sequence, mutated_sequence)
#     mutations, positions = annotate_mutations(aligned_seq1, aligned_seq2)

#     mutated_codon_positions = list(set((i-1)//3 for i in positions))
#     mutated_codons = list(mutated_sequence_codons[i] for i in mutated_codon_positions)
    
#     save_to_pdf(aligned_seq1, aligned_seq2, mutations, mutated_codons, "aligned_sequences.pdf")



def findMutation(seq1, seq2, i=-1):

    original_sequence = seq1
    mutated_sequence = seq2

    original_sequence_codons = get_codons(original_sequence)
    mutated_sequence_codons = get_codons(mutated_sequence)

    aligned_seq1, aligned_seq2 = needleman_wunsch(original_sequence, mutated_sequence)
    aligned_length=len(aligned_seq1)
    mutations, positions, similarity = annotate_mutations(aligned_seq1, aligned_seq2)

    mutated_codon_positions = list(set((i-1)//3 for i in positions))
    mutated_codons = list(mutated_sequence_codons[i] for i in mutated_codon_positions)

    #test in pdf
    if i!=-1:
        save_to_pdf(aligned_seq1, aligned_seq2, mutations, mutated_codons, "aligned_sequences.pdf")
        for i in mutations:
            print(i)
        print()
        for i in mutated_codons:
            print(i)
    mutated_amino_acids = []
    for i in mutated_codons:
        mutated_amino_acids.append(codon_dict[i])

    return len(mutations), len(mutated_codons), aligned_length, similarity, mutated_amino_acids

    

def main():
    index = -1
    data = pd.read_csv('geneTableCsv.csv')

    for i in range(len(data['Patient_Id'])):
        seq1 = data['Original_Sequence'][i]
        seq2 = data['Mutated_Sequence'][i]
        if index == i:
            mutations, mutated_codons, alignment_length, similarity, mutated_amino_acids = findMutation(seq1,seq2,index)
        else:
            mutations, mutated_codons, alignment_length, similarity, mutated_amino_acids  = findMutation(seq1,seq2)
        print(mutations, mutated_codons,alignment_length,similarity)
        data['Calculated_Nucleotide'][i] = int(mutations)
        data['Calculated_Amino_Acid'][i] = int(mutated_codons)
        data['alignment_length'][i]      = alignment_length
        data['homology'][i]              = similarity/alignment_length
        data['Mutated_Amino_Acids'][i]   = mutated_amino_acids

        
    saveData = pd.DataFrame(data)
    saveData.to_csv("geneTableCsv2.csv", index=False)
    print(data)


if __name__ == "__main__":
    main()


