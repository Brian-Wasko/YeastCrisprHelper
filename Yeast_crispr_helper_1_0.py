# -*- coding: utf-8 -*-
"""Yeast CRISPR Helper 1.0.ipynb

Original file is located at
    https://colab.research.google.com/drive/16WtgNXZ3trjumLkAIvDuwVZqKPtoYoYH
"""

# Wasko Yeast CRISPR Helper
# Requires Biopython !pip install biopython
!pip install biopython

import requests
import re
from Bio.Seq import Seq  # Use Biopython to handle sequence manipulation

def resolve_gene_name(yeast_gene):
    """Resolve the common or systematic yeast gene name to its Ensembl gene ID."""
    url = f"https://rest.ensembl.org/xrefs/symbol/saccharomyces_cerevisiae/{yeast_gene}?content-type=application/json"
    response = requests.get(url)
    if response.status_code == 200 and len(response.json()) > 0:
        return response.json()[0]['id']  # Return the Ensembl ID (e.g., YBR127C)
    else:
        raise ValueError("Gene not found. Please check the gene name.")

def get_gene_sequence(yeast_gene):
    """Fetch the yeast gene sequence using Ensembl API for the resolved gene ID."""
    url = f"https://rest.ensembl.org/sequence/id/{yeast_gene}?content-type=text/plain"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text.strip()  # Remove any leading/trailing whitespace
    else:
        raise ValueError("Gene sequence not found in database.")

def find_cas9_sites(gene_sequence, amino_acid_position, window=100):
    """Identify Cas9 guide RNA target sites near the amino acid position, on both forward and reverse strands."""
    # Convert amino acid position to nucleotide position
    nucleotide_position = (amino_acid_position - 1) * 3  # Convert to 0-based index and nucleotide position

    # Extract a region of 'window' nucleotides around the amino acid position
    start = max(0, nucleotide_position - window // 2)
    end = min(len(gene_sequence), nucleotide_position + window // 2)
    region = gene_sequence[start:end]

    # Forward strand: look for NGG PAM
    forward_cas9_sites = []
    for m in re.finditer(r'(?=([ACGT]{20}GG))', region):
        pos = start + m.start()
        guide_pam_seq = m.group(1)  # 20 nt guide + 'GG' PAM
        forward_cas9_sites.append((pos, guide_pam_seq, 'forward'))

    # Reverse strand: look for CCN PAM (NGG in reverse complement)
    reverse_cas9_sites = []
    for m in re.finditer(r'(?=(CC[ACGT]{20}))', region):
        pos = start + m.start()
        guide_pam_seq = m.group(1)  # 'CC' PAM + 20 nt guide
        reverse_cas9_sites.append((pos, guide_pam_seq, 'reverse'))

    # Combine and sort
    all_cas9_sites = forward_cas9_sites + reverse_cas9_sites
    cas9_sites_sorted = sorted(all_cas9_sites, key=lambda site: abs(site[0] - nucleotide_position))
    return cas9_sites_sorted

def generate_repair_template(gene_sequence, cas9_site, amino_acid_position, new_amino_acid):
    """Create a repair template with the mutation centered as much as possible, ensuring a silent mutation in the PAM."""
    from Bio.Seq import Seq

    # Convert amino acid position to nucleotide position
    mutation_position = (amino_acid_position - 1) * 3

    # Set the minimum flanking distance
    min_flanking_length = 30

    # Determine the start and end of the repair template
    cas9_cut_position = cas9_site[0] + 17  # Cas9 cuts 3 bp upstream of PAM
    pam_start = cas9_site[0]  # PAM position

    # Ensure the repair template has at least 30 nt on both sides of both mutations
    homology_start = max(0, min(cas9_cut_position, mutation_position) - min_flanking_length)
    homology_end = min(len(gene_sequence), max(cas9_cut_position, mutation_position) + min_flanking_length)

    # Extract the homology region
    homology_region = gene_sequence[homology_start:homology_end]

    # Complete codon table for substitution
    codon_table = {
        'A': ['GCT', 'GCC', 'GCA', 'GCG'], 'C': ['TGT', 'TGC'], 'D': ['GAT', 'GAC'],
        'E': ['GAA', 'GAG'], 'F': ['TTT', 'TTC'], 'G': ['GGT', 'GGC', 'GGA', 'GGG'],
        'H': ['CAT', 'CAC'], 'I': ['ATT', 'ATC', 'ATA'], 'K': ['AAA', 'AAG'],
        'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'M': ['ATG'],
        'N': ['AAT', 'AAC'], 'P': ['CCT', 'CCC', 'CCA', 'CCG'],
        'Q': ['CAA', 'CAG'], 'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
        'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'T': ['ACT', 'ACC', 'ACA', 'ACG'],
        'V': ['GTT', 'GTC', 'GTA', 'GTG'], 'W': ['TGG'], 'Y': ['TAT', 'TAC']
    }

    # Convert to a list for easy mutation
    homology_list = list(homology_region)

    # Apply the desired missense mutation
    codon_start_in_homology = mutation_position - homology_start
    new_codon = codon_table[new_amino_acid.upper()][0].lower()  # Use the first codon
    for i in range(3):
        homology_list[codon_start_in_homology + i] = new_codon[i]

    # Apply silent mutation to the PAM site
    homology_list = mutate_pam(homology_list, cas9_site[1], cas9_site[2], pam_start, homology_start)

    # Convert the list back to a string
    mutated_homology = ''.join(homology_list)

    # Generate reverse complement while preserving case
    repair_seq = Seq(mutated_homology)
    reverse_complement = repair_seq.reverse_complement()

    return mutated_homology, str(reverse_complement)

def mutate_pam(homology_list, guide_pam_seq, strand, pam_start, homology_start):
    """Ensure that the mutation only targets the PAM (GG for forward or CC for reverse strand)."""
    from Bio.Seq import Seq

    # Adjust PAM site for forward or reverse strand
    if strand == 'forward':
        pam_site = guide_pam_seq[-3:]  # The NGG PAM sequence on the forward strand
        pam_position_in_homology = pam_start - homology_start
        first_G_idx = pam_position_in_homology + 1
        second_G_idx = pam_position_in_homology + 2
    else:
        # Reverse strand PAM (CCN corresponds to NGG)
        pam_site = str(Seq(guide_pam_seq).reverse_complement())[:3]
        pam_position_in_homology = pam_start - homology_start
        first_G_idx = pam_position_in_homology + 1  # First "C" of CCN (mutating C)
        second_G_idx = pam_position_in_homology + 2  # Second "C" of CCN (mutating C)

    # Ensure we target the GG/CC of the PAM site only
    mutated = False
    if strand == 'forward':
        # Mutate G's of NGG
        if homology_list[first_G_idx].upper() == 'G':
            homology_list[first_G_idx] = 'c'
            mutated = True
        elif homology_list[second_G_idx].upper() == 'G':
            homology_list[second_G_idx] = 'c'
            mutated = True
    else:
        # Mutate C's of CCN (NGG in reverse)
        if homology_list[first_G_idx].upper() == 'C':
            homology_list[first_G_idx] = 't'
            mutated = True
        elif homology_list[second_G_idx].upper() == 'C':
            homology_list[second_G_idx] = 't'
            mutated = True

    if not mutated:
        print("Error: PAM silent mutation could not be applied correctly.")

    return homology_list

def main():
    # User inputs
    gene_name = input("Enter yeast gene name (common or systematic): ")
    amino_acid_number = int(input("Enter amino acid number: "))
    desired_aa = input("Enter desired amino acid substitution (single-letter code): ").upper()

    # Step 1: Resolve gene name (common or systematic)
    try:
        resolved_gene_name = resolve_gene_name(gene_name)
        print(f"Resolved gene name: {resolved_gene_name}")
    except ValueError as e:
        print(e)
        return

    # Step 2: Fetch the gene sequence
    gene_seq = get_gene_sequence(resolved_gene_name)

    # Step 3: Find Cas9 target sites near the amino acid position
    cas9_sites = find_cas9_sites(gene_seq, amino_acid_number)

    if not cas9_sites:
        print("No nearby Cas9 target sites found.")
        return

    # Step 4: Choose up to 6 closest Cas9 sites (or fewer if not enough are found)
    chosen_sites = cas9_sites[:6]

    # Step 5: Generate repair templates for each selected Cas9 site
    repair_templates = []
    for idx, site in enumerate(chosen_sites):
        try:
            repair_template, reverse_complement = generate_repair_template(gene_seq, site, amino_acid_number, desired_aa)
            repair_templates.append((site, repair_template, reverse_complement))
        except ValueError as e:
            print(f"Error for sgRNA {idx + 1}: {e}")
            continue

    # Step 6: Display all found guide RNAs and corresponding repair templates
    for idx, (site, template, rev_comp) in enumerate(repair_templates):
        if site[2] == 'reverse':
            # For reverse strand, reverse complement the guide_pam_seq to get sgRNA sequence ending with NGG
            sgRNA_seq = str(Seq(site[1]).reverse_complement())
        else:
            sgRNA_seq = site[1]  # Forward strand, guide RNA with PAM at 3'

        # Ensure sgRNA ends with NGG
        if not sgRNA_seq.endswith('GG'):
            print(f"Warning: sgRNA {idx + 1} does not end with NGG.")

        sgRNA_rev_comp = str(Seq(sgRNA_seq).reverse_complement())

        # Remove the PAM for the modified sgRNA (only for the modified sgRNA for cloning)
        modified_sgRNA_seq = sgRNA_seq[:-3]  # Remove PAM

        # Add cloning sequences to the modified sgRNA
        modified_sgRNA = 'gatc' + modified_sgRNA_seq + 'GTTTTAGAGCTAG'.lower()
        modified_sgRNA_rev_comp = 'ctagctctaaaac'.lower() + sgRNA_rev_comp[:-3]  # Also remove PAM for reverse complement

        strand_direction = "forward" if site[2] == 'forward' else "reverse"
        print(f"\nsgRNA {idx + 1}: {sgRNA_seq} (Strand: {strand_direction}, Position in gene: {site[0]})")
        if site[2] == 'reverse':
            print(f"Reverse complement of sgRNA {idx + 1}: {sgRNA_rev_comp}")
        print(f"pML104 cloning ready sgRNA {idx + 1}: {modified_sgRNA}")
        print(f"pML104 cloning ready Reverse Complement of sgRNA {idx + 1}: {modified_sgRNA_rev_comp}")
        print("     lowercase on cloning ready sgRNA sequences indicates sequence needed for cloning")
        print(f"Repair Template {idx + 1}: {template}")
        print(f"Reverse Complement of Repair Template {idx + 1}: {rev_comp}")
        print("     lowercase on repair templates indicates mutations")

    if len(chosen_sites) < 6:
        print(f"\nOnly {len(chosen_sites)} Cas9 target sites were found.")

if __name__ == "__main__":
    main()