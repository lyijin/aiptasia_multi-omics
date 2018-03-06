#!/usr/bin/env python3

"""
> compute_fc.py <

Compile the transcript and protein fold changes.
"""
import csv
import math

import natural_sort

strains =  ['CC7', 'H2', 'RS', 'SSB01']

# read permissible gene IDs: genes that are not isoforms/duplicates of others
valid_genes = []
for line in open('../ultimate_universe.txt'):
    if not line: continue    
    valid_genes.append(line.strip())

# read protein data
prot_data = {}
tsv_reader = csv.reader(open('prot.fold_changes.tsv'), delimiter='\t')

header_row = next(tsv_reader)[1:]
for row in tsv_reader:
    if not row: continue
    
    # change of plans: if gene_id contains multiple AIPGENEs chained together
    # (i.e. ambiguous protein groups) then ignore the entire group; not 
    # assign the same protein value to all those genes, nor pick leftmost gene.
    # this is tested AFTER removing proteins that are NOT in the universe.
    #   e.g. AIPGENE3199;AIPGENE3120 passes the filter because AIPGENE3120 is
    #        not a valid gene.
    gene_ids = row[0].split(';')
    gene_ids = [x for x in gene_ids if x in valid_genes]
    if len(gene_ids) != 1: continue
    gene_id = gene_ids[0]
    
    prot_data[gene_id] = {x:y for x, y in zip(header_row, row[1:])}

# read tx data
tx_data = {}
for s in strains:
    tsv_reader = csv.reader(
        open('{}_temp32_vs_temp25_spliced_rt.csv'.format(s.lower())))
    
    header_row = next(tsv_reader)[1:]
    for row in tsv_reader:
        if not row: continue
        
        gene_id = row[0]
        if gene_id not in valid_genes: continue
        
        if gene_id not in tx_data:
            tx_data[gene_id] = {}
        
        tx_data[gene_id][s] = row[3]

all_possible_genes = list(tx_data.keys()) + list(prot_data.keys())
all_possible_genes = natural_sort.natural_sort(set(all_possible_genes))

# print stuff out
output_header = ['gene_id'] + \
                [x + '_' + y for x in strains for y in ['tx', 'prot']]
with open('fold_changes.tsv', 'w') as f:
    print (*output_header, sep='\t', file=f)
    for g in all_possible_genes:
        output_line = [g]
        for s in strains:
            temp = ['', '']
            if g in tx_data:
                if s in tx_data[g]:
                    if tx_data[g][s] != 'NA':
                        beta = float(tx_data[g][s])
                        log2_beta = beta / math.log(2)
                        temp[0] = round(log2_beta, 5)
            
            if g in prot_data:
                if s in prot_data[g]:
                    if prot_data[g][s]:
                        temp[1] = round(float(prot_data[g][s]), 5)
            
            output_line += temp
        
        # only print lines that has at least a value
        if any(output_line[1:]):
            print (*output_line, sep='\t', file=f)
