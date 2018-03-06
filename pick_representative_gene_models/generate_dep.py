#!/usr/bin/env python3

"""
> generate_dep.py <

Reads in the protein fold changes .tsv file, then pick out proteins that 
are up/downexpressed.

NOTE: for GO term analysis, the LEFTMOST & VALID gene is picked for ambiguous 
protein groups--this avoids biasing the universe/DEP list with multiple
occurrences of genes with the same function.
"""
import csv
import math
import re

import natural_sort

# read permissible gene IDs: genes that are not isoforms/duplicates of others
valid_genes = []
for line in open('../ultimate_universe.txt'):
    if not line: continue    
    valid_genes.append(line.strip())

# read protein data
#   prot_data[gene][sample] = prot_value (no transformation required)
prot_data = {}

tsv_reader = csv.reader(open('prot.fold_changes.tsv'), delimiter='\t')
strains = next(tsv_reader)[1:]
for row in tsv_reader:
    if not row: continue
    
    # change of plans: if gene_id contains multiple AIPGENEs chained together,
    # then assign the same protein value to the leftmost valid gene).
    gene_ids = re.findall('AIPGENE\d+', row[0])
    gene_ids = [x for x in gene_ids if x in valid_genes]
    if not gene_ids: continue
    
    log2_tpms = [float(x) if x else 0 for x in row[1:]]
    prot_data[gene_ids[0]] = {x:y for x, y in zip(strains, log2_tpms)}

# print stuff out
for s in strains:
    up   = [x for x in prot_data if prot_data[x][s] > math.log(1.2, 2)]
    down = [x for x in prot_data if prot_data[x][s] < -math.log(1.2, 2)]
    combined = up + down
    universe = [x for x in prot_data if prot_data[x][s] != 0]
    
    with open('{}_prot_diff.txt'.format(s.lower()), 'w') as f:
        for x in natural_sort.natural_sort(combined):
            print (x, file=f)
    
    with open('{}_prot_up.txt'.format(s.lower()), 'w') as f:
        for x in natural_sort.natural_sort(up):
            print (x, file=f)
    
    with open('{}_prot_down.txt'.format(s.lower()), 'w') as f:
        for x in natural_sort.natural_sort(down):
            print (x, file=f)
    
    with open('{}_prot_universe.txt'.format(s.lower()), 'w') as f:
        for x in natural_sort.natural_sort(universe):
            print (x, file=f)
