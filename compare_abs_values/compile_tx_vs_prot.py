#!/usr/bin/env python3

"""
> compile_tx_vs_prot.py <

Based on the transcript and protein data, merge across the files to produce
something that can be plotted as a scatterplot.

Uses hardcoded filenames, as it's not really reusable anyway *shrug*.
"""
import csv
import math
import statistics

import natural_sort

def force_rounding(float_string):
    try:
        return round(float(float_string), 4)
    except:
        return ''

# read permissible gene IDs: genes that are not isoforms/duplicates of others
valid_genes = []
for line in open('../ultimate_universe.txt'):
    if not line: continue    
    valid_genes.append(line.strip())

# read transcript data
#   tx_data[gene][sample] = log2(tpm + 1)
tx_data = {}

tsv_reader = csv.reader(open('kallisto.tpms.tsv'), delimiter='\t')
tx_samples = next(tsv_reader)[1:]
for row in tsv_reader:
    if not row: continue
    
    gene_id = row[0]
    if gene_id not in valid_genes: continue
    
    log2_tpms = [math.log(float(x) + 1, 2) for x in row[1:]]
    tx_data[gene_id] = {x:y for x, y in zip(tx_samples, log2_tpms)}

# read protein data
#   prot_data[gene][sample] = prot_value (no transformation required)
prot_data = {}

tsv_reader = csv.reader(open('prot.scpms.tsv'), delimiter='\t')
prot_samples = next(tsv_reader)[1:]
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
    
    log2_scpms = [math.log(float(x) + 1, 2) if x != '0' else '' for x in row[1:]]
    prot_data[gene_id] = {x:y for x, y in zip(prot_samples, log2_scpms)}

# determine strains & temperature combinations
strains = sorted(list(set([x.split('-')[0] for x in tx_samples])))
temps = sorted(list(set([x.split('-')[1] for x in tx_samples])))

for s in strains:
    for t in temps:
        st = '-'.join([s, t])
        
        # create sublists of samples on a per-sample, per-temp basis
        tx_st = []
        prot_st = []
        for u in tx_samples:
            if st in u: tx_st.append(u)
        
        for v in prot_samples:
            if st in v: prot_st.append(v)
        
        # only consider the intersection of tx_st and prot_st: ignore samples
        # where we have transcript data but not protein!
        prot_tx_st = sorted(list(set(tx_st) & set(prot_st)))
        
        # grab list of valid protein IDs--at least one non-blank value across 
        # samples of a particular strain/temp combination
        valid_prots = []
        for gene_id in prot_data:
            prot_values = [prot_data[gene_id][x] for x in prot_tx_st]
            if any(prot_values): valid_prots.append(gene_id)
        
        prot_tx_genes = set(valid_prots) & set(tx_data.keys())
        prot_tx_genes = natural_sort.natural_sort(prot_tx_genes)
        
        # print RNA values first, then proteins
        # BONUS: calculate means for RNA and proteins, and print them as
        #        subsequent columns
        output_filename = open(st + '.tsv', 'w')
        with output_filename as f:
            header = ['gene_id'] + ['tx_' + x for x in prot_tx_st] + \
                     ['prot_' + x for x in prot_tx_st]
            print (*header, sep='\t', file=f)
            
            for gene_id in prot_tx_genes:
                gene_tx_values = [force_rounding(tx_data[gene_id][x]) for x in prot_tx_st]
                gene_prot_values = [force_rounding(prot_data[gene_id][x]) for x in prot_tx_st]
                
                print (gene_id, *gene_tx_values, *gene_prot_values, sep='\t', file=f)
