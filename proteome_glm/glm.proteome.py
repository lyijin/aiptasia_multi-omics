#!/usr/bin/env python3

"""
> glm.proteome.py <

Perform a GLM analysis on a per-gene basis to detect significant changes in
proteome that can be attributed to strain (CC7, SSB01, H2) or 
temperature (25, 32), or through the interaction of both variables!
"""
import collections

import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

import correct_p_values

# read data
data = pd.read_table('sheet_7_standalone.tsv', index_col=0)
data = data.drop('Majority protein IDs', axis=1)

# for now, remove rows that contain A SINGLE NA (i.e. all remaining rows
# have protein values across ALL replicates
data = data.dropna(how='all')

# convert each row into a dataframe that can be used for subsequent GLM
pval_table = pd.DataFrame()
for gene, row in data.iterrows():
    data_per_gene = pd.DataFrame({'strain': [], 'temperature': [], 
                                  'prot_expr': []})
    for r in row.iteritems():
        strain = r[0].split(' ')[0]
        temp = r[0].split(' ')[1]
        
        # if prot_expr is np.nan, skip to the next item
        if pd.isnull(r[1]): continue
        
        prot_expr = float(r[1])
        temp = pd.DataFrame({'strain': strain, 'temperature': temp,
                             'prot_expr': prot_expr}, index=[0])
        data_per_gene = pd.concat([data_per_gene, temp], ignore_index=True)
        
    # remove any combination of strain/temp that only occur once/twice
    data_per_gene['sttemp'] = data_per_gene['strain'] + data_per_gene['temperature']
    counter = collections.Counter(data_per_gene['sttemp'])
    if len(counter) < 6 or min(counter.values()) < 2: continue
    print (len(counter), counter.values())
    
    
    model = smf.glm('prot_expr ~ strain * temperature', data=data_per_gene,
                    family=sm.families.Gaussian())
    results = model.fit()
    
    # hardcode parsing of results' p values
    pvals = results.pvalues.values
    
    temp = pd.DataFrame(results.pvalues).transpose()
    temp.index = [gene]
    
    # temp = pd.DataFrame({'strain_h2': pvals[1],
                         # 'strain_rs': pvals[2],
                         # 'temp': pvals[3],
                         # 'strain_h2:temp': pvals[4],
                         # 'strain_rs:temp': pvals[5]}, index=[gene])
    
    pval_table = pd.concat([pval_table, temp])

# correct for multiple testing
for cols in pval_table.columns:
    x = pd.Series(correct_p_values.correct_p_values(pval_table[cols]),
                  name='corr.{}'.format(cols))
    pval_table = pd.concat([pval_table, x], axis=1)

# save output as file
pval_table.to_csv('glm_output.tsv', sep='\t')
