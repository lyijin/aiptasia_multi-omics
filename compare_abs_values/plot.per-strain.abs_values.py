#!/usr/bin/env python3

"""
> plot.per-strain.abs_values.py <

Uses matplotlib to plot a linear regressions for CC7, H2 and RS.
"""
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import statsmodels.api as sm

def melt_table_pairwise(pd_dataframe):
    """
    The input dataframes are in the pattern:
        tx_X; tx_Y; tx_Z; prot_X; prot_Y; prot_Z
        
    This def melts the dataframe into
        tx_X; prot_X
        ...
        tx_Y; prot_Y
        ...
        tx_Z; prot_Z
        ...
    """
    pd_dataframe = pd_dataframe.drop('gene_id', axis=1)
    skip_cols = int(len(pd_dataframe.columns) / 2)
    
    temp = pd.DataFrame({'tx': [], 'prot': []})
    for n in range(0, skip_cols):
        # num_cols is definitely even. if it isn't, there's something wrong
        # with the upstream code that tabulated the values
        temp = pd.concat(
            [temp,
             pd.DataFrame({'tx': pd_dataframe.iloc[:, n],
                           'prot': pd_dataframe.iloc[:, n+skip_cols]})],
            ignore_index=True)
    
    temp = temp.dropna()
    return temp

# data is read (and subsequently plotted) in a loop--therefore define seaborn
# stuff first
sns.set_style('white')
sns.set_style('ticks')
fig, ax = plt.subplots(2, 3, figsize=(13, 7), sharey=True, sharex=True)
colours = {'CC7 25': '#9acd32', 'CC7 32' : '#548b54',
           'H2 25': '#009acd', 'H2 32': '#000080',
           'RS 25': '#e9967a', 'RS 32': '#8b3a3a'}

for n, strain in enumerate(['CC7', 'H2', 'RS']):
    # read data from both xx-25.tsv and xx-32.tsv files
    data_25 = pd.read_table('{}-25.tsv'.format(strain))
    data_25 = melt_table_pairwise(data_25)
    
    data_32 = pd.read_table('{}-32.tsv'.format(strain))
    data_32 = melt_table_pairwise(data_32)
    
    # calculate coefficient of determination (R squared) with statsmodels
    # model as first-order: y = mx + c
    y = data_25['prot']
    x = data_25['tx']
    X = sm.add_constant(x)
    model = sm.OLS(y, X)
    results = model.fit()
    rsq_1 = results.rsquared
    rmaxp_1 = '{:.0e}'.format(max([x for x in results.pvalues]))
    if rmaxp_1 == '0e+00': rmaxp_1 = '1e-300'
    
    y = data_32['prot']
    x = data_32['tx']
    X = sm.add_constant(x)
    model = sm.OLS(y, X)
    results = model.fit()
    rsq_2 = results.rsquared
    rmaxp_2 = '{:.0e}'.format(max([x for x in results.pvalues]))
    if rmaxp_2 == '0e+00': rmaxp_2 = '1e-300'
    
    # subtle but important: 32 is darker in colour than 25. plot 32 first, then
    # overlay 25 on top of it
    sns.regplot(x='tx', y='prot', data=data_25, ci=None, 
                color=colours['{} 25'.format(strain)], label='', marker='^',
                scatter_kws={'s': 5, 'alpha': 0.15},
                line_kws={'lw': 2, 'color': colours['{} 25'.format(strain)],
                          'alpha': 0.8},
                ax=ax[0][n])
    sns.regplot(x='tx', y='prot', data=data_32, ci=None, 
                color=colours['{} 32'.format(strain)], label='', marker='s',
                scatter_kws={'s': 5, 'alpha': 0.15},
                line_kws={'lw': 2, 'color': colours['{} 32'.format(strain)],
                          'alpha': 0.8},
                ax=ax[1][n])
    
    ax[0][n].text(14, 0,
                  'r^2 = {}\np < {}'.format(round(rsq_1, 3), rmaxp_1),
                  ha='right', va='bottom')
    ax[1][n].text(14, 0,
                  'r^2 = {}\np < {}'.format(round(rsq_2, 3), rmaxp_2),
                  ha='right', va='bottom')
    
    ax[0][n].set_title(strain)
    
    ax[0][n].set_xlim(0, 14)
    ax[0][n].set_ylim(0, 14)
    ax[1][n].set_xlim(0, 14)
    ax[1][n].set_ylim(0, 14)
    
    ax[0][n].set_xlabel('')
    ax[0][n].set_ylabel('')
    ax[1][n].set_xlabel('')
    ax[1][n].set_ylabel('')

ax[0][0].set_ylabel('Protein abundance (log2 SCPM)')
ax[1][1].set_xlabel('Transcript abundance (log2 TPM)')

sns.despine(offset=5, trim=True)

# save figure
fig = plt.gcf()

# without bbox_inches, the saved figure has truncated axes.
output_filename = 'per-strain.abs_values.pdf'
fig.savefig(output_filename, bbox_inches='tight')
