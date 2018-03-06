#!/usr/bin/env python3

"""
> plot.per-strain.abs_values.py <

Uses matplotlib to plot a linear regressions for CC7, H2 and RS.
"""
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import statsmodels.api as sm

# data is read (and subsequently plotted) in a loop--therefore define seaborn
# stuff first
sns.set_style('white')
sns.set_style('ticks')
fig, ax = plt.subplots(1, 3, figsize=(13, 4), sharey=True, sharex=True)
colours = {'CC7 25': '#9acd32', 'CC7 32' : '#548b54',
           'H2 25': '#009acd', 'H2 32': '#000080',
           'RS 25': '#e9967a', 'RS 32': '#8b3a3a'}

data = pd.read_table('fold_changes.tsv')
for n, strain in enumerate(['CC7', 'H2', 'RS']):
    # read data from both xx-25.tsv and xx-32.tsv files
    sub_data = data.filter(regex=strain)
    sub_data = sub_data.dropna()
    
    # calculate coefficient of determination (R squared) with statsmodels
    # model as first-order: y = mx + c
    y = sub_data['{}_prot'.format(strain)]
    x = sub_data['{}_tx'.format(strain)]
    X = sm.add_constant(x)
    model = sm.OLS(y, X)
    results = model.fit()
    rsq_1 = results.rsquared
    rmaxp_1 = '{:.0e}'.format(max([x for x in results.pvalues]))
    if rmaxp_1 == '0e+00': rmaxp_1 = '1e-300'
    ylab_1 = results.params[0] + 2 * results.params[1] + 0.1
    
    # subtle but important: 32 is darker in colour than 25. plot 32 first, then
    # overlay 25 on top of it
    sns.regplot(x='{}_tx'.format(strain), y='{}_prot'.format(strain),
                data=sub_data, ci=None, 
                color=colours['{} 25'.format(strain)], label='',
                scatter_kws={'s': 5, 'alpha': 0.15},
                line_kws={'lw': 2, 'color': colours['{} 32'.format(strain)],
                          'alpha': 0.8},
                ax=ax[n])
    
    ax[n].text(2, -1.2,
               'r^2 = {}\np < {}'.format(round(rsq_1, 3), rmaxp_1),
               ha='right', va='bottom')
    
    ax[n].set_title(strain)
    
    ax[n].set_xlim(-2, 2)
    ax[n].set_ylim(-1.2, 1.2)
    
    ax[n].set_xlabel('')
    ax[n].set_ylabel('')

ax[0].set_ylabel('Delta protein abundance (log2 FC)')
ax[1].set_xlabel('Delta transcript abundance (log2 FC)')

sns.despine(offset=5, trim=True)

# save figure
fig = plt.gcf()

# without bbox_inches, the saved figure has truncated axes.
output_filename = 'per-strain.fc_values.pdf'
fig.savefig(output_filename, bbox_inches='tight')
