#!/usr/bin/python3

"""
run meta analysis on the data
"""

import pandas as pd
import matplotlib.pyplot as plt
import math

df = pd.read_csv('output.csv')
ks_p = df['r_ks_p']
df['mlogp'] = [-1*math.log(10**-20,10) if x==0 else -1*math.log(x,10) for x in ks_p]

fig1, ax = plt.subplots(1,1)
ax.scatter(df['n_all'],df['mlogp'])
ax.set_xlabel('total positions')
ax.set_ylabel('-log10(ks pvalue)')
plt.savefig('figures/total_vs_ksp.png')

fig2, ax = plt.subplots(1,1)
ax.scatter((df['n_pos']/(df['n_pos']+df['n_impos'])),df['mlogp'])
ax.set_xlabel('% possible (pos/(pos+impos))')
ax.set_ylabel('-log10(ks pvalue)')
plt.savefig('figures/centpos_vs_ksp.png')

fig3, ax = plt.subplots(1,1)
df.boxplot(ax=ax, column = 'mlogp', by = 'assay')
plt.subplots_adjust(bottom=0.4, top=0.85)
ax.set_title('Assay Type vs -log10(ks p-value)')
plt.xticks(rotation=90)
ax.set_xlabel('Assay Type')
ax.set_ylabel('-log10(ks p-value)')
plt.savefig('figures/assaytype_ksp.png')

#group_df = df.groupby('manual_annotation_assay').mean()
#print(group_df['n_all'])
