import pandas as pd
import matplotlib.pyplot as plt
import math

df = pd.read_csv('output.csv')

fig, ax = plt.subplots(1,1)
ks_p = df['ks_p']
logp = [-1*math.log(x,10) for x in ks_p]
ax.scatter(df['n_all'],logp)
ax.set_xlim(0,20000)
ax.set_xlabel('total positions')
ax.set_ylabel('-log10(ks pvalue)')
plt.savefig('total_vs_ksp.png')