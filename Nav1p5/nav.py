""" analysis of DMS of SCN5A Voltage Sensor data """

import pandas as pd
import matplotlib.pyplot as plt

pos_codons = pd.read_csv('NM_198056.2_codons.csv', header = None)
pos_codons_list = pos_codons[0].tolist()
df = pd.read_csv('deep_mut_scan.csv')
df['possible'] = df['mutation'].isin(pos_codons_list)
print('number of variants: ' + str(len(df)))
posdf = df[df['possible'] == True]
print('number of possible variants: ' + str(len(posdf)))
imposdf = df[df['possible'] == False]
print('number of impossible variants: ' + str(len(imposdf)))

from scipy.stats import ks_2samp
from scipy.stats import ttest_ind

a = posdf['dms'].tolist()
b = imposdf['dms'].tolist()

print(ks_2samp(a,b))
print(ttest_ind(a,b,nan_policy='omit'))

fig,(ax1,ax2,ax3) = plt.subplots(3,1)
ax1.hist(df['dms'], bins = 40)
ax1.set_title('a and b')
ax2.hist(posdf['dms'], bins = 40)
ax2.set_title('a only')
ax3.hist(imposdf['dms'], bins = 40)
ax3.set_title('b only')
plt.savefig('DMS_histograms.png')
plt.close(fig)

posdf_sample = posdf.sample(n=50, random_state=42)
imposdf_sample = imposdf.sample(n=50, random_state=42)

c = posdf_sample['dms'].tolist()
d = imposdf_sample['dms'].tolist()

print('With random sampling of N=50 each')
print(ks_2samp(c,d))
print(ttest_ind(c,d,nan_policy='omit'))

fig,(ax1,ax2,ax3) = plt.subplots(3,1)
fig.suptitle('Sampled N=50')
ax1.hist(df['dms'], bins = 40)
ax1.set_title('All DMS scores')
ax2.hist(posdf_sample['dms'], bins = 40)
ax2.set_title('DMS scores in possible variants')
ax3.hist(imposdf_sample['dms'], bins = 40)
ax3.set_title('DMS scores in impossible variants')
plt.savefig('sampled_DMS_histograms.png')
plt.close(fig)
