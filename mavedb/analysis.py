#!/usr/bin/python3

"""
runs batch analysis on mavedb data 

Goal: compare the distributions of variants possible with single amino acid changes
to those that are not possible, for differences in functional parameters
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re

# testing ks test from r package
from rpy2.robjects.vectors import FloatVector
from rpy2.robjects.packages import importr
r_stats = importr('stats')

from scipy.stats import ks_2samp
#from scipy.stats import ttest_ind

output_list = [['target','urn','n_multi','n_syn','mean_score_all','std_score_all','n_all','mean_score_pos',
                'std_score_pos','n_pos','mean_score_impos','std_score_impos',
                'n_impos','delta_pos_impos','ks_stat','ks_p','r_ks_stat','r_ks_p']]

# iterate through targets and files
for target in os.listdir('scoresets'):
    target_dir = os.path.join('scoresets',target)
    if not os.path.isdir(target_dir):
        continue
    for urn in os.listdir(target_dir):
        urn_dir = os.path.join(target_dir, urn)
        if not os.path.isdir(urn_dir):
            continue
        for filename in os.listdir(urn_dir):
            if filename.endswith('_data.csv'):
                # the urn number that identifies the data set
                urn = filename[:-9]
                print('Analyzing data for: ' + urn)
                # find the corresponding posAA file
                posAA_filename = urn + '_posAA.csv'

                # Begin analysis
                pos_codons = pd.read_csv(os.path.join(urn_dir,posAA_filename), header = None)
                pos_codons_list = pos_codons[0].tolist()

                # open data file
                df = pd.read_csv(os.path.join(urn_dir,filename), skiprows = 4)
                df['single_aa'] = df.apply(lambda x: len(x['hgvs_pro'].split(';')) == 1, axis=1)
                df['syn'] = df.apply(lambda x: len(re.findall('=', x['hgvs_pro']))==1, axis = 1)
                n_multi = len(df[df['single_aa'] == False])
                n_syn = len(df[df['syn'] == True])
                df = df[df['single_aa'] == True]
                df = df[df['syn'] == False]

                # calculate possible variants
                df['possible'] = df['hgvs_pro'].isin(pos_codons_list)
                n_all = len(df)
                #print('number of variants: ' + str(pos_var))
                posdf = df[df['possible'] == True]
                n_pos = len(posdf)
                #print('number of possible variants: ' + str(len(posdf)))
                imposdf = df[df['possible'] == False]
                n_impos = len(imposdf)
                #print('number of impossible variants: ' + str(len(imposdf)))

                posdf.to_csv(os.path.join(urn_dir,urn + '_possible.csv'))
                imposdf.to_csv(os.path.join(urn_dir,urn + '_impossible.csv'))

                mean_score_all = df['score'].mean()
                mean_score_pos = posdf['score'].mean()
                mean_score_impos = imposdf['score'].mean()

                std_score_all = df['score'].std()
                std_score_pos = posdf['score'].std()
                std_score_impos = imposdf['score'].std()

                delta_pos_impos = mean_score_pos - mean_score_impos

                poslist = posdf['score'].tolist()
                imposlist = imposdf['score'].tolist()

                ks_stat, ks_p = ks_2samp(poslist,imposlist,mode='exact')

                r_poslist = FloatVector(poslist)
                r_imposlist = FloatVector(imposlist)

                r_result = r_stats.ks_test(r_poslist, r_imposlist, exact=True)
                #print(r_result)
                r_ks_stat = r_result.rx('statistic')[0][0]
                r_ks_p = r_result.rx('p.value')[0][0]

                #ks_stat_exact, ks_p_exact = ks_2samp(poslist,imposlist,mode='exact')
                #print(ttest_ind(a,b,nan_policy='omit'))

                # new histogram design
                fig,ax = plt.subplots(1,1)
                bins = 40
                a = posdf['score']
                b = imposdf['score']
                a_w = np.empty(a.shape)
                a_w.fill(1/a.shape[0])
                b_w = np.empty(b.shape)
                b_w.fill(1/b.shape[0])
                ax.hist([a,b], bins, color = ['blue','red'], weights = [a_w,b_w], label = ['possible','impossible'])
                #ax.hist(posdf['score'], bins, color ='blue', alpha=0.5, label='possible', density=True)
                #ax.hist(imposdf['score'], bins, color ='yellow', alpha=0.5, label='impossible', density=True)
                plt.xlabel('Score')
                plt.ylabel('Probability')
                plt.suptitle('Normalized Histogram of ' + target + ' (' + urn + ')')
                plt.title('K-S p value = ' + str(r_ks_p))
                plt.legend()
                plt.savefig(os.path.join(urn_dir, urn+'_histograms.png'))
                plt.savefig('figures/histograms/' + target + '_' + urn + '.png')
                plt.close(fig)
                """
                # old histogram code
                fig,(ax1,ax2,ax3) = plt.subplots(3,1, sharex=True)
                ax1.hist(df['score'], bins = 40)
                ax1.set_title('all data')
                ax2.hist(posdf['score'], bins = 40)
                ax2.set_title('possible variants only')
                ax3.hist(imposdf['score'], bins = 40)
                ax3.set_title('impossible variants only')
                plt.savefig(os.path.join(urn_dir,urn + '_histograms.png'))
                plt.close(fig)
                """

                output_row = [target, urn, n_multi, n_syn, mean_score_all, std_score_all, n_all, mean_score_pos,
                          std_score_pos, n_pos, mean_score_impos, std_score_impos,
                          n_impos, delta_pos_impos, ks_stat, ks_p, r_ks_stat, r_ks_p]

                output_list.append(output_row)

with open('output.csv', 'w') as f:
    for row in output_list:
        f.write(','.join(str(e) for e in row)+'\n')