#/usr/bin/python3

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

from add_info import add_info

# first we need to load our file that includes the annotations for each urn.
# this is because we have annotated some urns that we want to exclude from our analysis


anno_df = pd.read_csv('annotations.csv')
info_df = add_info(download=False)

all_info_df = pd.merge(info_df, anno_df, how = 'outer', on ='urn')
removed_df = all_info_df[all_info_df['remove'] == 'REMOVE']
removed_list = removed_df['urn'].tolist()
non_remove = all_info_df['remove'] != 'REMOVE'
all_info_df = all_info_df[non_remove]
all_info_df = all_info_df.drop(['remove','reason'],axis=1)

human_df = all_info_df[all_info_df['organism'] == 'Homo sapiens']
human_list = human_df['urn'].tolist()
print(set(human_df['gene'].tolist()))

output_list = [['target','urn','n_multi','n_syn','n_stop','mean_score_all','std_score_all','n_all','mean_score_pos',
                'std_score_pos','n_pos','mean_score_impos','std_score_impos',
                'n_impos','delta_pos_impos','ks_stat','ks_p','r_ks_stat','r_ks_p']]

from matplotlib.backends.backend_pdf import PdfPages

# create a pdf for our histograms to all go in
pp = PdfPages('figures/histograms.pdf')
pp_delta = PdfPages('figures/delta_scores.pdf')
pp_syn = PdfPages('figures/delta_syn.pdf')
pp_gmd = PdfPages('figures/isgnomAD.pdf')

# iterate through targets and files
for target in os.listdir('scoresets'):
    target_dir = os.path.join('scoresets',target)
    if not os.path.isdir(target_dir):
        continue
    for urn in os.listdir(target_dir):
        # we want to check to make sure this urn is not in our removal list
        if urn in removed_list:
            continue
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

                # this will give all of the info for the given urn
                urn_df = all_info_df[all_info_df['urn'] == urn]
                gene = urn_df.iloc[0]['gene']

                # Begin analysis
                pos_codons = pd.read_csv(os.path.join(urn_dir,posAA_filename), header = None)
                pos_codons_list = pos_codons[0].tolist()

                # open data file
                df = pd.read_csv(os.path.join(urn_dir,filename), skiprows = 4)

                # figure out what type of variant each of these are
                df['stop'] = df.apply(lambda x: len(re.findall('Ter', x['hgvs_pro']))==1, axis = 1)
                df['single_aa'] = df.apply(lambda x: len(x['hgvs_pro'].split(';')) == 1, axis=1)
                df['syn'] = df.apply(lambda x: len(re.findall('=', x['hgvs_pro']))==1, axis = 1)
                n_stop = len(df[df['stop'] == True])
                n_multi = len(df[df['single_aa'] == False])
                n_syn = len(df[df['syn'] == True])
                df = df[df['stop'] == False]
                df = df[df['single_aa'] == True]
                df = df[df['syn'] == False]

                # get the amino acid number
                df['aa_num'] = df.apply(lambda x: re.findall('[0-9]+', x['hgvs_pro'])[0] if len(re.findall('[0-9]+', x['hgvs_pro'])) >= 1 else None, axis = 1)

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

                # figure out the difference in mean possible - impossible per site (aaNum)
                grouped_posdf = posdf.groupby('aa_num').agg({'score':['mean']})
                grouped_posdf.columns = ['_'.join(col).strip() for col in grouped_posdf.columns.values]
                grouped_posdf = grouped_posdf.reset_index()

                grouped_imposdf = imposdf.groupby('aa_num').agg({'score':['mean']})
                grouped_imposdf.columns = ['_'.join(col).strip() for col in grouped_imposdf.columns.values]
                grouped_imposdf = grouped_imposdf.reset_index()

                mergedf = grouped_posdf.merge(grouped_imposdf, on = 'aa_num', how = 'outer', suffixes = ('_pos','_impos'))
                mergedf['delta_score'] = mergedf['score_mean_pos'] - mergedf['score_mean_impos']

                #print(mergedf)
                fig1,ax1 = plt.subplots(1,1)
                non_null_df = mergedf[mergedf['delta_score'].notnull()]
                ax1.hist(non_null_df['delta_score'])
                plt.suptitle('Difference in score per codon Possible - Impossible for ' + target)
                plt.savefig(os.path.join(urn_dir, urn+'_delta.png'))
                plt.savefig('figures/delta_scores/' + target + '_' + urn + '.png')
                pp_delta.savefig(fig1)
                plt.close(fig1)

                if isinstance(gene, str): # this checks that there is a human gene associated with this data
                    """ only gets urns that will perform analysis based on gnomAD data """
                    # i am not sure why but missing numbers are qualified as floats right now
                    # this should get rid of missing values in this category
                    gmd_file = '../gnomAD/gnomADv2_' + gene + '.csv'
                    gmd_df = pd.read_csv(gmd_file)
                    syn_df = gmd_df[gmd_df['Annotation'] == 'synonymous_variant']
                    syn_df['aa_num'] = syn_df.apply(lambda x: re.findall('[0-9]+', x['Consequence'])[0] if len(re.findall('[0-9]+', x['Consequence'])) >= 1 else None, axis = 1)
                    syn_df['is_syn'] = syn_df.apply(lambda x: 1 if x['Allele Count'] > 0 else 0, axis=1)
                    grouped_syn_df = syn_df.groupby('aa_num').agg({'is_syn':['sum']})
                    grouped_syn_df.columns = ['_'.join(col).strip() for col in grouped_syn_df.columns.values]
                    grouped_syn_df = grouped_syn_df.reset_index()
                    #print(grouped_syn_df)
                    grouped_syn_df = grouped_syn_df[['aa_num','is_syn_sum']]
                    merged_gmd_syn_df = grouped_syn_df.merge(mergedf, on = 'aa_num', how = 'outer')
                    merged_gmd_syn_df = merged_gmd_syn_df.fillna(0)

                    # delta_score boxplots
                    fig, ax = plt.subplots(1,1)
                    merged_gmd_syn_df.boxplot(column = 'delta_score', by = 'is_syn_sum', ax=ax)
                    plt.suptitle(gene)
                    plt.savefig('figures/delta_score_syn/' + gene + '_syn.png')
                    pp_syn.savefig(fig)
                    plt.close(fig)

                    # figure out the functional distribution of gnomAD data
                    mis_df = gmd_df[gmd_df['Annotation'] == 'missense_variant']
                    mis_df['aa_num'] = mis_df.apply(lambda x: re.findall('[0-9]+', x['Consequence'])[0] if len(re.findall('[0-9]+', x['Consequence'])) >= 1 else None, axis = 1)
                    mis_df['is_mis'] = mis_df.apply(lambda x: 1 if x['Allele Count'] > 0 else 0, axis=1)
                    mis_df = mis_df[['Protein Consequence','Allele Count','is_mis']]
                    merged_gmd_mis_df = mis_df.merge(df, right_on = 'hgvs_pro', left_on = 'Protein Consequence', how = 'right')
                    merged_gmd_mis_df['Allele Count'] = merged_gmd_mis_df['Allele Count'].fillna(0)
                    merged_gmd_mis_df['in_gnomAD'] = merged_gmd_mis_df.apply(lambda x: True if x['Allele Count'] > 0 else False, axis=1)

                    # isgnomAD boxplots
                    fig, ax = plt.subplots(1,1)
                    merged_gmd_mis_df.boxplot(column = 'score', by = 'in_gnomAD', ax=ax)
                    plt.suptitle(gene)
                    plt.savefig('figures/isgnomAD/' + gene + '_isgnomAD.png')
                    pp_gmd.savefig(fig)
                    plt.close(fig)

                    gmd_len = len(merged_gmd_mis_df[merged_gmd_mis_df['in_gnomAD'] == True])
                    print(gmd_len)

                    # isgnomAD histograms
                    if gmd_len > 0:
                        fig,ax = plt.subplots(1,1)
                        bins = 40
                        a = merged_gmd_mis_df[merged_gmd_mis_df['in_gnomAD'] == True]['score']
                        print(a)
                        b = merged_gmd_mis_df[merged_gmd_mis_df['in_gnomAD'] == False]['score']
                        print(b)
                        a_w = np.empty(a.shape)
                        a_w.fill(1/a.shape[0])
                        b_w = np.empty(b.shape)
                        b_w.fill(1/b.shape[0])
                        ax.hist(a, bins, color = 'blue', weights = a_w, label = 'in gnomAD', alpha = 0.5)
                        ax.hist(b, bins, color = 'red', weights = b_w, label = 'not in gnomAD', alpha = 0.5)
                        plt.xlabel('Score')
                        plt.ylabel('Probability')
                        plt.suptitle(gene)
                        plt.savefig('figures/isgnomAD/' + gene + '_isgnomAD_hist.png')
                        plt.close(fig)

                # write them out to file for later investigation/debugging
                posdf.to_csv(os.path.join(urn_dir,urn + '_possible.csv'))
                imposdf.to_csv(os.path.join(urn_dir,urn + '_impossible.csv'))

                # perform some basic statistics
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

                # create a histogram to show to difference in scores
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
                if r_ks_p < 0.001:
                    plt.title('K-S p value < 0.001')
                else:
                    plt.title('K-S p value = ' + str(r_ks_p))
                plt.legend()
                plt.savefig(os.path.join(urn_dir, urn+'_histograms.png'))
                plt.savefig('figures/histograms/' + target + '_' + urn + '.png')
                pp.savefig(fig)
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

                output_row = [target, urn, n_multi, n_syn, n_stop, mean_score_all, std_score_all, n_all, mean_score_pos,
                          std_score_pos, n_pos, mean_score_impos, std_score_impos,
                          n_impos, delta_pos_impos, ks_stat, ks_p, r_ks_stat, r_ks_p]

                output_list.append(output_row)

pp.close()
pp_delta.close()
pp_syn.close()
pp_gmd.close()

output_df = pd.DataFrame(output_list[1:], columns = output_list[0])
final_output = pd.merge(output_df, all_info_df, how='outer',on='urn')
final_output.to_csv('output.csv',index=False)

