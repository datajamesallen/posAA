#!/usr/bin/python3

import sys
import requests
import time
import os
import pandas as pd 

def add_info(download = False):
    """
    add in extra information to the output.csv file from the mavedb API

    also adding longer abstracts, etc into the scoresets folders
    """
    anno_df = pd.read_csv('annotations.csv')
    removed_df = anno_df[anno_df['remove'] == 'REMOVE']
    anno_df = anno_df[anno_df['remove'] != 'REMOVE']
    removed_list = removed_df['urn'].tolist()

    if download:
        print('Downloading mavedb scoreset information')
        res = requests.get('https://www.mavedb.org/api/scoresets/?format=json')

        target_list = [['urn','keywords','title','short_desc','organism','pubmed_ids']]
        if res:
            scoreset_info = res.json()
            print('iterating through scoreset results')
            for scoreset in scoreset_info:
                # we need limit our analysis to protein coding data
                target_type = scoreset['target']['type']
                if target_type != 'Protein coding':
                    continue
                urn = scoreset['urn']
                # we want to make sure the urn is not in our removed list
                if urn in removed_list:
                    continue

                # create a directory for each target
                target_name = scoreset['target']['name']
                if not os.path.exists('scoresets/' + target_name):
                    os.makedirs(target_name)

                textfile = os.path.join('scoresets', target_name, urn, urn + '_experiment.txt')

                pubmed_urls = []
                pubmed_ids = []
                title = scoreset['title']
                short_desc = scoreset['short_description']
                organism = scoreset['target']['reference_maps'][0]['genome']['organism_name']
                keywords = []
                with open(textfile, 'w') as f:
                    f.write(title + '\n')
                    f.write('keywords:\n')
                    for item in scoreset['keywords']:
                        if item:
                            keyword = item['text']
                            f.write(keyword+'\n')
                            keywords.append(keyword)
                    f.write('Description\n')
                    f.write(short_desc + '\n')
                    f.write('Abstract\n')
                    f.write(scoreset['abstract_text'] + '\n')
                    f.write('Method\n')
                    f.write(scoreset['method_text'] + '\n')
                    f.write('Organism\n')
                    f.write(organism + '\n')
                    f.write('Pubmed\n')
                    for item in scoreset['pubmed_ids']:
                        if item:
                            pmurl = item['url']
                            f.write(pmurl+'\n')
                            pubmed_urls.append(pmurl)
                            pubmed_ids.append(item['identifier'])
                pubmed_string = '/'.join(pubmed_urls)
                keywords_string = '/'.join(keywords)
                target_row = [urn, keywords, title, short_desc, organism, pubmed_ids]
                target_list.append(target_row)

        info_df = pd.DataFrame(target_list[1:], columns = target_list[0])
        info_df.to_csv('auto_info.csv', index=False)

    else: # if not downloading
        info_df = pd.read_csv('auto_info.csv')

    return info_df

if __name__ == "__main__":
    add_info(download=True)

