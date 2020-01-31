""" Download the proper data required from mavedb for this projects analysis """

import requests
import time
import os

print('Downloading mavedb scoreset information')
res = requests.get('https://www.mavedb.org/api/scoresets/?format=json')

if res:
    scoreset_info = res.json()
    print('iterating through scoreset results')
    for scoreset in scoreset_info:
        # we need limit our analysis to protein coding data
        target_type = scoreset['target']['type']
        if target_type != 'Protein coding':
            continue
        urn = scoreset['urn']
        # create a directory for each target
        target_name = scoreset['target']['name']
        if not os.path.exists('scoresets/' + target_name):
            os.makedirs(target_name)

        # get the sequence and save a file for it in the target dir
        seq = scoreset['target']['reference_sequence']['sequence']
        seq_filename = 'scoresets/' + target_name + '/'+ urn + '_seq.fasta'
        with open(seq_filename, 'w') as f:
            f.write(seq)

        # this is the url where we can download the scoreset
        scoreset_url = 'https://www.mavedb.org/scoreset/' + urn + '/scores/'
        res = requests.get(scoreset_url)
        data_filename = 'scoresets/' + target_name + '/' + urn + '_data.csv'
        with open(data_filename,'wb') as f:
            f.write(res.content)
        print(urn + ' -- data written to file')

        # to limit the amount of strain we put on their servers
        # and prevent time out errors
        time.sleep(10)