# helper functions to fix up directories
# will likely get removed when the project is finished
# only necessary for organization

import os

for target in os.listdir('scoresets'):
    target_dir = os.path.join('scoresets', target)
    if target == ".DS_Store":
        continue
    for item in os.listdir(target_dir):
        item_dir = os.path.join(target_dir, item)
        if os.path.isdir(item_dir):
            for subfile in os.listdir(item_dir):
                subfile_dir = os.path.join(item_dir, subfile)
                last = subfile[24:]
                if last in ("possible", "impossible"):
                    os.rename(subfile_dir, subfile_dir+'.csv')
            continue
        if item == " experiment_info.txt":
            os.remove(item_dir)
            continue
        if item == ".DS_Store":
            continue
        urn = item[:23] # get the first 23 items
        urn_dir = os.path.join(target_dir, urn)
        os.rename(item_dir, os.path.join(urn_dir, item))