import numpy as np
import pandas as pd
import skimage.io
import matplotlib.pyplot as plt
from organoid import Organ, stat_object
import trace_parent as tp

if __name__ == '__main__':
    import argparse
    import re
    import os
    parser = argparse.ArgumentParser()
    parser.add_argument('f', type=str, help='DAPI image folder')
    parser.add_argument('--n', type=str, nargs='+', help='name for the fluoroscence channels')
    args = parser.parse_args()

    all_files = os.listdir(args.f)
    regex = re.compile ('series[0-9]+_cyto_seg.tiff')
    regex2 = re.compile ('cyto_seg')
    regex3 = re.compile ('cyto_seg.tiff')

    log_f = open (args.f+'/'+'series_log.txt', 'w')
    log_f.write ('starting a new session \n')
    log_f.close ()
    log_f = open (args.f+'/'+'series_log.txt', 'a')
    for one_file in all_files:
        if regex.match (one_file):
            log_f.write ('segmenting '+one_file+'\n')

            seg_org = skimage.io.imread (args.f+'/'+one_file)
            img_one_file = regex2.sub ('nuc_seg', one_file)
            seg_img = skimage.io.imread (args.f+'/'+img_one_file)

            log_f.write ('image dimension= {}, {} \n'.format(*seg_img.shape))
            org = Organ (seg_org[:,:,np.newaxis], log_f)
            org.get_centers ()
            org.get_index ()

            img = Organ (seg_img[:,:,np.newaxis], log_f)
            img.get_centers ()
            img.get_index ()

            stats1 = stat_object (img.lab_mask)
            stats2 = stat_object (org.lab_mask)

            tree = tp.find_parent_fast (img, org)
            tree_df = pd.DataFrame (tree, columns=['child', 'parent'])

            # remove objects without parents
            tree_df = tree_df.loc [tree_df.iloc[:,1].values!=0]
            # need to minus 1 from index because python starts from 0 whereas the numbering
            # of 'child' and 'parent' fro 1.
            tree_df['cell_volume'] = np.array([img.pcounts [np.where (img.index_id==i)[0]] 
                    for i in tree_df.child.values.astype(int) ])
            tree_df['structure_volume']=np.array([org.pcounts[np.where (org.index_id==i)[0]] 
                    for i in tree_df.parent.values.astype(int)])

            for chan in args.n:
                log_f.write ('working on channel '+chan+'\n')
                if chan not in ['DAPI']:
                    img_fl_filename= regex2.sub (chan, one_file)
                    img_fl = skimage.io.imread (args.f+'/'+img_fl_filename)
                    if img_fl.shape[2] < img_fl.shape[1]:
                        img_fl = np.transpose (img_fl, axes=[2, 0, 1])
                    img_fl = img_fl.max (0)[:,:,np.newaxis]
                    # use the nuclear mask for TFAP2C:
                    if chan in ['TFAP2C']:
                        tree_df = tp.get_fluoro_fast (img, img_fl, tree_df, chan)
                    # use the entire cell border for the following channels:
                    if chan in ['HLAG', 'CGB']:
                        tree_df = tp.get_fluoro_fast (org, img_fl, tree_df, chan)

            tree_df['condition'] = args.f
            tree_df['series'] = one_file
            csv_filename = regex3.sub ('.csv', one_file)
            tree_df.to_csv (args.f+'/'+csv_filename)

            stats1['condition'] = args.f
            stats1['series'] = one_file
            stats1.to_csv (args.f+'/'+regex3.sub ('nuc_seg.csv', one_file))
            stats2['condition'] = args.f
            stats2['series'] = one_file
            stats2.to_csv (args.f+'/'+regex3.sub ('cyto_seg.csv', one_file))
    log_f.close ()
