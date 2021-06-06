# combine all csv files
# create a folder called 'proc_data' inside
def bind_csv (args, regex, label=''):
    save_dir = args.f+'/'+args.c
    all_files = os.listdir(save_dir)
    all_csv = []
    for one_file in all_files:
        if regex.match (one_file):
            print ('adding '+one_file)
            all_csv.append (pd.read_csv (save_dir+'/'+one_file) )
    all_csv = pd.concat (all_csv, axis=0)
    csv_path = args.f+'/proc_data/'
    if not os.path.exists (csv_path): os.makedirs(csv_path)
    all_csv.to_csv (csv_path+args.c+'_'+label+'.csv')

if __name__ == '__main__':
    import argparse
    import re
    import os
    import pandas as pd

    parser = argparse.ArgumentParser()
    parser.add_argument('f', type=str, help='DAPI image folder')
    parser.add_argument('c', type=str, help='condition name')
    args = parser.parse_args ()

    regex = re.compile ('^.*nuc_seg.csv$')
    bind_csv (args, regex)
