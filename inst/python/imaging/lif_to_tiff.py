#root='/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/results/imaging/'
#filename=root+'110221_ADAMredo_HLAG_CGB_AP2G_QUANT.lif'

import bioformats as bf
import argparse
import read_lif as rl
import os

parser = argparse.ArgumentParser()
parser.add_argument('f', type=str, help='image filename')
parser.add_argument('n', type=str, help='folder name')
parser.add_argument('--c', nargs='+', help='list of channels in the image')
args = parser.parse_args()

if __name__ == '__main__':
    rdr = bf.ImageReader(args.f, perform_init=True)
    for i in range (100):
        try:
            zxy_img = rl.obtain_z_stack (rdr, series_num=i)
            #display_z_stack (zxy_img, at_index=2, at_axis=0, along_axis=1, color=None)
            new_folder = os.path.dirname (args.f)+'/data/'+args.n
            if not os.path.exists (new_folder): os.makedirs (new_folder)
            rl.save_tif (zxy_img, new_folder, save_3D=True,
                    channel_name=args.c, save_name='series'+str(i))
        except: 
            break

    print ('finishing extracting tiff files')
    os._exit (0)
