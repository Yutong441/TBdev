''' 
Read lif files using python-bioformats. Do not use read-lif or readlif. The
former is inaccurate. The latter does not work with lifext files. Bioformats
however have very poor documentation.

Reference:
https://ilovesymposia.com/2014/08/10/read-microscopy-images-to-numpy-arrays-with-python-bioformats/
'''
import numpy as np
import matplotlib.pyplot as plt
import tifffile as tif
import os
import bioformats as bf
import javabridge
javabridge.start_vm(class_path=bf.JARS)

def obtain_z_stack (bio_reader, series_num, time_index=0, verbose=True):
    z_stack, zxyc_img  = 0, []
    while True:
        if verbose: print ('processing image {}'.format (z_stack))
        try: 
            zxyc_img.append (bio_reader.read(z=z_stack, t=0, series=series_num,
                rescale=False))
            z_stack += 1
        except: break
    return np.transpose (np.stack (zxyc_img, axis=0), (3, 0, 1, 2) )

def display_z_stack (czxy_img, at_index=0, at_axis=0, along_axis=1, ncol=3,
        color='gray'):
    '''
    Args:
        `czxy_img`: numpy array in (channel, z, x, y) order
        `at_index`: select which index
        `at_axis`: select which axis
        `along_axis`: plot 2D plots along which axis
        `ncol`: number of columns of images
    '''
    s = [slice(None) for i in range(czxy_img.ndim)]
    s[at_axis] = slice(at_index, at_index+1)
    if along_axis > at_axis: along_axis -=1
    zxy_img = np.squeeze (czxy_img [s])
    show_img = np.moveaxis (zxy_img, along_axis, 0)

    print (show_img.shape)
    print (zxy_img.shape)
    num_chan = show_img.shape[0]
    nrow = np.ceil (num_chan/ncol )
    for index, one_slice in enumerate (show_img):
        plt.subplot (ncol, nrow, index+1)
        plt.imshow (one_slice, cmap=color)
        plt.axis('off')
    plt.show ()

def save_tif (czxy_img, root, save_3D=True, save_name='practice',
        channel_name=None):
    chan_num =czxy_img.shape[0] 
    if channel_name is None: 
        channel_name = ['chan'+str(i) for i in range(chan_num)]
    else:
        if len (channel_name) != chan_num:
            raise Exception ('Ensure the number of channel names match the \
                    number of channels in the supplied image')

    for index, chan in enumerate (channel_name):
        print ('saving '+chan)
        if save_3D:
            filename = os.path.join (root, save_name+'_'+chan+'.tiff') 
            tif.imsave (filename, czxy_img[index],compress=6)
        else:
            for z_index, xy_img in enumerate (zxy_img[index]):
                root_chan = os.path.join (root, chan)
                if not os.path.exists (root_chan): os.makedirs (root_chan)
                filename = os.path.join (root_chan,
                        save_name+'_'+str(z_index)+'.tiff') 
                tif.imsave (filename, xy_img,compress=6)
