import numpy as np
import matplotlib.pyplot as plt
import skimage.io

root='/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/results/imaging/data/all/'
seg_org = skimage.io.imread (root+'/'+'series41_cyto_seg.tiff')
cgb = skimage.io.imread (root+'/'+'series39_CGB.tiff')
hlag = skimage.io.imread (root+'/'+'series41_HLAG.tiff')

cgb = cgb.max (0)
hlag = hlag.max(0)

plt.imshow (cgb)
plt.imshow (seg_org, alpha=0.5)
plt.show()

lab_seg = skimage.color.label2rgb (seg_org, cgb, bg_label=0)
plt.imshow (lab_seg)
plt.show()

lab_seg = skimage.color.label2rgb (seg_org, hlag, bg_label=0)
plt.imshow (lab_seg)
plt.show()
plt.imshow (seg_org)

hlag = skimage.io.imread (root+'/'+'series0_HLAG.tiff')
hlag = hlag.max(2)
hlagc = hlag.clip (0, 50)
plt.imshow (hlagc)
plt.show()
