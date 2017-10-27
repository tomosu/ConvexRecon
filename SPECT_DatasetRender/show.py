from scipy import misc
import numpy as np
from PIL import Image
from glob import glob

names = glob('spect_rec_e0.01_a0.0003/0000*.npy')
for name in names:
    img = np.load(name)
    img_name_0 =name.split(".")[-2]
    img_id =img_name_0.split("/")[-1]
    misc.imsave('./ims/{}.png'.format(img_id), img)
