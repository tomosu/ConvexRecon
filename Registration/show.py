from scipy import misc
import numpy as np
from PIL import Image
from glob import glob

names = glob('reg_images/cnn_result_20171017/*.npy')
for name in names:

    img = np.load(name)[0,:,:,0]
    img_name_0 =name.split(".")[0]
    img_id =img_name_0.split("/")[-1]
    misc.imsave('./ims_result/cnn{}.png'.format(img_id), img)


names = glob('reg_images/spect_rec/*.npy')
for i,name in enumerate(names):
    if i ==500:
        break
    img = np.load(name)
    img_name_0 =name.split(".")[0]
    img_id =img_name_0.split("/")[-1]
    misc.imsave('./ims_ref/spect_rec{}.png'.format(img_id), img)
