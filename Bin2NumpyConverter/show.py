from scipy import misc
import numpy as np
from PIL import Image

img = np.load('./npy_output/PANCREAS_0001_000109.npy')
misc.imsave('test.png', img)
