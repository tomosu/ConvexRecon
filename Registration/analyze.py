from scipy import misc
import numpy as np
from PIL import Image
from glob import glob

#originals2 =glob("/Users/tomohiro_suzuki/Work/LowDoseML/noise0.6_dataset/iter100/*.npy")
targets1 =glob("/Users/tomohiro_suzuki/Work/LowDoseML/mlem100_fbp_120000iter_noise0.6/*.npy")
targets2 =glob("/Users/tomohiro_suzuki/Work/LowDoseML/input_iter100_to_fbponly_120000iter_noise0.6/*.npy")
answer =glob("/Users/tomohiro_suzuki/Work/LowDoseML/input_iter100_to_fbponly_120000iter_noise0.6/*.npy")

error_sum_0 =0.0
error_sum_1 =0.0
error_sum_2 =0.0
denomi =0.0
for i in range(len(originals1)):

    #original1 =(np.load(originals1[i]))
    #original2 =np.load(originals2[i])
    target1 =np.load(targets1[i])[0,:,:,0]
    target2 =np.load(targets2[i])[0,:,:,0]

    print("o1mean:", np.mean(original1))
    print("o1mean:", np.mean(original2))
    print("t1mean:", np.mean(target1))
    print("t2mean:", np.mean(target2))

    diff0 =original1 -original2
    error0 =np.linalg.norm(diff0)
    error_sum_0 =error_sum_0 +error0

    diff1 =original1 -target1
    error1 =np.linalg.norm(diff1)
    error_sum_1 = error_sum_1 +error1

    diff2 =original1 -target2
    error2 =np.linalg.norm(diff2)
    error_sum_2 = error_sum_2 +error2

    denomi = denomi +1.0

print("mlem_recon:", error_sum_0/denomi)
print("mlem_and_fbp:", error_sum_1/denomi)
print("fbp:", error_sum_2/denomi)
