import os
import dicom
import math
import numpy as np
from scipy import misc


def find_modality_images_path(path, modality_type):
    cp_list =os.listdir(path)
    if '.DS_Store' in cp_list:
        cp_list.remove('.DS_Store')
    if './' in cp_list:
        cp_list.remove('./')
    if '../' in cp_list:
        cp_list.remove('../')

    if '000000.dcm' in cp_list:
        d =dicom.read_file(path +'/000000.dcm')
        if d[0x0008, 0x0060].value == modality_type:
            return(path +'/')
        else:
            return None
    else:
        for cp in cp_list:
            if find_modality_images_path(path +'/' +cp, modality_type) ==None:
                continue
            else:
                return find_modality_images_path(path +'/' +cp, modality_type)




def show_pet_info(path):
    d =dicom.read_file(path)
    if d[0x0008, 0x0060].value =='PT':
        print(d[0x0018, 0x0050])
        print(d[0x0008, 0x0060])
        print(d[0x0028, 0x0010])
        print(d[0x0028, 0x0011])
        print(d[0x0028, 0x0100])
        print(d[0x0028, 0x0030])
        print(d[0x0054, 0x1310])
        print(d[0x0020, 0x0013])
        print(d[0x0028, 0x1052])
        print(d[0x0028, 0x1053])
        print(d[0x0020, 0x0032])
        print(d[0x0020, 0x1041])
        print(d[0x0008, 0x1090])
    else:
        print(d[0x0008, 0x0060])
    return d



def show_ct_info(path):
    d =dicom.read_file(path)
    if d[0x0008, 0x0060].value =='CT':
        print(d[0x0008, 0x0060])
        print(d[0x0018, 0x0050])
        print(d[0x0018, 0x0090])
        print(d[0x0018, 0x1100])
        print(d[0x0018, 0x1110])
        print(d[0x0018, 0x1111])
        print(d[0x0018, 0x1120])
        print(d[0x0018, 0x1140])
        print(d[0x0018, 0x1210])
        print(d[0x0020, 0x0013])
        print(d[0x0020, 0x0032])
        print(d[0x0020, 0x0037])
        print(d[0x0020, 0x1041])
        print(d[0x0028, 0x0010])
        print(d[0x0028, 0x0011])
        print(d[0x0028, 0x0030])
        print(d[0x0028, 0x0100])
        print(d[0x0028, 0x0101])
        print(d[0x0028, 0x0102])
        print(d[0x0028, 0x1050])
        print(d[0x0028, 0x1051])
        print(d[0x0028, 0x1052])
        print(d[0x0028, 0x1053])
    else:
        print(d[0x0008, 0x0060])
    return d



def get_dicom_images_sorted_by_slice(path):
    img_list =os.listdir(path)
    dl =[]
    for iname in img_list:
        dl.append(dicom.read_file(path +'/' +iname))

    dl.sort(key=lambda x:x[0x0020, 0x0013].value)
    return dl



def change_slope_intercept(image, slope, intercept):
    im =image *slope +intercept
    print("max:",im.max(), "min:",im.min())
    return im



def rescale_and_align_images(pet_path, ct_path):
    #resize and align
    pdl =get_dicom_images_sorted_by_slice(pet_path)
    cdl =get_dicom_images_sorted_by_slice(ct_path)

    pet_images =[]
    ct_images =[]
    for sl in range(min(len(pdl), len(cdl))):

        print('process:',sl)
        #rescale and align size
        pd =pdl[sl]
        cd =cdl[sl]

        scale_row =cd[0x0028, 0x0030].value[0] / pd[0x0028, 0x0030].value[0]
        scale_col =cd[0x0028, 0x0030].value[1] / pd[0x0028, 0x0030].value[1]

        pd_size =( pd[0x0028, 0x0010].value, pd[0x0028, 0x0010].value )
        cd_size =( cd[0x0028, 0x0010].value, cd[0x0028, 0x0010].value )
        resized_cd_size =( math.floor(cd_size[0]*scale_row), math.floor(cd_size[1]*scale_col) )

        cd_pix =cd.pixel_array
        resized_cd_pix =np.ones(shape =resized_cd_size)

        for r in range(cd_size[0]):
            for c in range(cd_size[1]):
                newr = min(resized_cd_size[0]-1, round(r *scale_row))
                newc = min(resized_cd_size[1]-1, round(c *scale_col))
                resized_cd_pix[newr, newc] =cd_pix[r, c]

        final_cd_pix = np.full(shape =pd_size, fill_value =-1000.0)

        ofst_r = int((pd_size[0] -resized_cd_size[0])/2)
        ofst_c = int((pd_size[1] -resized_cd_size[1])/2)

        c_slope =cd[0x0028, 0x1053].value
        c_intercept =cd[0x0028, 0x1052].value
        final_cd_pix[ofst_c : int(ofst_c+resized_cd_size[0]), ofst_r : int(ofst_r+resized_cd_size[1])] =change_slope_intercept(resized_cd_pix, c_slope, c_intercept)
        ct_images.append(final_cd_pix)

        p_slope =pd[0x0028, 0x1053].value
        p_intercept =pd[0x0028, 0x1052].value
        pet_images.append(change_slope_intercept(pd.pixel_array, p_slope, p_intercept))


    return pet_images, ct_images




#search path
studies_path ="../DICOM/PETCT/TCIA_soft_tissue_Sarcoma/DOI"
study_list =os.listdir(studies_path)
if '.DS_Store' in study_list:
    study_list.remove('.DS_Store')
if './' in study_list:
    study_list.remove('./')
if '../' in study_list:
    study_list.remove('../')

for i,sts in enumerate(study_list):
    root_path =studies_path +"/" +sts +"/"

    ct_path =find_modality_images_path(root_path, modality_type ='CT')
    pet_path =find_modality_images_path(root_path, modality_type ='PT')
    print(ct_path)
    show_ct_info(ct_path+'000000.dcm' )
    print(pet_path)
    show_pet_info(pet_path+'000000.dcm' )
    print('')

    pet_images, ct_images =rescale_and_align_images(pet_path, ct_path)
    for j in range(len(pet_images)):
        sdir ="reg_images/"
        np.save(sdir +'pet/{0:04d}{1:06d}.npy'.format(i,j), pet_images[j])
        np.save(sdir +'ct/{0:04d}{1:06d}.npy'.format(i,j), ct_images[j])
