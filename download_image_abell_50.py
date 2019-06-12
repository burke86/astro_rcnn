# generate 50 512x512 image at one of the abell cluster
# example: ACO 1656 : ra=194.9529 dec=27.9806

from download_image import *

if __name__ == "__main__":

    ra_cen=197.860331
    dec_cen=-1.336351
    shift_ra=512*0.27/3600
    shift_dec=512*0.27/3600

    output_dir = "testset"
    make_dir(output_dir)

    for n_ra in [-2,-1,0,1,2]:
        for n_dec in [-4,-3,-2,-1,0,1,2,3,4,5]:
            ra = ra_cen+n_ra*shift_ra
            dec = dec_cen+n_dec*shift_dec 
            size = 512
            number = check_set_number(output_dir)
            get_multi_images(ra,dec,size,"%s/set_%s" %(output_dir,number))

