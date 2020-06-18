import numpy as np

NUM_CUTOUTS = 50
CUTOUT_SEMI = 0.01193

def gen_cutouts_repl(num_cutouts, cutout_semi):
    tract_data = [9812, 147.92194 + 0.01194, 1.6576 - 0.01194, 1.40294 + 0.01194, 1.65613 - 0.01194]
    ra_start = tract_data[1]
    ra_interval = tract_data[2]
    dec_start = tract_data[3]
    dec_interval = tract_data[4]
    tract_number = 9812

    fid = open('cutout_list.txt', 'a')
    filters = ['HSC-R', 'HSC-G', 'HSC-Z']
    for i in range(0, NUM_CUTOUTS):
        ra_seed = np.random.rand()
        dec_seed = np.random.rand()
        ra = (ra_seed * ra_interval) + ra_start
        dec = (dec_seed * dec_interval) + dec_start
        for j in range(0, len(filters)):
            fid.write('\n')
            fid.write('pdr2_dud ')
            fid.write(str(tract_number))
            fid.write(' ')
            fid.write(str(ra))
            fid.write(' ')
            fid.write(str(dec))
            fid.write(' ')
            fid.write(str(cutout_semi))
            fid.write(' ')
            fid.write(str(cutout_semi))
            fid.write(' ')
            fid.write(filters[j])
            fid.write(' ')
            fid.write('true ')
            fid.write('true ')
            fid.write('false ')
            fid.write('coadd')




    fid.close()
        
        

        





if __name__ == "__main__":
    gen_cutouts_repl(NUM_CUTOUTS, CUTOUT_SEMI)
