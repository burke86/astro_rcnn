import numpy as np
import subprocess

# spiral galaxy

def make_spiral():
    # disk
    box = 2 # arceconds
    ra = np.random.normal(-box/3600,box/3600) # deg
    dec = np.random.normal(-box/3600,box/3600) # deg
    mag = np.random.normal(21.0,1)
    ax1 = np.clip(np.random.normal(3,2),0,8) # size of axis 1 in arcseconds
    ax2 = ax1 # size of axis 2 in arcseconds
    ax3 = ax1/8 # size of axis 3 in arcseconds
    phi = np.random.uniform(0,360) # polar angle in degrees
    theta = np.random.uniform(0,360) # position angle in degrees
    n = np.clip(np.random.normal(1,0.1),.1,4) # Sersic index
    fclump = np.clip(np.random.normal(.2,.08),0,.4) # fraction of light in clumps
    nclump = np.clip(np.random.normal(50,20),0,1000) # number of clumps
    rclump = np.clip(np.random.normal(0.3,0.1),0,1) # gaussian clump size in arcseconds
    fspiral = np.clip(np.random.normal(0.5,0.2),0,0.8) # fraction of light in spiral
    waspiral = np.clip(np.random.normal(20,10),0,90) # winding angle of spiral in degrees
    rbar = np.clip(np.random.normal(0.5,0.3),0,1) # spiral bar size in arcseconds
    rspiral = np.clip(np.random.normal(0.5,0.2),0,1) # spiral gaussian width in arcseconds
    nodoc = np.clip(np.random.normal(0.3,0.2),0,1) # ?
    thetaspiral = np.random.uniform(0,360) # position angle of spiral
    pars = (ra,dec,mag,ax1,ax2,ax3,phi,theta,n,fclump,nclump,rclump,fspiral,waspiral,rbar,rspiral,nodoc,thetaspiral)
    disk = ('object 0 %f %f %f ../sky/sed_flat.txt  0 0 0 0 0 0 sersicComplex %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f none none \n' % pars)
    # bulge
    mag = mag + 1
    ax1 = np.random.normal(1,.1) # size of axis 1 in arcseconds
    ax2 = ax1 # size of axis 2 in arcseconds
    ax3 = ax1 # size of axis 3 in arcseconds
    phi = 0 # polar angle in degrees
    theta = 0 # position angle in degrees
    n = np.clip(np.random.normal(4.0,0.1),1,10) # Sersic index
    fclump = 0 # fraction of light in clumps
    nclump = 0 # number of clumps
    rclump = 0 # gaussian clump size in arcseconds
    fspiral = 0 # fraction of light in spiral
    waspiral = 0 # winding angle of spiral in degrees
    rbar = 0 # spiral bar size in arcseconds
    rspiral = 0 # spiral gaussian width in arcseconds
    thetaspiral = 0 # position angle of spiral
    pars = (ra,dec,mag,ax1,ax2,ax3,phi,theta,n,fclump,nclump,rclump,fspiral,waspiral,rbar,rspiral,thetaspiral)
    bulge = ('object 0 %f %f %f ../sky/sed_flat.txt  0 0 0 0 0 0 sersicComplex %f %f %f %f %f %f %f %f %f %f %f %f %f %f none none \n' % pars)
    return header+disk+bulge
    # elliptical
for i in range(0,10):
    header = ('rightascension 0\ndeclination 0\nazimuth 0\naltitude 89\nrottelpos 0\nfilter 2\nseeing 0.67\nvistime 90.0\nnsnap 1\n')
    name = './examples/spiral%d' % i
    with open(name,'w+') as f:
        f.write(header)
        f.write(make_spiral())
        # make second galaxy 20 % of the time
        if np.random.uniform() > 0.2:
            f.write(make_spiral())
            # make second galaxy 10 % of the time
            if np.random.uniform() > 0.5:
                f.write(make_spiral())
    bashCommand = './phosim %s -c examples/training -i decam --ds9' % name
    subprocess.call(bashCommand.split())
