Description:

Plot telescope, instruments, and photon raytraces.

Installation:

On unix systems with pip, ensure you configure in your root python (typing "which python" should give /usr/local/bin/python.)

./configure

Usage:

python phosim_visualizer.py ISCDirectory EventFile1Path ...

Usage with PhoSim integration:

./phosim examples/star -c examples/nobackground --visualize

Canopy Configuration Option:

This solution is robust but cumbersome:
Step 1.  Download "Enthought Canopy", and install it (available for Windows, MacOSX, and Linux)
Step 2.  Open Canopy, go to "Tools"----> "Package Manager".  Search for "astropy" and click "free install"
Step 3.  Follow the same procedure, search for "mayavi" and click "free install" 
Step 4.  Open Canopy, go to "Tools"----> "Canopy Terminal" and follow usage.

Notes:
PhoSim checks the transmission of the next surface before the raytrace. Many rays will appear to stop at the surface before the filter and detector.