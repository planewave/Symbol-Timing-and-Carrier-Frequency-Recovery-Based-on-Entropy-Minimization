This repository stores the source codes for my paper: 

    Liu, X., & Bousquet, J.-F. (2018). Symbol Timing and Carrier Frequency Recovery Based on Entropy Minimization. IEEE Access. 

It can be found in https://ieeexplore.ieee.org/document/8468150/

Codes are tested under Matlab 2018a (with communication toolbox), but should be OK for early version. 

I sorry that the repository is not well orgnized.
You can find the core algorithms, MRE and BMRE, from `mre.m` and `bmre.m`.
Other figures can be recreated using the following scripts.

* Fig. 1 to Fig. 4: `paper_plot.m`. **Make sure to run it section by section**
* Fig. 5: `carrier_entropy.m`
* Fig. 6: `MCRB_timing_nda5.m`

Hope it helps, and welcome to cite my paper. 
