# SMILE
This folder contains codes and data for reproducing result from paper Accelerated Simultaneous Multislice Imaging via Linear Phase Modulated Extended Field of View (SMILE).
Source codes for several algorithms are inlcuded:
- SG, SPSG: https://mchiew.github.io/Tools.html
- ROCK-SPIRiT: https://github.com/obdemirel/ROCK_SPIRiT
- GRO-CAVA: https://github.com/OSU-CMR/GRO-CAVA
- SENSE: https://www.mathworks.com/help/matlab/ref/lsqr.html
- HICU: https://github.com/OSU-CMR/HICU



The brain dataset is from fastMRI: https://fastmri.med.nyu.edu/  
The perfusion dataset can be downloaded from 
- 3 slices single band perfusion: https://shorturl.at/tMas2
- 5 slices single band perfusion: https://shorturl.at/TOJmL
- SMILE MB = 3 perfusion: https://shorturl.at/16eYD
- SMILE MB = 5 perfusion: https://shorturl.at/OR1td, https://shorturl.at/EJTAU 

We provide the detailed recosntruction code for brain dataset, the preamble code can be adapted for 2D + t perfusion as well, where the tuned parameters can be found in the reference.

# Reference
Zhao, Shen, et al. "Whole Heart Perfusion with High-Multiband Simultaneous Multislice Imaging via Linear Phase Modulated Extended Field of View (SMILE)." arXiv preprint arXiv:2409.04353 (2024). (https://arxiv.org/abs/2409.04353)
