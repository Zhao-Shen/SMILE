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
- 3 slices single band perfusion: https://drive.google.com/file/d/1dnYLOkC_2J-eDtJdnSwv6OMmU-4eD13v/view?usp=drive_link
- 5 slices single band perfusion: https://drive.google.com/file/d/1Jvx2vEPt7Y4V_X4BZft-JSkFW4mT42XE/view?usp=drive_link
- SMILE MB = 3 https://drive.google.com/file/d/1TNcwleiISbD-L9aNC701ApO2qoxaI1xu/view?usp=drive_link
- SMILE MB = 5 perfusion: https://drive.google.com/file/d/12X-vSvNC3vRAL89Q2CN82Qx9zWUQt-KA/view?usp=drive_link, https://drive.google.com/file/d/1_o_ocOLyuenRXEDHbDeYqonJLvLh6HCJ/view?usp=drive_link

We provide the detailed recosntruction code for brain dataset, the preamble code can be adapted for 2D + t perfusion as well, where the tuned parameters can be found in the reference.

# Reference
- Zhao, Shen, et al. "Whole heart perfusion with high‐multiband Simultaneous Multislice Imaging via Linear phase modulated Extended field of view (SMILE)." Magnetic Resonance in Medicine (2025). (https://onlinelibrary.wiley.com/doi/10.1002/mrm.30541?af=R)
- Zhao, Shen, et al. "Whole Heart Perfusion with High-Multiband Simultaneous Multislice Imaging via Linear Phase Modulated Extended Field of View (SMILE)." arXiv preprint arXiv:2409.04353 (2024). (https://arxiv.org/abs/2409.04353)
