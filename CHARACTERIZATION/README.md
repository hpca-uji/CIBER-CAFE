<p align="center">
  <img src="https://github.com/hpca-uji/CIBER-CAFE/blob/main/LOGOS/CIBER-CAFE_logo.jpeg" width="100" height="100">
</p>

## CHARACTERIZATION OF HOMOMORPHIC ENCRYPTION LIBRARIES

As a starting point of the CIBER-CAFE project, we have conducted an evaluation of the performance of the **matrix multiplication** in the context of Fully Homomorphic Encryption (FHE).

### EVALUATED LIBRARIES

- SEAL (https://github.com/microsoft/SEAL)
- PALISADE (http://palisade-crypto.org)
- HElib (https://homenc.github.io/HElib/)
- TFHE (https://tfhe.github.io/tfhe/)

### TESTS CONDUCTED

We offer an exhaustive evaluation that covers all the possible combinations of:
- Four different matrix dimensions: 4x4, 8x8, 16x16, 32x32
- Two schemes: BFV for matrices with integer values, CKKS for matrices with real values (except for TFHE)
- Four different libraries (SEAL, PALISADE, HElib, TFHE)
- HE key size of 128 and 256

### PLATFORMS

The tests have been conducted in such a way that we characterize the performance on eithers systems that can be equivalent to server nodes, or hardware that can form client/edge nodes.

This is the description of the server-like nodes tested:
- c00: equipped with 2 Intel Xeon 4210R processors, 128GB RAM memory.   
- c01: equipped with 2 Intel Xeon E5-2610 (8 cores each) processors, 32GB RAM memory.   
- c02: equipped with 2 AMD EPYC 7282 (16 cores each) processors, 256GB RAM memory.   
- c03: equipped with 2 Intel Xeon 6418H (24 cores each) processors, 2TB RAM memory.
  
This is the description of the client/edge-like nodes tested:
- NVIDIA Jetson AGX Xavier: NVIDIA Carmel Platform (8 cores) at 2.26 GHz, 32GB of DDR4 RAM. 
- NVIDIA Jetson AGX Orin: ARM A78AE (12 cores) at 2.20 GHz, 64GB of DDR5 RAM. 
- NVIDIA Jetson Nano: NVIDIA Cortex A57 (4 cores) at 1.48 GHz, 4GB of LPDDR4 RAM. 

### Contact information

For more information regarding this characterization, please contact:

- Manuel F. Dolz Zaragozá (dolzm@uji.es)

- Sandra Catalán Pallarés (catalans@uji.es)

- Rocío Carratalá-Sáez (rcarrata@uji.es)

### Acknowledgements

<p align="center">
  <img src="https://github.com/hpca-uji/CIBER-CAFE/blob/main/LOGOS/Banner_logos_funding.jpg" width="500">
</p>

<p align="center">
  <img src="https://github.com/hpca-uji/CIBER-CAFE/blob/main/LOGOS/UJI_logo.png" width="100">
</p>
