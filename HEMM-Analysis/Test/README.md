
<p align="center">
  <img src="https://github.com/hpca-uji/CIBER-CAFE/blob/main/LOGOS/CIBER-CAFE_logo.jpeg" width="100" height="100">
</p>

# CIBER-CAFE: HEMM Analysis

## List of tests prepared

Libraries:

- SEAL or EVA: HE matrix-matrix product

- Plain (non encrypted) matrix-matrix product

Matrix dimensions:

- 1x1 

- 2x2

- 4x4

- 8x8

- 16x16

- 32x32

- 64x64

- 128x128

Key size: 128, 256

Precision: 40

Number of repetitions of each test: 10

## How to run the tests

```
cd Test
./run_tests.sh
```

## How to generate the plots in SC24-Poster

```
./parse_results.sh
```

Note that the plots will be stored in `CIBER-CAFE/HEMM-Analysis/Test/` with the names `Time128.png` and `Time256.png`

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

