<p align="center">
  <img src="https://github.com/rociocarratalasaez/CIBER-CAFE/blob/main/LOGOS/UJI_logo.png" width="100">
</p>

<p align="center">
  <img src="https://github.com/rociocarratalasaez/CIBER-CAFE/blob/main/LOGOS/CIBER-CAFE_logo.png" width="100" height="100">
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

## How to obtain a CSV file with all the data

```
python3 parse.py
```

Note that the combined data will be stored in `Results/combined_results.csv`

### Contact information

For more information regarding this characterization, please contact:

- Manuel F. Dolz Zaragozá (dolzm@uji.es)

- Sandra Catalán Pallarés (catalans@uji.es)

- Rocío Carratalá-Sáez (rcarrata@uji.es)

### Acknowledgements

<p align="center">
  <img src="https://github.com/rociocarratalasaez/CIBER-CAFE/blob/main/LOGOS/Banner_logos_funding.jpg" width="500">
</p>
