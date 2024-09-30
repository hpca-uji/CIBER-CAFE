<p align="center">
  <img src="https://github.com/hpca-uji/CIBER-CAFE/blob/main/LOGOS/CIBER-CAFE_logo.jpeg" width="200" height="200">
</p>

# CIBER-CAFE: HEMM Analysis

We provide tools to evaluate the Homomorphic Encryption Matrix-Matrix Product (HEMM).

In particular, we use the SEAL library, and compare its performance to that of EVA.

## Software Dependencies

[Microsoft EVA version 1.0.1](https://github.com/microsoft/EVA/releases/tag/v1.0.1)

[Microsoft SEAL version 3.6](https://github.com/microsoft/SEAL/tree/3.6.4)

[SEAL-Python version 3.6.2](https://github.com/Huelse/SEAL-Python/tree/3.6.2)

Note that, for coherency, SEAL-Python should be linked to the previously installed Microsoft SEAL version 3.6, instead of that included when obtaining SEAL-Python. Thus, it is better to remove SEAL-Python/SEAL.

### Additional software dependencies

- Python3, pip3
  
- Python3 modules: numpy, psutil, columnar, pandas

- Pybind11

- Git

- CMake

## Compilation/Installation - Option 1: Docker

**STEP 1) Obtain the Docker image**

Download the docker image tar file from:

[https://drive.google.com/file/d/1jZ7zmd19ScIlcnN3ARXNRGw0zRfqceYg/view?usp=share_link](https://drive.google.com/file/d/1jZ7zmd19ScIlcnN3ARXNRGw0zRfqceYg/view?usp=share_link)

**STEP 2) Load the Docker image**

```
docker load --input docker-image.tar
```

**STEP 3) Run the Docker image**

```
docker run -it ciber-cafe:sc24-hemm bash
```

## Compilation/Installation - Option 2: From scratch

After clonning this repository, complete the following steps.

Everything is described assuming that, in the same folder, these will be available:

```
HEMM-Analysis
├── SEAL        - Microsoft SEAL version 3.6.4
├── SEAL-Python - SEAL-Python version 3.6.2
├── EVA         - EVA version 1.0.1
└── Test
    ├── SEAL_Test
    ├── EVA_Test
    └── Plain_computation_Test
```


**STEP 1) Make sure the dependencies are satisfied, otherwise install what is needed:**

```
apt update
apt install python3
apt install python3-pip
pip3 install numpy
pip3 install psutil
pip3 install columnar
pip3 install pandas
apt-get install pybind11
apt-get install pybind11-dev
apt-get install git
apt-get install cmake
```

**STEP 2) Install Microsoft SEAL version 3.6:**

```
git clone -b v3.6.4 https://github.com/microsoft/SEAL.git
cd SEAL
cmake -S . -B build -DSEAL_USE_MSGSL=OFF -DSEAL_USE_ZLIB=OFF -DSEAL_USE_ZSTD=OFF -DSEAL_THROW_ON_TRANSPARENT_CIPHERTEXT=OFF
cmake --build build
cmake --install build
cd ..
```

If `cmake --build build` fails as described in [this issue](https://github.com/microsoft/SEAL/issues/674), solve the problems reated to this issue by adding `#include <mutex>` in locks.h (located in SEAL/native/src/seal/util/locks.h)

**STEP 3) Install SEAL-Python**

```
git clone --branch 3.6.2 https://github.com/Huelse/SEAL-Python.git
cd SEAL-Python
rm -rf SEAL
sed -i 's/SEAL/..SEAL/g' setup.py 
python3 setup.py build_ext -i
cp seal.*.so ../Test/SEAL_Test
cd ..
```

**STEP 4) Install EVA**

```
apt install cmake libboost-all-dev libprotobuf-dev protobuf-compiler clang
update-alternatives --install /usr/bin/cc cc /usr/bin/clang 100
update-alternatives --install /usr/bin/c++ c++ /usr/bin/clang++ 100
git clone https://github.com/microsoft/EVA.git
cd EVA
git submodule update --init
cd third_party/pybind11
git pull origin master
cd ../../
cmake .
make -j
python3 -m pip install -e python/
cd ..
```

*Observation 1*

If issues with uint_32.t happen, add (before the other includes) in eva/ckks/ckks_config.cpp:

`#include <stdint.h>`

*Observation 2*

If issues with incomplete type ‘PyFrameObject’ happen as in this [issue](https://github.com/numpy/numpy/issues/21422), solve them by:

```
python3 -m pip install cython==3.0.0a10 
python3 -m pip install git+https://github.com/cython/cython@master
```

## Test the HEMM

Check the README in the Test folder

## Contact information

For more information regarding this characterization, please contact:

- Manuel F. Dolz Zaragozá (dolzm@uji.es)

- Sandra Catalán Pallarés (catalans@uji.es)

- Rocío Carratalá-Sáez (rcarrata@uji.es)


## Acknowledgements

<p align="center">
  <img src="https://github.com/hpca-uji/CIBER-CAFE/blob/main/LOGOS/Banner_logos_funding.jpg" width="500">
</p>

<p align="center">
  <img src="https://github.com/rociocarratalasaez/CIBER-CAFE/blob/main/LOGOS/UJI_logo.png" width="100">
</p>

