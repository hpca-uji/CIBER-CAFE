
# HEMM-Optimization - Benchmarking Matrix Multiplication with Homomorphic Encryption

## 📋 Objective

**HEMM-Optimization** is a comprehensive benchmarking framework designed to evaluate and compare various matrix multiplication algorithms implemented over homomorphic encryption (HE) schemes. The project supports multiple HE libraries (OpenFHE, HElib) and implements advanced algorithms, including:
- Traditional matrix multiplication
- Rows × Columns decomposition
- Diagonals algorithm (with Halevi-Shoup packing inspiration)
- Jiang-Kim-Lauter-Song (HE-MatMult) algorithm ([paper](https://yongsoosong.github.io/files/papers/matrix.pdf))
- Rizomiliotis-Triakosias algorithm ([paper](https://dl.acm.org/doi/pdf/10.1145/3560810.3564267))
- Strassen technique variants
  - Strassen1x1x1: Strassen's method applied to 1×1×1 submatrices (single elements)
  - StrassenRizomiliotisTriakosia: Strassen's method applied to submatrices processed with the Rizomiliotis-Triakosias algorithm

The goal is to provide performance metrics, error analysis, and resource consumption data for encrypted matrix operations across different cryptographic parameters.

## 📁 Project Structure

```
hemm-optimization/
├── Tests/                         # Test implementations
│   ├── OpenFHE_Tests/             # OpenFHE benchmark executables
│   └── HElib_Tests/               # HElib benchmark executables
├── Analisis/                      # Analysis scripts and tools
├── Resultados/                    # Results directory (auto-generated)
├── ejecTests.sh                   # Main benchmark execution script
├── mat1_<size>.csv                # Input matrix 1 (auto-generated)
├── mat2_<size>.csv                # Input matrix 2 (auto-generated)
└── README.md                      # This file
```

## 📦 Installation

### Prerequisites

- **Linux/Unix** environment (for automated scripts and compatibility)
- **CMake** ≥ 3.12
- **GCC/G++** ≥ 9.0 or **Clang++** ≥ 10.0
- **Python 3.8+** (for result analysis)
- **wget** (for downloading dependencies)

### Manual Installation

#### Dependencies
- **GMP**, **NTL**, and **HEXL** libraries
- **OpenFHE**:
  - **OpenFHE-development**: [openfheorg/openfhe-development](https://github.com/openfheorg/openfhe-development)
  - **OpenFHE-HEXL** (If you want to use HEXL optimizations): [openfheorg/openfhe-hexl](https://github.com/openfheorg/openfhe-hexl)
- **HElib**: [homenc/HElib](https://github.com/homenc/HElib)
- **OpenMP** support for parallel execution (optional but recommended)

Test executables are located in `Tests/OpenFHE_Tests/` and `Tests/HElib_Tests/`. You can compile them using CMake:
```bash
# Compile HElib tests
cmake -S ./Tests/HElib_Tests -B .build/HElib_Tests/ -DHELIB_ROOT=$libInstall -DENABLE_OPENMP=$openmpInAlg
cmake --build .build/HElib_Tests/ -j$(nproc)

# Compile OpenFHE tests
cmake -S ./Tests/OpenFHE_Tests -B .build/OpenFHE_Tests/ -DOpenFHE_ROOT=$libInstall -DENABLE_OPENMP=$openmpInAlg
cmake --build .build/OpenFHE_Tests/ -j$(nproc)
```

## 📖 Usage

### Running Benchmarks

#### Automated Full Suite (Recommended)

```bash
bash ejecTests.sh
```

Results are saved in `Resultados/HEMM-Optimization_<timestamp>/` directory.

#### Automated Script Features

The `ejecTests.sh` script provides:

##### System Detection
- Automatic detection of CPU cores, RAM, and compiler versions
- SLURM cluster support for distributed testing
- Fallback to hostname-based configuration

##### Logging
- Timestamped, color-coded log output
- Real-time progress tracking
- Comprehensive error reporting

##### Configuration Options

Edit `ejecTests.sh` to customize:

###### Library Configuration

```bash
gmpVers="6.2.1"                                     # GMP version
gmpRuta=~/gmp-$gmpVers                              # GMP installation path

gf2x=False                                          # Set to True to install GF2X
gf2xVers="1.3.0"                                    # GF2X version
gf2xRuta=~/gf2x-$gf2xVers                           # GF2X installation path

ntlVers="11.5.1"                                    # NTL version
ntlRuta=~/ntl-$ntlVers                              # NTL installation path

hexlBranch="main"                                   # HEXL branch
hexlRuta=~/.git/HEXL                                # HEXL installation path

openFHEBranch=""                                    # OpenFHE branch
openFHERuta=~/.git/OpenFHEConfig                    # OpenFHE installation path
openFHEArgs=("-DNATIVE_SIZE=64" "-DNATIVE_SIZE=128") # OpenFHE build arguments

helibBranch=""                                      # HElib branch
helibRuta=~/.git/HElib                              # HElib installation path
helibArgs=("")                                      # HElib build arguments

check=True                                          # Set to True to run 'make check'
remove=True                                         # Set to True to remove source files after installation
```

###### Matrix Generation Parameters

```bash
precision=16                                        # Number of fractional digits
numNegativos=True/False                             # Allow negative numbers in matrices
ficMat1=mat1                                        # Name of matrix file 1
ficMat2=mat2                                        # Name of matrix file 2
```

###### Benchmark Parameters

```bash
copiarfichResult=False                              # Set to True to copy a previous results file if it exists
timeoutSec=3600                                     # Timeout per test in seconds (3600s = 1 hour)

# Test details
openmpInAlg=OFF                                     # Enable OpenMP in the algorithms
repeticiones=10                                     # Number of repetitions per test
reiteraciones=(1)                                   # Number of multiplication reiterations
tamsMatriz=(2 3 4 6 8 12 16 24 32 48 64 96 128 192 256 384 512 768 1024 1536 2048 3072 4096)  # Matrix sizes to test

# OpenFHE test parameters
OpenFHEAlg=("HEMatMult" "RizomiliotisTriakosia" "RowsXCols" "HaleviShoup" "Traditional" "Strassen1x1x1" "StrassenRizomiliotisTriakosia")
OpenFHERingDim=(8192 16384 32768 65536 131072 262144 524288 1048576 2097152 4194304 8388608 16777216 33554432)
OpenFHEMultDepth=(1 2 3 4 5 6 7 8 9 10)

# HElib test parameters
HElibAlg=("HEMatMult" "RizomiliotisTriakosia" "RowsXCols" "HaleviShoup" "Traditional" "Strassen1x1x1" "StrassenRizomiliotisTriakosia")
HElibParams=("16384 119 2" "32768 358 6" "32768 299 3" "32768 239 2" "65536 725 8" "65536 717 6" "65536 669 4" "65536 613 3" "65536 558 2" "131072 1445 8" "131072 1435 6" "131072 1387 5" "131072 1329 4" "131072 1255 3" "131072 1098 2" "262144 2940 8" "262144 2870 6" "262144 2763 5" "262144 2646 4" "262144 2511 3" "262144 2234 2")
HElibPrecisions=(40)
```

#### Individual Test Execution
Test executables require two input matrix files in CSV format. The automated script generates these files automatically (e.g., `mat1_4096.csv` and `mat2_4096.csv` for 4096×4096 matrices) in the project root directory if they don't exist. You can customize the matrix size, precision, and other parameters in the "Matrix Generation Parameters" section.

Example commands:
```bash
# OpenFHE test
./Resultados/HEMM-Optimization_<timestamp>/OpenFHE_Tests/<config>/matMultTest \
    --mat1 mat1_4096.csv \
    --mat2 mat2_4096.csv \
    --res results.csv \
    --alg Traditional \
    --m 256 --k 256 --n 256 \
    --ringDim 65536 \
    --multDepth 3 \
    --reIter 10

# HElib test
./Resultados/HEMM-Optimization_<timestamp>/HElib_Tests/<config>/matMultTest \
    --mat1 mat1_4096.csv \
    --mat2 mat2_4096.csv \
    --res results.csv \
    --alg Traditional \
    --m 256 --k 256 --n 256 \
    --param_m 16384 \
    --bits 119 \
    --precision 40 \
    --reIter 10
```

#### Help

```bash
./matMultTest --help
```

##### Expected Output

```bash
Usage: matMultTest [OPTIONS]
Options:
  --mat1 PATH           Path to matrix 1 CSV file
  --mat2 PATH           Path to matrix 2 CSV file
  --res PATH            Path to output results CSV file
  --alg TEXT            Algorithm to test (e.g., Traditional, RowsXCols, etc.)
  --m INTEGER           Number of rows in matrix 1
  --k INTEGER           Number of columns in matrix 1 / rows in matrix 2
  --n INTEGER           Number of columns in matrix 2
  # OpenFHE-specific parameters
  --ringDim INTEGER     Ring dimension (OpenFHE)
  --multDepth INTEGER   Multiplicative depth (OpenFHE)
  # HElib-specific parameters
  --param_m INTEGER     Parameter m (HElib)
  --bits INTEGER        Bits of security (HElib)
  --precision INTEGER   Precision for fixed-point encoding (HElib)
  # Common parameters
  --reIter INTEGER      Number of multiplication reiterations
```

## 📊 Output Structure

```
Resultados/
├── HEMM-Optimization_<date>_<time>/
│   ├── OpenFHE_Tests/
│   │   ├── <config1>/
│   │   │   └── matMultTest (executable)
│   │   └── <config2>/
│   ├── HElib_Tests/
│   ├── results.csv
│   └── ejecucion.log
```

## 📊 Output Format

Results are saved in CSV format with columns including:

- **Encryption Parameters**: `ringDim`, `param_m`, `multDepth`, etc.
- **Matrix Dimensions**: `m`, `k`, `n`
- **Performance**: `dur_MatrixMultiplication`, `dur_All`, etc.
- **Error Metrics**: `errorAbs_mean`, `errorRel_mean`, `errorSmape_mean`
- **Memory Usage**: `ram_CryptoContext`, `ram_MatrixMultiplication`, etc.

## 🐛 Troubleshooting

| Issue | Solution |
|-------|----------|
| dependency installation fails | Check write permissions in `~/.local` or change installation directory |
| CMake not found | Install: `apt install cmake` |
| Out of memory | Reduce `tamsMatriz` sizes or `repeticiones` count |
| Timeout errors | Increase `timeoutSec` in script |


## 📝 Example Workflow

```bash
# 1. Clone and setup
cd /home/usuario/.git/hemm-optimization

# 2. Configure for quick test
sed -i 's/repeticiones=10/repeticiones=2/' ejecTests.sh
sed -i 's/tamsMatriz=(.*)/tamsMatriz=(2 4 8)/' ejecTests.sh

# 3. Run tests
bash ejecTests.sh

# 4. Analyze results
cat Resultados/HEMM-Optimization_*/ejecucion.log
grep "RESULT\|ERROR" Resultados/HEMM-Optimization_*/ejecucion.log
```
