#include <helib/helib.h>
#include <omp.h>
#include <map>

#include "../utils.cpp"
#include "utils_HElib.cpp"
#include "matMultTradicional.cpp"
#include "matMultRowsXCols.cpp"
#include "matMultHaleviShoup.cpp"
#include "matMultHE.cpp"
#include "matMultRizomiliotisTriakosia.cpp"
#include "matMultStrassen.cpp"

#ifndef ENABLE_DEBUG
#define ENABLE_DEBUG 0
#endif

using namespace helib;

int TestMatMult(const std::vector<std::vector<double>>& A,
                 const std::vector<std::vector<double>>& B,
                 Algorithm alg,
                 size_t param_m,
                 size_t bits,
                 size_t c,
                 size_t precision,
                 size_t reiteraciones,
                 std::map<std::string, std::string>& statsResults) {
    measureBlock("Test", statsResults, [&]() {});
    int ret = 0;
    std::vector<std::vector<std::vector<double>>> result;
    if (!isRectangular(alg) && (param_m/4) < (1 << ((int)ceil(log2(std::max({A.size(), A[0].size(), B.size(), B[0].size()})))*2))){
        newException(statsResults, "Not enough slots for the given matrices");
        return 3;
    }
    measureBlock("All", statsResults, [&]() {
        ContextBuilder<CKKS> contextBuilder;
        measureBlock("Params", statsResults, [&]() {
            contextBuilder = ContextBuilder<CKKS>()
                .m(param_m)
                .bits(bits)
                .precision(precision)
                .c(c);
        });
        statsResults["param_m"] = to_string_precise(param_m);
        statsResults["bits"] = to_string_precise(bits);
        statsResults["c"] = to_string_precise(c);
        statsResults["precision"] = to_string_precise(precision);
        
        #if ENABLE_DEBUG
        std::cout << "Starting CryptoContext Generation..." << std::endl;
        #endif
        std::unique_ptr<Context> cc_ptr;
        try {
            measureBlock("CryptoContext", statsResults, [&]() {
                // Construir el Context dinámicamente usando buildPtr() (evita move/assign)
                cc_ptr.reset(contextBuilder.buildPtr());
            });
        } catch (const std::exception& e) {
            newException(e, statsResults, "Error during CryptoContext generation");
            // Indicar error y salir del lambda. La variable 'ret' está capturada por referencia desde el
            // scope exterior (TestMatMult), así que asignamos el código de error aquí y retornamos desde
            // el lambda (tipo void). Después del measureBlock comprobaremos 'ret' y devolveremos.
            ret = 3;
            return;
        }
        if (!cc_ptr) {
            newException(std::runtime_error("Failed to build Context"), statsResults, "Context construction returned null");
            ret = 3;
            return;
        }
        statsResults["securityLevel"] = cc_ptr->securityLevel();
        statsResults["nSlots"] = cc_ptr->getNSlots();
        std::cout << "Number of Slots: " << cc_ptr->getNSlots() << std::endl;
        result = std::vector<std::vector<std::vector<double>>>(reiteraciones, std::vector<std::vector<double>>(A.size(), std::vector<double>(B[0].size(), 0.0)));
        switch (alg) {
            case Traditional:
                ret = MatMultTradicional(*cc_ptr, A, B, result, reiteraciones, statsResults);
                break;
            case RowsXCols:
            case RowsXColsSum:
                ret = MatMultRowsXCols(*cc_ptr, A, B, result, reiteraciones, statsResults, alg==RowsXColsSum);
                break;
            case HaleviShoup:
            case HaleviShoupShift:
                ret = MatMultHaleviShoup(*cc_ptr, A, B, result, reiteraciones, statsResults, alg==HaleviShoupShift);
                break;
            case HEMatMult:
                ret = MatMultHE(*cc_ptr, A, B, result, reiteraciones, statsResults);
                break;
            case RizomiliotisTriakosia:
                ret = MatMultRizomiliotisTriakosia(*cc_ptr, A, B, result, reiteraciones, statsResults);
                break;
            case Strassen1x1x1:
            //case StrassenHEMatMult:
            case StrassenRizomiliotisTriakosia:
                ret = MatMultStrassen(*cc_ptr, A, B, result, reiteraciones, statsResults, 0, StrassenSubAlgorithm.at(alg));
                break;
            default:
                newException(statsResults, "Algorithm " + AlgorithmToString.at(alg) + " not implemented.");
                ret = -1;
                break;
        }
    });
    // Si ocurrió un error dentro del bloque (por ejemplo, durante la generación del contexto),
    // 'ret' habrá sido actualizado dentro del lambda; propagar el retorno inmediatamente.
    if (ret != 0) return ret;
    return calcError(A, B, result, reiteraciones, statsResults);
}

void printHelp (){
    std::cout << "Usage: matmult\n\t--mat1\t<path>\n\t--mat2\t<path>\n\t--res\t<path>\n\t--alg\t<algorithm:";
    for (const auto& alg : AlgorithmToString) {
        std::cout << " " << alg.second << "(" << static_cast<int>(alg.first) << ")";
    }
    std::cout << ">\n\t--m\t<rows of A>\n\t--k\t<cols of A/rows of B>\n\t--n\t<cols of B>\n\t--param_m\t<param_m, power of 2>\n\t--bits\t<bits size, at least 1>\n\t--c\t<c>\n\t--precision\t<precision>\n\t--reIter\t<reiterate the multiplication, at least 1>\n\t--rep\t<repetition count, at least 0>\n\t-h|--help\tShow this help message" << std::endl;
}

int main(int argc, char* argv[]) {
    std::string mat1_path, mat2_path, res_path;
    std::vector<std::string> resHeader = {
        // Parámetros de paralelismo
        "algParallelism",
        // Dimensiones de matrices
        "m", "k", "n",
        // Identificación y configuración general
        "algorithm", "reIteraciones", "repeticion",
        // Parámetros de cifrado
        "param_m", "bits", "c", "precision",
        "securityLevel", "nSlots",
        // Estadísticas de contexto
        // Estadísticas de Ciphertext
        "ctxA_size",
        "ctxA_capacity_min", "ctxA_capacity_max", "ctxA_capacity_mean", "ctxA_capacity_median",
        "ctxA_capacity_p10", "ctxA_capacity_p25", "ctxA_capacity_p75", "ctxA_capacity_p90", "ctxA_capacity_std",
        "ctxA_errorBound_min", "ctxA_errorBound_max", "ctxA_errorBound_mean", "ctxA_errorBound_median",
        "ctxA_errorBound_p10", "ctxA_errorBound_p25", "ctxA_errorBound_p75", "ctxA_errorBound_p90", "ctxA_errorBound_std",
        "ctxB_size",
        "ctxB_capacity_min", "ctxB_capacity_max", "ctxB_capacity_mean", "ctxB_capacity_median",
        "ctxB_capacity_p10", "ctxB_capacity_p25", "ctxB_capacity_p75", "ctxB_capacity_p90", "ctxB_capacity_std",
        "ctxB_errorBound_min", "ctxB_errorBound_max", "ctxB_errorBound_mean", "ctxB_errorBound_median",
        "ctxB_errorBound_p10", "ctxB_errorBound_p25", "ctxB_errorBound_p75", "ctxB_errorBound_p90", "ctxB_errorBound_std",
        "ctxC_size",
        "ctxC_capacity_min", "ctxC_capacity_max", "ctxC_capacity_mean", "ctxC_capacity_median",
        "ctxC_capacity_p10", "ctxC_capacity_p25", "ctxC_capacity_p75", "ctxC_capacity_p90", "ctxC_capacity_std",
        "ctxC_errorBound_min", "ctxC_errorBound_max", "ctxC_errorBound_mean", "ctxC_errorBound_median",
        "ctxC_errorBound_p10", "ctxC_errorBound_p25", "ctxC_errorBound_p75", "ctxC_errorBound_p90", "ctxC_errorBound_std",
        // Estadísticas de error absoluto
        "errorAbs_min", "errorAbs_max", "errorAbs_mean", "errorAbs_median",
        "errorAbs_p10", "errorAbs_p25", "errorAbs_p75", "errorAbs_p90", "errorAbs_std",
        // Estadísticas de error relativo
        "errorRel_min", "errorRel_max", "errorRel_mean", "errorRel_median",
        "errorRel_p10", "errorRel_p25", "errorRel_p75", "errorRel_p90", "errorRel_std",
        // Estadísticas de error SMAPE
        "errorSmape_min", "errorSmape_max", "errorSmape_mean", "errorSmape_median",
        "errorSmape_p10", "errorSmape_p25", "errorSmape_p75", "errorSmape_p90", "errorSmape_std",
        // Medidas de tiempo
        "dur_Test", "dur_Params", "dur_CryptoContext", "dur_KeyGeneration",
        "dur_Preprocess", "dur_Encode", "dur_Encrypt",
        "dur_MatrixMultiplication",
        "dur_Decrypt", "dur_Decode", "dur_Postprocess",
        "dur_All",
        // Medidas de RAM
        "ram_Test", "ram_Params", "ram_CryptoContext", "ram_KeyGeneration",
        "ram_Preprocess", "ram_Encode", "ram_Encrypt",
        "ram_MatrixMultiplication",
        "ram_Decrypt", "ram_Decode", "ram_Postprocess",
        "ram_All",
        // Excepciones
        "exceptions"
    };
    uint m = 2,
         k = 2,
         n = 2,
         param_m = 1<<14,
         bits = 119,
         c = 2,
         precision = 40,
         reiteraciones = 1,
         repeticion = 0;
    Algorithm alg = Traditional;
    // Parse command line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--mat1") {
            if (i + 1 >= argc) { std::cerr << "--mat1 requires a path\n"; printHelp(); return -1; }
            mat1_path = argv[++i];
        } else if (arg == "--mat2") {
            if (i + 1 >= argc) { std::cerr << "--mat2 requires a path\n"; printHelp(); return -1; }
            mat2_path = argv[++i];
        } else if (arg == "--res") {
            if (i + 1 >= argc) { std::cerr << "--res requires a path\n"; printHelp(); return -1; }
            res_path = argv[++i];
        } else if (arg == "--alg") {
            if (i + 1 >= argc) { std::cerr << "--alg requires an index\n"; printHelp(); return -1; }
            std::string val = argv[++i];
            try {
                alg = static_cast<Algorithm>(std::stoul(val));
            } catch (const std::invalid_argument&) {
                auto it = std::find_if(AlgorithmToString.begin(), AlgorithmToString.end(),
                    [&val](const std::pair<Algorithm, std::string>& pair) { return pair.second == val; });
                if (it != AlgorithmToString.end()) {
                    alg = it->first;
                } else {
                    std::cerr << "Unknown algorithm: " << val << std::endl;
                    printHelp();
                    return -1;
                }
            } catch (const std::out_of_range&) {
                std::cerr << "Algorithm index out of range: " << val << std::endl;
                printHelp();
                return -1;
            }
        } else if (arg == "--m") {
            if (i + 1 >= argc) { std::cerr << "--m requires a value\n"; printHelp(); return -1; }
            m = std::stoul(argv[++i]);
            if (m == 0) {
                std::cerr << "m must be at least 1." << std::endl;
                printHelp();
                return -1;
            }
        } else if (arg == "--k") {
            if (i + 1 >= argc) { std::cerr << "--k requires a value\n"; printHelp(); return -1; }
            k = std::stoul(argv[++i]);
            if (k == 0) {
                std::cerr << "k must be at least 1." << std::endl;
                printHelp();
                return -1;
            }
        } else if (arg == "--n") {
            if (i + 1 >= argc) { std::cerr << "--n requires a value\n"; printHelp(); return -1; }
            n = std::stoul(argv[++i]);
            if (n == 0) {
                std::cerr << "n must be at least 1." << std::endl;
                printHelp();
                return -1;
            }
        } else if (arg == "--param_m") {
            if (i + 1 >= argc) { std::cerr << "--param_m requires a value\n"; printHelp(); return -1; }
            param_m = std::stoul(argv[++i]);
            if (param_m & (param_m - 1)) {
                std::cerr << "param_m must be a power of two." << std::endl;
                printHelp();
                return -1;
            }
        } else if (arg == "--bits") {
            if (i + 1 >= argc) { std::cerr << "--bits requires a value\n"; printHelp(); return -1; }
            bits = std::stoul(argv[++i]);
            if (bits < 1) {
                std::cerr << "bits must be at least 1." << std::endl;
                printHelp();
                return -1;
            }
        } else if (arg == "--c") {
            if (i + 1 >= argc) { std::cerr << "--c requires a value\n"; printHelp(); return -1; }
            c = std::stoul(argv[++i]);
        } else if (arg == "--precision") {
            if (i + 1 >= argc) { std::cerr << "--precision requires a value\n"; printHelp(); return -1; }
            precision = std::stoul(argv[++i]);
        } else if (arg == "--reIter") {
            if (i + 1 >= argc) { std::cerr << "--reIter requires a value\n"; printHelp(); return -1; }
            reiteraciones = std::stoul(argv[++i]);
            if (reiteraciones < 1) {
                std::cerr << "reiteraciones must be at least 1." << std::endl;
                printHelp();
                return -1;
            }
        } else if (arg == "--rep") {
            if (i + 1 >= argc) { std::cerr << "--rep requires a value\n"; printHelp(); return -1; }
            repeticion = std::stoul(argv[++i]);
            if (repeticion < 0) {
                std::cerr << "repeticion must be at least 0." << std::endl;
                printHelp();
                return -1;
            }
        } else if (arg == "-h" || arg == "--help") {
            printHelp();
            return 0;
        } else {
            std::cerr << "Unknown argument: " << arg << std::endl;
            printHelp();
            return -1;
        }
    }

    if (!isRectangular(alg) && (m != k || m != n)) {
        std::cerr << "For square algorithms, m, n, and k must be equal." << std::endl;
        printHelp();
        return -1;
    }

    if (reiteraciones > 1 && isRectangular(alg) && k != n) {
        std::cerr << "For multiple iterations with rectangular algorithms, matrix 2 must be square (k must equal n)." << std::endl;
        printHelp();
        return -1;
    }


    if (mat1_path.empty() || mat2_path.empty()) {
        std::cerr << "Both --mat1 and --mat2 must be provided.\n";
        printHelp();
        return -1;
    }

    // Read the required submatrices
    #if ENABLE_DEBUG
    std::cout << "Reading matrices from files..." << std::endl;
    #endif
    std::vector<std::vector<double>> A = readMatrix(mat1_path, m, k),
                                     B = readMatrix(mat2_path, k, n);
    #if ENABLE_DEBUG
    std::cout << "Loaded matrices: A(" << m << "x" << k << "):" << std::endl;
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < k; ++j) {
            std::cout << A[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "Loaded matrices: B(" << k << "x" << n << "):" << std::endl;
    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < n; ++j) {
            std::cout << B[i][j] << " ";
        }
        std::cout << std::endl;
    }
    #endif
    
    std::map<std::string, std::string> statsResults;
    statsResults["param_m"] = to_string_precise(param_m);
    statsResults["bits"] = to_string_precise(bits);
    statsResults["c"] = to_string_precise(c);
    statsResults["precision"] = to_string_precise(precision);
    statsResults["algorithm"] = AlgorithmToString.at(alg);
    statsResults["reIteraciones"] = to_string_precise(reiteraciones);
    statsResults["m"] = to_string_precise(m);
    statsResults["n"] = to_string_precise(n);
    statsResults["k"] = to_string_precise(k);
    statsResults["repeticion"] = to_string_precise(repeticion);
    #if ENABLE_OPENMP
    #if USE_OMP_TASKLOOP
    statsResults["algParallelism"] = "MultiThreadTask(" + to_string_precise(omp_get_max_threads())+")";
    #else
    statsResults["algParallelism"] = "MultiThreadFor(" + to_string_precise(omp_get_max_threads())+")";
    #endif
    #else
    statsResults["algParallelism"] = "SingleThread";
    #endif
    
    int ret = TestMatMult(A, B, alg, param_m, bits, c, precision, reiteraciones, statsResults);
    
    #if ENABLE_DEBUG
    std::cout << "Results:" << std::endl;
    for (const auto &p : statsResults) {
        std::cout << "\t" << p.first << ":\t" << p.second << std::endl;
    }
    std::cout << "Writing results to file..." << std::endl;
    #endif

    writeStats(res_path, resHeader, statsResults);

    return ret;
}