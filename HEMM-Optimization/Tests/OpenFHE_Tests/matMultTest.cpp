#include <openfhe.h>
#include <omp.h>
#include <map>

#include "../utils.cpp"
#include "utils_OpenFHE.cpp"
#include "matMultTradicional.cpp"
#include "matMultRowsXCols.cpp"
#include "matMultHaleviShoup.cpp"
#include "matMultHE.cpp"
#include "matMultRizomiliotisTriakosia.cpp"
#include "matMultStrassen.cpp"

#ifndef ENABLE_DEBUG
#define ENABLE_DEBUG 0
#endif

using namespace lbcrypto;

static const std::map<SecurityLevel, std::string> SecurityLevelToString = {
    {HEStd_128_classic, "128_classic"},
    {HEStd_192_classic, "192_classic"},
    {HEStd_256_classic, "256_classic"},
    {HEStd_128_quantum, "128_quantum"},
    {HEStd_192_quantum, "192_quantum"},
    {HEStd_256_quantum, "256_quantum"},
    {HEStd_NotSet, "NotSet"},
};

static const std::map<ScalingTechnique, std::string> ScalingTechniqueToString = {
    {FIXEDMANUAL, "FIXEDMANUAL"},
    {FIXEDAUTO, "FIXEDAUTO"},
    {FLEXIBLEAUTO, "FLEXIBLEAUTO"},
    {FLEXIBLEAUTOEXT, "FLEXIBLEAUTOEXT"},
    {NORESCALE, "NORESCALE"},
    {INVALID_RS_TECHNIQUE, "INVALID_RS_TECHNIQUE"}
};

static const std::map<ProxyReEncryptionMode, std::string> ProxyReEncryptionModeToString = {
    {NOT_SET, "NOT_SET"},
    {INDCPA, "INDCPA"},
    {FIXED_NOISE_HRA, "FIXED_NOISE_HRA"},
    {NOISE_FLOODING_HRA, "NOISE_FLOODING_HRA"},
};

static const std::map<MultipartyMode, std::string> MultipartyModeToString = {
    {INVALID_MULTIPARTY_MODE, "INVALID_MULTIPARTY_MODE"},
    {FIXED_NOISE_MULTIPARTY, "FIXED_NOISE_MULTIPARTY"},
    {NOISE_FLOODING_MULTIPARTY, "NOISE_FLOODING_MULTIPARTY"},
};

static const std::map<ExecutionMode, std::string> ExecutionModeToString = {
    {EXEC_EVALUATION, "EXEC_EVALUATION"},
    {EXEC_NOISE_ESTIMATION, "EXEC_NOISE_ESTIMATION"},
};

static const std::map<DecryptionNoiseMode, std::string> DecryptionNoiseModeToString = {
    {FIXED_NOISE_DECRYPT, "FIXED_NOISE_DECRYPT"},
    {NOISE_FLOODING_DECRYPT, "NOISE_FLOODING_DECRYPT"},
};

static const std::map<KeySwitchTechnique, std::string> KeySwitchTechniqueToString = {
    {INVALID_KS_TECH, "INVALID_KS_TECH"},
    {BV, "BV"},
    {HYBRID, "HYBRID"},
};

static const std::map<EncryptionTechnique, std::string> EncryptionTechniqueToString = {
    {STANDARD, "STANDARD"},
    {EXTENDED, "EXTENDED"},
};

static const std::map<MultiplicationTechnique, std::string> MultiplicationTechniqueToString = {
    {BEHZ, "BEHZ"},
    {HPS, "HPS"},
    {HPSPOVERQ, "HPSPOVERQ"},
    {HPSPOVERQLEVELED, "HPSPOVERQLEVELED"}
};

static const std::map<PlaintextEncodings, std::string> PlaintextEncodingsToString = {
    {INVALID_ENCODING, "INVALID_ENCODING"},
    {COEF_PACKED_ENCODING, "COEF_PACKED_ENCODING"},
    {PACKED_ENCODING, "PACKED_ENCODING"},
    {STRING_ENCODING, "STRING_ENCODING"},
    {CKKS_PACKED_ENCODING, "CKKS_PACKED_ENCODING"},
};

int TestMatMult(const std::vector<std::vector<double>>& A,
                 const std::vector<std::vector<double>>& B,
                 Algorithm alg,
                 size_t ringDim,
                 size_t multDepth,
                 usint dcrtBits,
                 usint firstMod,
                 size_t reiteraciones,
                 std::map<std::string, std::string>& statsResults) {
    measureBlock("Test", statsResults, [&]() {});
    int ret = 0;
    std::vector<std::vector<std::vector<double>>> result;
    measureBlock("All", statsResults, [&]() {
        CCParams<CryptoContextCKKSRNS> parameters;
        measureBlock("Params", statsResults, [&]() {
            parameters.SetRingDim(ringDim);
            parameters.SetMultiplicativeDepth(multDepth);
            #if NATIVEINT == 128 && !defined(__EMSCRIPTEN__)
                // Currently, only FIXEDMANUAL and FIXEDAUTO modes are supported for 128-bit CKKS bootstrapping.
                ScalingTechnique rescaleTech = FIXEDAUTO;
                if (dcrtBits > 78) dcrtBits = 78;
                if (firstMod > 89) firstMod = 89;
            #else
                // All modes are supported for 64-bit CKKS bootstrapping.
                ScalingTechnique rescaleTech = FLEXIBLEAUTO;
                if (dcrtBits > 59) dcrtBits = 59;
                if (firstMod > 60) firstMod = 60;
            #endif
            parameters.SetScalingModSize(dcrtBits);
            parameters.SetScalingTechnique(rescaleTech);
            parameters.SetFirstModSize(firstMod);
            if (!isRectangular(alg)) {
                size_t lado = 1 << (int)ceil(log2(std::max({A.size(), A[0].size(), B.size(), B[0].size()})));
                // Evitar pow(double) y usar multiplicación entera, luego castear a usint
                parameters.SetBatchSize(lado*lado);
            }
        });
        statsResults["ringDim"] = to_string_precise(parameters.GetRingDim());
        statsResults["batchSize"] = to_string_precise(parameters.GetBatchSize());
        statsResults["multDepth"] = to_string_precise(parameters.GetMultiplicativeDepth());
        statsResults["dcrtBits"] = to_string_precise(parameters.GetScalingModSize());
        statsResults["firstMod"] = to_string_precise(parameters.GetFirstModSize());
        statsResults["securityLevel"] = SecurityLevelToString.at(parameters.GetSecurityLevel());
        statsResults["scalingTechnique"] = ScalingTechniqueToString.at(parameters.GetScalingTechnique());
        statsResults["proxyReEncryptionMode"] = ProxyReEncryptionModeToString.at(parameters.GetPREMode());
        statsResults["multipartyMode"] = MultipartyModeToString.at(parameters.GetMultipartyMode());
        statsResults["executionMode"] = ExecutionModeToString.at(parameters.GetExecutionMode());
        statsResults["decryptionNoiseMode"] = DecryptionNoiseModeToString.at(parameters.GetDecryptionNoiseMode());
        statsResults["keySwitchTechnique"] = KeySwitchTechniqueToString.at(parameters.GetKeySwitchTechnique());
        statsResults["encryptionTechnique"] = EncryptionTechniqueToString.at(parameters.GetEncryptionTechnique());
        statsResults["multiplicationTechnique"] = MultiplicationTechniqueToString.at(parameters.GetMultiplicationTechnique());
        
        #if ENABLE_DEBUG
        std::cout << "Starting CryptoContext Generation..." << std::endl;
        #endif
        CryptoContext<DCRTPoly> cc;
        try {
            measureBlock("CryptoContext", statsResults, [&]() {
                // Usar la fábrica genérica de CryptoContext. Si su versión de OpenFHE requiere otra función,
                // modifique esta línea conforme al error del compilador.
                cc = GenCryptoContext(parameters);
                cc->Enable(PKE);
                cc->Enable(KEYSWITCH);
                cc->Enable(LEVELEDSHE);
                cc->Enable(ADVANCEDSHE);
            });
        } catch (const std::exception& e) {
            newException(e, statsResults, "Error during CryptoContext generation");
            // Indicar error y salir del lambda. La variable 'ret' está capturada por referencia desde el
            // scope exterior (TestMatMult), así que asignamos el código de error aquí y retornamos desde
            // el lambda (tipo void). Después del measureBlock comprobaremos 'ret' y devolveremos.
            ret = 3;
            return;
        }
        // Imprimir modulus elemento a elemento
        auto params = cc->GetCryptoParameters()->GetElementParams()->GetParams();
        statsResults["modulus"] = "";
        std::vector<std::string> modulus(params.size());
        for (size_t i = 0; i < params.size(); ++i) {
            modulus[i] = params[i]->GetModulus().ToString();
            statsResults["modulus"] += modulus[i];
            if (i < params.size() - 1) statsResults["modulus"] += ";";
        }
        #if ENABLE_DEBUG
        std::cout << "Modulus:";
        for (size_t i = 0; i < modulus.size(); ++i) std::cout << " " << modulus[i];
        std::cout << std::endl;
        #endif

        result = std::vector<std::vector<std::vector<double>>>(reiteraciones, std::vector<std::vector<double>>(A.size(), std::vector<double>(B[0].size(), 0.0)));
        switch (alg) {
            case Traditional:
                ret = MatMultTradicional(cc, A, B, result, reiteraciones, statsResults);
                break;
            case RowsXCols:
            case RowsXColsSum:
                ret = MatMultRowsXCols(cc, A, B, result, reiteraciones, statsResults, alg==RowsXColsSum);
                break;
            case HaleviShoup:
            case HaleviShoupFastRot:
                ret = MatMultHaleviShoup(cc, A, B, result, reiteraciones, statsResults, alg==HaleviShoupFastRot);
                break;
            case HEMatMult:
            case HEMatMultFastRot:
                ret = MatMultHE(cc, A, B, result, reiteraciones, statsResults, alg==HEMatMultFastRot || alg==RHEMatMultFastRot);
                break;
            case RizomiliotisTriakosia:
            case RizomiliotisTriakosiaFastRot:
            case RizomiliotisTriakosiaSums:
            case RizomiliotisTriakosiaSumsFastRot:
                ret = MatMultRizomiliotisTriakosia(cc, A, B, result, reiteraciones, statsResults,
                    alg==RizomiliotisTriakosiaSums || alg==RizomiliotisTriakosiaSumsFastRot,
                    alg==RizomiliotisTriakosiaFastRot || alg==RizomiliotisTriakosiaSumsFastRot);
                break;
            case Strassen1x1x1:
            //case StrassenHEMatMult:
            case StrassenRizomiliotisTriakosia:
            case StrassenRizomiliotisTriakosiaFastRot:
            case StrassenRizomiliotisTriakosiaSums:
            case StrassenRizomiliotisTriakosiaSumsFastRot:
                ret = MatMultStrassen(cc, A, B, result, reiteraciones, statsResults, 0, StrassenSubAlgorithm.at(alg));
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
    std::cout << ">\n\t--m\t<rows of A>\n\t--k\t<cols of A/rows of B>\n\t--n\t<cols of B>\n\t--ringDim\t<ring dimension, power of 2>\n\t--multDepth\t<multiplicative depth, at least 1>\n\t--dcrtBits\t<bit size of each CRT modulus, at most 78 for 128-bit and 59 for 64-bit>\n\t--firstMod\t<bit size of the first CRT modulus, at most 89 for 128-bit and 60 for 64-bit>\n\t--reIter\t<reiterate the multiplication, at least 1>\n\t--rep\t<repetition count, at least 0>\n\t-h|--help\tShow this help message" << std::endl;
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
        "ringDim", "batchSize", "multDepth", "dcrtBits", "firstMod",
        "securityLevel", "scalingTechnique", "encryptionTechnique",
        "keySwitchTechnique", "multiplicationTechnique", "proxyReEncryptionMode",
        "multipartyMode", "executionMode", "decryptionNoiseMode",
        // Modulus y otros detalles de contexto
        "modulus",
        // Estadísticas de contexto
        "ctxA_size", "ctxB_size", "ctxC_size", "ctx_noiseScaleDeg", "ctx_scalingFactor", "ctx_hopLevel",
        "ctx_level_min", "ctx_level_max", "ctx_level_mean", "ctx_level_median",
        "ctx_level_p10", "ctx_level_p25", "ctx_level_p75", "ctx_level_p90", "ctx_level_std",
        // Estadísticas de plaintext
        "ptx_length",
        "ptx_logError_min", "ptx_logError_max", "ptx_logError_mean", "ptx_logError_median",
        "ptx_logError_p10", "ptx_logError_p25", "ptx_logError_p75", "ptx_logError_p90", "ptx_logError_std",
        "ptx_logPrecision_min", "ptx_logPrecision_max", "ptx_logPrecision_mean", "ptx_logPrecision_median",
        "ptx_logPrecision_p10", "ptx_logPrecision_p25", "ptx_logPrecision_p75", "ptx_logPrecision_p90", "ptx_logPrecision_std",
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
         ringDim = 1<<13,
         multDepth = 1,
         dcrtBits = 78,
         firstMod = 89,
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
        } else if (arg == "--ringDim") {
            if (i + 1 >= argc) { std::cerr << "--ringDim requires a value\n"; printHelp(); return -1; }
            ringDim = std::stoul(argv[++i]);
            if (ringDim & (ringDim - 1)) {
                std::cerr << "ringDim must be a power of two." << std::endl;
                printHelp();
                return -1;
            }
        } else if (arg == "--multDepth") {
            if (i + 1 >= argc) { std::cerr << "--multDepth requires a value\n"; printHelp(); return -1; }
            multDepth = std::stoul(argv[++i]);
            if (multDepth < 1) {
                std::cerr << "multDepth must be at least 1." << std::endl;
                printHelp();
                return -1;
            }
        } else if (arg == "--dcrtBits") {
            if (i + 1 >= argc) { std::cerr << "--dcrtBits requires a value\n"; printHelp(); return -1; }
            dcrtBits = std::stoul(argv[++i]);
        } else if (arg == "--firstMod") {
            if (i + 1 >= argc) { std::cerr << "--firstMod requires a value\n"; printHelp(); return -1; }
            firstMod = std::stoul(argv[++i]);
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
    statsResults["ringDim"] = to_string_precise(ringDim);
    statsResults["multDepth"] = to_string_precise(multDepth);
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
    
    int ret = TestMatMult(A, B, alg, ringDim, multDepth, 78, 89, reiteraciones, statsResults);
    
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