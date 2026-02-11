#include <openfhe.h>
#include <omp.h>
#include <map>

using namespace lbcrypto;

std::set<int> rotRowXCols(const int colsA) {
    std::set<int> rot;
    for (int i = 0; (1 << i) < colsA; ++i) {
        rot.insert(1 << i);
        rot.insert(-(1 << i));
    }
    return rot;
}

template <class Element>
void matMultRowsXCols(const std::vector<Ciphertext<Element>>& encA,
                      const std::vector<Ciphertext<Element>>& encB,
                      std::vector<Ciphertext<Element>>& encC,
                      std::map<std::string, std::string>& statsResults,
                      size_t& reiteraciones,
                      size_t k = 0,
                      const size_t reitActual=0,
                      const bool sum = false) {
    bool error = false;
    bool isCKKS = (encA[0]->GetEncodingType() == CKKS_PACKED_ENCODING);
    CryptoContext<Element> cc = encA[0]->GetCryptoContext();
    if (k==0)
        k = encA[0]->GetSlots();
    std::vector<Ciphertext<Element>> encTemp = std::vector<Ciphertext<Element>>(encA.size()*encB.size());
    encC = std::vector<Ciphertext<Element>>(encA.size());
    std::vector<Plaintext> ptxtMasks = std::vector<Plaintext>(encB.size());
    #if USE_OMP_TASKLOOP
    #pragma omp taskloop shared(ptxtMasks)
    #else
    #pragma omp parallel for shared(ptxtMasks) schedule(dynamic)
    #endif
    for (size_t j = 0; j < encB.size(); j++)
        ptxtMasks[j] = makePlaintext<Element>(cc, isCKKS, encB.size(), [j](int ell){return ell == j;});
    #pragma omp parallel for collapse(2) shared(encTemp, reiteraciones, statsResults, error) schedule(dynamic)
    for (size_t i = 0; i < encA.size(); i++)
        for (size_t j = 0; j < encB.size(); j++)
            try {
                if (sum) {
                    encTemp[i*encB.size()+j] = cc->EvalMult(cc->EvalSum(cc->EvalMult(encA[i], encB[j]), encA[0]->GetSlots()), ptxtMasks[j]);
                } else {
                    encTemp[i*encB.size()+j] = cc->EvalMult(encA[i], encB[j]);
                    for (size_t kk = 0; (1<<kk) < k; ++kk)
                        encTemp[i*encB.size()+j] = cc->EvalAdd(encTemp[i*encB.size()+j], cc->EvalRotate(encTemp[i*encB.size()+j], (1<<kk)*(j&(1<<kk) ? -1 : 1)));
                    encTemp[i*encB.size()+j] = cc->EvalMult(encTemp[i*encB.size()+j], ptxtMasks[j]);
                }
            } catch (const std::exception& e) {
                newException(e, statsResults, "MatrixMultiplication_" + to_string_precise(reitActual) + ".1[" + to_string_precise(i) + "][" + to_string_precise(j) + "]");
                reiteraciones = reitActual;
                error = true;
            }
    if (error)
        return;
    #pragma omp parallel for shared(encC, reiteraciones, statsResults) schedule(dynamic)
    for (size_t i = 0; i < encA.size(); i++)
        try {
            encC[i] = encTemp[i*encB.size()];
            for (size_t j = 1; j < encB.size(); j++)
                encC[i] += encTemp[i*encB.size()+j];
        } catch (const std::exception& e) {
            newException(e, statsResults, "MatrixMultiplication_" + to_string_precise(reitActual) + ".2[" + to_string_precise(i) + "]");
            reiteraciones = reitActual;
        }
}

template <class Element, typename T = double>
int MatMultRowsXCols(const CryptoContext<Element>& cc,
                    const std::vector<std::vector<T>>& A,
                    const std::vector<std::vector<T>>& B,
                    std::vector<std::vector<std::vector<T>>>& result,
                    size_t reiteraciones,
                    std::map<std::string, std::string>& statsResults,
                    const bool sum = false) {
    #if ENABLE_DEBUG
    std::cout << "Starting Key Generation..." << std::endl;
    #endif
    KeyPair<Element> keys;
    std::set<int> rotations;
    measureBlock("KeyGeneration", statsResults, [&](){
        keys = cc->KeyGen();
        cc->EvalMultKeyGen(keys.secretKey);
        if (sum) {
            cc->EvalSumKeyGen(keys.secretKey);
        } else {
            rotations = rotRowXCols(A[0].size());
            cc->EvalRotateKeyGen(keys.secretKey, std::vector<int>(rotations.begin(), rotations.end()));
        }
    });
    #if ENABLE_DEBUG
    std::cout << "Rotations(RowXCols): ";
    for (int r : rotations)
        std::cout << r << ", ";
    std::cout << std::endl;
    #endif

    #if ENABLE_DEBUG
    std::cout << "Starting Encryption..." << std::endl;
    #endif
    std::vector<std::vector<T>> transB;
    measureBlock("Preprocess", statsResults, [&](){
        transB = std::vector<std::vector<T>>(B[0].size(), std::vector<T>(B.size()));
        #pragma omp parallel for collapse(2) shared(transB) schedule(dynamic)
        for (size_t i = 0; i < B.size(); ++i)
            for (size_t j = 0; j < B[0].size(); ++j)
                transB[j][i] = B[i][j];
    });
    std::vector<Plaintext> ptxA(A.size()), ptxB(transB.size());
    #if USE_OMP_TASKLOOP
    #pragma omp parallel
    {
    #pragma omp single nowait
    {
    #endif
    if constexpr (std::is_floating_point<T>::value) {
        measureBlock("Encode", statsResults, [&](){
            #if USE_OMP_TASKLOOP
            #pragma omp taskgroup
            {
            // Encriptar A
            #pragma omp task
            {
            #pragma omp taskloop shared(ptxA)
            #else
            #pragma omp parallel for shared(ptxA) schedule(dynamic)
            #endif
            for (size_t i = 0; i < A.size(); ++i)
                ptxA[i] = cc->MakeCKKSPackedPlaintext(A[i]);
            #if USE_OMP_TASKLOOP
            }
            #pragma omp task
            {
            #pragma omp taskloop shared(ptxB)
            #else
            #pragma omp parallel for shared(ptxB) schedule(dynamic)
            #endif
            for (size_t i = 0; i < transB.size(); ++i)
                ptxB[i] = cc->MakeCKKSPackedPlaintext(transB[i]);
            #if USE_OMP_TASKLOOP
            }
            }
            #endif
        });
    } else {
        measureBlock("Encode", statsResults, [&](){
            #if USE_OMP_TASKLOOP
            #pragma omp taskgroup
            {
            // Encriptar A
            #pragma omp task
            {
            #pragma omp taskloop shared(ptxA)
            #else
            #pragma omp parallel for shared(ptxA) schedule(dynamic)
            #endif
            for (size_t i = 0; i < A.size(); ++i)
                ptxA[i] = cc->MakePackedPlaintext(A[i]);
            #if USE_OMP_TASKLOOP
            }
            #pragma omp task
            {
            #pragma omp taskloop shared(ptxB)
            #else
            #pragma omp parallel for shared(ptxB) schedule(dynamic)
            #endif
            for (size_t i = 0; i < transB.size(); ++i)
                ptxB[i] = cc->MakePackedPlaintext(transB[i]);
            #if USE_OMP_TASKLOOP
            }
            }
            #endif
        });
    }
    std::vector<Ciphertext<Element>> encA(A.size()), encB(transB.size());
    measureBlock("Encrypt", statsResults, [&](){
        #if USE_OMP_TASKLOOP
        #pragma omp taskgroup
        {
        // Encriptar A
        #pragma omp task
        {
        #pragma omp taskloop shared(encA)
        #else
        #pragma omp parallel for shared(encA) schedule(dynamic)
        #endif
        for (size_t i = 0; i < A.size(); ++i)
            encA[i] = cc->Encrypt(keys.publicKey, ptxA[i]);
        #if USE_OMP_TASKLOOP
        }
        #pragma omp task
        {
        #pragma omp taskloop shared(encB)
        #else
        #pragma omp parallel for shared(encB) schedule(dynamic)
        #endif
        for (size_t i = 0; i < transB.size(); ++i)
            encB[i] = cc->Encrypt(keys.publicKey, ptxB[i]);
        #if USE_OMP_TASKLOOP
        }
        }
        #endif
    });
    addToMap(statsResults, "ctxA_size", to_string_precise(encA.size()));
    addToMap(statsResults, "ctxB_size", to_string_precise(encB.size()));

    std::vector<std::vector<Ciphertext<Element>>> encC(reiteraciones);
    std::vector<std::vector<Plaintext>> ptxC(reiteraciones);
    std::vector<std::vector<std::vector<T>>> decC(reiteraciones);
    for (size_t r = 0; r < reiteraciones; ++r) {
        #if ENABLE_DEBUG
        std::cout << "Starting Matrix Multiplication " << r+1 << "/" << reiteraciones << "..." << std::endl;
        #endif
        if (r > 0) {
            measureBlock("MatrixMultiplication", statsResults, [&](){
                matMultRowsXCols(encC[r-1], encB, encC[r], statsResults, reiteraciones, A[0].size(), r, sum);
            });
        } else {
            measureBlock("MatrixMultiplication", statsResults, [&](){
                matMultRowsXCols(encA, encB, encC[0], statsResults, reiteraciones, A[0].size(), 0, sum);
            });
        }
        if(reiteraciones > r) {
            std::vector<size_t> level = std::vector<size_t>(encA.size());
            #pragma omp parallel for shared(level) schedule(dynamic)
            for (size_t i = 0; i < encA.size(); ++i)
                level[i] = encC[r][i]->GetLevel();
            calculateStats("ctx_level", statsResults, level);
            addToMap(statsResults, "ctxC_size", to_string_precise(encC[r].size()));
            addToMap(statsResults, "ctx_hopLevel", to_string_precise(encC[r][0]->GetHopLevel()));
            addToMap(statsResults, "ctx_noiseScaleDeg", to_string_precise(encC[r][0]->GetNoiseScaleDeg()));
            addToMap(statsResults, "ctx_scalingFactor", to_string_precise(encC[r][0]->GetScalingFactor()));
            measureBlock("Decrypt", statsResults, [&](){
                ptxC[r] = std::vector<Plaintext>(encC[0].size());
                #pragma omp parallel for shared(ptxC, statsResults, reiteraciones) schedule(dynamic)
                for (size_t i = 0; i < encC[0].size(); ++i)
                    try {
                        cc->Decrypt(keys.secretKey, encC[r][i], &ptxC[r][i]);
                    } catch (const std::exception& e) {
                        newException(e, statsResults, "Decrypt_" + to_string_precise(r) + "[" + to_string_precise(i) + "]");
                        reiteraciones = r; // Stop further repetitions
                    }
            });
            if(reiteraciones > r) {
                std::vector<double> logError(ptxC[0].size()),
                                    logPrecision(ptxC[0].size());
                #pragma omp parallel for shared(logError, logPrecision) schedule(dynamic)
                for (size_t i = 0; i < ptxC[0].size(); ++i) {
                    logError[i] = ptxC[r][i]->GetLogError();
                    logPrecision[i] = ptxC[r][i]->GetLogPrecision();
                }
                calculateStats("ptx_logError", statsResults, logError);
                calculateStats("ptx_logPrecision", statsResults, logPrecision);
                addToMap(statsResults, "ptx_length", to_string_precise(ptxC[r][0]->GetLength()));
                if constexpr (std::is_floating_point<T>::value) {
                    measureBlock("Decode", statsResults, [&](){
                        decC[r] = std::vector<std::vector<T>>(encC[0].size());
                        #pragma omp parallel for shared(decC, statsResults, reiteraciones) schedule(dynamic)
                        for (size_t i = 0; i < encC[0].size(); ++i)
                            try {
                                decC[r][i] = ptxC[r][i]->GetRealPackedValue();
                            } catch (const std::exception& e) {
                                newException(e, statsResults, "Decode_" + to_string_precise(r) + "[" + to_string_precise(i) + "]");
                                reiteraciones = r; // Stop further repetitions
                            }
                    });
                } else {
                    measureBlock("Decode", statsResults, [&](){
                        decC[r] = std::vector<std::vector<T>>(encC[0].size());
                        #pragma omp parallel for shared(decC, statsResults, reiteraciones) schedule(dynamic)
                        for (size_t i = 0; i < encC[0].size(); ++i)
                            try {
                                decC[r][i] = ptxC[r][i]->GetPackedValue();
                            } catch (const std::exception& e) {
                                newException(e, statsResults, "Decode_" + to_string_precise(r) + "[" + to_string_precise(i) + "]");
                                reiteraciones = r; // Stop further repetitions
                            }
                    });
                }
                if(reiteraciones > r)
                    measureBlock("Postprocess", statsResults, [&](){
                        #pragma omp parallel for shared(decC) schedule(dynamic)
                        for (size_t i = 0; i < decC[0].size(); ++i)
                            result[r][i] = std::vector<T>(decC[r][i].begin(), decC[r][i].begin() + B[0].size());
                    });
            }
        }
    }
    #if USE_OMP_TASKLOOP
    }
    }
    #endif
    return reiteraciones == 0 ? 2 : 0;
}
