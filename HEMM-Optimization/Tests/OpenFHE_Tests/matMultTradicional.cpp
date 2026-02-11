#include <openfhe.h>
#include <omp.h>
#include <map>

using namespace lbcrypto;

template <class Element>
void matMultTradicional(const std::vector<std::vector<Ciphertext<Element>>>& encA,
                        const std::vector<std::vector<Ciphertext<Element>>>& encB,
                        std::vector<std::vector<Ciphertext<Element>>>& encC,
                        std::map<std::string, std::string>& statsResults,
                        size_t& reiteraciones,
                        const size_t reitActual=0) {
    CryptoContext<Element> cc = encA[0][0]->GetCryptoContext();
    encC = std::vector<std::vector<Ciphertext<Element>>>(encA.size(), std::vector<Ciphertext<Element>>(encB[0].size()));
    #pragma omp parallel for collapse(2) shared(encC, statsResults, reiteraciones) schedule(dynamic)
    for (size_t i = 0; i < encA.size(); ++i)
        for (size_t j = 0; j < encB[0].size(); ++j)
            try {
                encC[i][j] = cc->EvalMult(encA[i][0], encB[0][j]);
                for (size_t k = 1; k < encA[0].size(); ++k)
                    encC[i][j] = cc->EvalAdd(encC[i][j], cc->EvalMult(encA[i][k], encB[k][j]));
            } catch (const std::exception& e) {
                newException(e, statsResults, "MatrixMultiplication_" + to_string_precise(reitActual) + "[" + to_string_precise(i) + "][" + to_string_precise(j) + "]");
                reiteraciones = reitActual; // Stop further repetitions
            }
}

template <class Element, typename T = double>
int MatMultTradicional(const CryptoContext<Element>& cc,
                       const std::vector<std::vector<T>>& A,
                       const std::vector<std::vector<T>>& B,
                       std::vector<std::vector<std::vector<T>>>& result,
                       size_t reiteraciones,
                       std::map<std::string, std::string>& statsResults) {
    #if ENABLE_DEBUG
    std::cout << "Starting Key Generation..." << std::endl;
    #endif
    KeyPair<Element> keys;
    measureBlock("KeyGeneration", statsResults, [&](){
        keys = cc->KeyGen();
        cc->EvalMultKeyGen(keys.secretKey);
    });

    #if ENABLE_DEBUG
    std::cout << "Starting Encryption..." << std::endl;
    #endif
    std::vector<std::vector<Ciphertext<Element>>> encA(A.size(), std::vector<Ciphertext<Element>>(A[0].size())), encB(B.size(), std::vector<Ciphertext<Element>>(B[0].size()));
    #if USE_OMP_TASKLOOP
    #pragma omp parallel
    {
    #pragma omp single nowait
    {
    #endif
    std::vector<std::vector<std::vector<T>>> preA(A.size(), std::vector<std::vector<T>>(A[0].size(), std::vector<T>(1))), preB(B.size(), std::vector<std::vector<T>>(B[0].size(), std::vector<T>(1)));
    measureBlock("Preprocess", statsResults, [&](){
        #if USE_OMP_TASKLOOP
        #pragma omp taskgroup
        {
        #pragma omp task
        {
        #pragma omp taskloop collapse(2) shared(preA)
        #else
        #pragma omp parallel for collapse(2) shared(preA) schedule(dynamic)
        #endif
        for (size_t i = 0; i < A.size(); ++i)
            for (size_t j = 0; j < A[0].size(); ++j)
                preA[i][j][0] = A[i][j];
        #if USE_OMP_TASKLOOP
        }
        #pragma omp task
        {
        #pragma omp taskloop collapse(2) shared(preB)
        #else
        #pragma omp parallel for collapse(2) shared(preB) schedule(dynamic)
        #endif
        for (size_t i = 0; i < B.size(); ++i)
            for (size_t j = 0; j < B[0].size(); ++j)
                preB[i][j][0] = B[i][j];
        #if USE_OMP_TASKLOOP
        }
        }
        #endif
    });
    std::vector<std::vector<Plaintext>> ptxA(A.size(), std::vector<Plaintext>(A[0].size())), ptxB(B.size(), std::vector<Plaintext>(B[0].size()));
    if constexpr (std::is_floating_point<T>::value) {
        measureBlock("Encode", statsResults, [&](){
            #if USE_OMP_TASKLOOP
            #pragma omp taskgroup
            {
            #pragma omp task
            {
            // Encriptar A
            #pragma omp taskloop collapse(2) shared(ptxA)
            #else
            #pragma omp parallel for collapse(2) shared(ptxA) schedule(dynamic)
            #endif
            for (size_t i = 0; i < A.size(); ++i)
                for (size_t j = 0; j < A[0].size(); ++j)
                    ptxA[i][j] = cc->MakeCKKSPackedPlaintext(preA[i][j]);
            #if USE_OMP_TASKLOOP
            }
            #pragma omp task
            {
            #pragma omp taskloop collapse(2) shared(ptxB)
            #else
            #pragma omp parallel for collapse(2) shared(ptxB) schedule(dynamic)
            #endif
            for (size_t i = 0; i < B.size(); ++i)
                for (size_t j = 0; j < B[0].size(); ++j)
                    ptxB[i][j] = cc->MakeCKKSPackedPlaintext(preB[i][j]);
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
            #pragma omp task
            {
            // Encriptar A
            #pragma omp taskloop collapse(2) shared(ptxA)
            #else
            #pragma omp parallel for collapse(2) shared(ptxA) schedule(dynamic)
            #endif
            for (size_t i = 0; i < A.size(); ++i)
                for (size_t j = 0; j < A[0].size(); ++j)
                    ptxA[i][j] = cc->MakePackedPlaintext(A[i][j]);
            #if USE_OMP_TASKLOOP
            }
            #pragma omp task
            {
            #pragma omp taskloop collapse(2) shared(ptxB)
            #else
            #pragma omp parallel for collapse(2) shared(ptxB) schedule(dynamic)
            #endif
            for (size_t i = 0; i < B.size(); ++i)
                for (size_t j = 0; j < B[0].size(); ++j)
                    ptxB[i][j] = cc->MakePackedPlaintext(B[i][j]);
            #if USE_OMP_TASKLOOP
            }
            }
            #endif
        });
    }
    measureBlock("Encrypt", statsResults, [&](){
        #if USE_OMP_TASKLOOP
        #pragma omp taskgroup
        {
        #pragma omp task
        {
        // Encriptar A
        #pragma omp taskloop collapse(2) shared(encA)
        #else
        #pragma omp parallel for collapse(2) shared(encA) schedule(dynamic)
        #endif
        for (size_t i = 0; i < A.size(); ++i)
            for (size_t j = 0; j < A[0].size(); ++j)
                encA[i][j] = cc->Encrypt(keys.publicKey, ptxA[i][j]);
        #if USE_OMP_TASKLOOP
        }
        #pragma omp task
        {
        #pragma omp taskloop collapse(2) shared(encB)
        #else
        #pragma omp parallel for collapse(2) shared(encB) schedule(dynamic)
        #endif
        for (size_t i = 0; i < B.size(); ++i)
            for (size_t j = 0; j < B[0].size(); ++j)
                encB[i][j] = cc->Encrypt(keys.publicKey, ptxB[i][j]);
        #if USE_OMP_TASKLOOP
        }
        }
        #endif
    });
    #if USE_OMP_TASKLOOP
    }
    }
    #endif
    addToMap(statsResults, "ctxA_size", to_string_precise(encA.size() * encA[0].size()));
    addToMap(statsResults, "ctxB_size", to_string_precise(encB.size() * encB[0].size()));

    std::vector<std::vector<std::vector<Ciphertext<Element>>>> encC(reiteraciones);
    std::vector<std::vector<std::vector<Plaintext>>> ptxC(reiteraciones);
    std::vector<std::vector<std::vector<std::vector<T>>>> decC(reiteraciones);
    for (size_t r = 0; r < reiteraciones; ++r) {
        #if ENABLE_DEBUG
        std::cout << "Starting Matrix Multiplication " << r+1 << "/" << reiteraciones << "..." << std::endl;
        #endif
        if (r > 0) {
            measureBlock("MatrixMultiplication", statsResults, [&](){
                matMultTradicional(encC[r-1], encB, encC[r], statsResults, reiteraciones, r);
            });
        } else {
            measureBlock("MatrixMultiplication", statsResults, [&](){
                matMultTradicional(encA, encB, encC[0], statsResults, reiteraciones);
            });
        }
        if(reiteraciones > r) {
            std::vector<size_t> level = std::vector<size_t>(encA.size() * encB[0].size());
            #pragma omp parallel for collapse(2) shared(level) schedule(dynamic)
            for (size_t i = 0; i < encA.size(); ++i)
                for (size_t j = 0; j < encB[0].size(); ++j)
                    level[i * encB[0].size() + j] = encC[r][i][j]->GetLevel();
            calculateStats("ctx_level", statsResults, level);
            addToMap(statsResults, "ctxC_size", to_string_precise(encC[r].size() * encC[r][0].size()));
            addToMap(statsResults, "ctx_hopLevel", to_string_precise(encC[r][0][0]->GetHopLevel()));
            addToMap(statsResults, "ctx_noiseScaleDeg", to_string_precise(encC[r][0][0]->GetNoiseScaleDeg()));
            addToMap(statsResults, "ctx_scalingFactor", to_string_precise(encC[r][0][0]->GetScalingFactor()));
            measureBlock("Decrypt", statsResults, [&](){
                ptxC[r] = std::vector<std::vector<Plaintext>>(encC[0].size(), std::vector<Plaintext>(encC[0][0].size()));
                #pragma omp parallel for collapse(2) shared(ptxC, statsResults, reiteraciones) schedule(dynamic)
                for (size_t i = 0; i < encC[0].size(); ++i)
                    for (size_t j = 0; j < encC[0][0].size(); ++j)
                        try {
                            cc->Decrypt(keys.secretKey, encC[r][i][j], &ptxC[r][i][j]);
                        } catch (const std::exception& e) {
                            newException(e, statsResults, "Decrypt_" + to_string_precise(r) + "[" + to_string_precise(i) + "][" + to_string_precise(j) + "]");
                            reiteraciones = r; // Stop further repetitions
                        }
            });
            if(reiteraciones > r) {
                std::vector<double> logError(ptxC[0].size() * ptxC[0][0].size()),
                                    logPrecision(ptxC[0].size() * ptxC[0][0].size());
                #pragma omp parallel for collapse(2) shared(logError, logPrecision) schedule(dynamic)
                for (size_t i = 0; i < ptxC[0].size(); ++i)
                    for (size_t j = 0; j < ptxC[0][0].size(); ++j) {
                        logError[i * ptxC[0][0].size() + j] = ptxC[r][i][j]->GetLogError();
                        logPrecision[i * ptxC[0][0].size() + j] = ptxC[r][i][j]->GetLogPrecision();
                    }
                calculateStats("ptx_logError", statsResults, logError);
                calculateStats("ptx_logPrecision", statsResults, logPrecision);
                addToMap(statsResults, "ptx_length", to_string_precise(ptxC[r][0][0]->GetLength()));
                if constexpr (std::is_floating_point<T>::value) {
                    measureBlock("Decode", statsResults, [&](){
                        decC[r] = std::vector<std::vector<std::vector<T>>>(encC[0].size(), std::vector<std::vector<T>>(encC[0][0].size()));
                        #pragma omp parallel for collapse(2) shared(decC, statsResults, reiteraciones) schedule(dynamic)
                        for (size_t i = 0; i < encC[0].size(); ++i)
                            for (size_t j = 0; j < encC[0][0].size(); ++j)
                                try {
                                    decC[r][i][j] = ptxC[r][i][j]->GetRealPackedValue();
                                } catch (const std::exception& e) {
                                    newException(e, statsResults, "Decode_" + to_string_precise(r) + "[" + to_string_precise(i) + "][" + to_string_precise(j) + "]");
                                    reiteraciones = r; // Stop further repetitions
                                }
                    });
                } else {
                    measureBlock("Decode", statsResults, [&](){
                        decC[r] = std::vector<std::vector<std::vector<T>>>(encC[0].size(), std::vector<std::vector<T>>(encC[0][0].size()));
                        #pragma omp parallel for collapse(2) shared(decC, statsResults, reiteraciones) schedule(dynamic)
                        for (size_t i = 0; i < encC[0].size(); ++i)
                            for (size_t j = 0; j < encC[0][0].size(); ++j)
                                try {
                                    decC[r][i][j] = ptxC[r][i][j]->GetPackedValue();
                                } catch (const std::exception& e) {
                                    newException(e, statsResults, "Decode_" + to_string_precise(r) + "[" + to_string_precise(i) + "][" + to_string_precise(j) + "]");
                                    reiteraciones = r; // Stop further repetitions
                                }
                    });
                }
                if (reiteraciones > r)
                    measureBlock("Postprocess", statsResults, [&](){
                        #pragma omp parallel for collapse(2) shared(result, statsResults, reiteraciones) schedule(dynamic)
                        for (size_t i = 0; i < A.size(); ++i)
                            for (size_t j = 0; j < B[0].size(); ++j)
                                result[r][i][j] = decC[r][i][j][0];
                    });
            }
        }
    }
    
    return reiteraciones == 0 ? 2 : 0;
}
