#include <openfhe.h>
#include <omp.h>
#include <map>

using namespace lbcrypto;

std::set<int> rotHE(const usint lado, const usint lvlStrassen=0) {
    #if ENABLE_DEBUG
    printf("RotationsHEMatMult(%d, %d): ", lado, lvlStrassen);
    #endif
    int ladoStr = (lado/(1 << lvlStrassen));
    std::set<int> ret;
    for (int i=1; i<ladoStr; i++) {
        ret.insert(i-ladoStr);
        ret.insert(i);
        ret.insert(i*ladoStr);
    }
    #if ENABLE_DEBUG
    for (const auto& r : ret)
        printf("%d ", r);
    printf("\n");
    #endif
    return ret;
}

template <class Element>
void matMultHE(const Ciphertext<Element>& encA,
               const Ciphertext<Element>& encB,
               Ciphertext<Element>& encC,
               std::map<std::string, std::string>& statsResults,
               size_t& reiteraciones,
               const int lado,
               const size_t reitActual=0,
               const bool fastRot=false,
               const usint lvlStrassen=0) {
    if (lvlStrassen > 0) {
        std::cerr << "HEMatMult does not support Strassen's algorithm" << std::endl;
        addToMap(statsResults, "exceptions", "HEMatMult does not support Strassen's algorithm");
        reiteraciones = reitActual;
        return;
    }
    bool error = false;
    CryptoContext<Element> cc = encA->GetCryptoContext();
    bool isCKKS = (encA->GetEncodingType() == CKKS_PACKED_ENCODING);
    //size_t slots = encA->GetSlots();
    uint nElem = lado*lado;
    Ciphertext<Element> encA0, encB0;
    std::map<std::string, std::shared_ptr<std::vector<Element>>> precomRot;
    if (fastRot) {
        precomRot["encA"] = cc->EvalFastRotationPrecompute(encA);
        precomRot["encB"] = cc->EvalFastRotationPrecompute(encB);
    }

    #if USE_OMP_TASKLOOP
    #pragma omp taskgroup
    {
    #pragma omp task shared(encA0)
    {
    try {
    #endif
        // Step 1-1 encA0
        #if ENABLE_DEBUG
        printf("encA0: START - Step 1-1\n");
        #endif
        if (-lado+1 != 0) {
            if (fastRot) {
                precomRot["encA"] = cc->EvalFastRotationPrecompute(encA);
                encA0 = cc->EvalFastRotation(encA, -lado+1, cc->GetCyclotomicOrder(), precomRot["encA"]);
            } else {
                encA0 = cc->EvalRotate(encA, -lado+1);
            }
        } else {
            encA0 = encA;
            if (fastRot)
                precomRot["encA"] = cc->EvalFastRotationPrecompute(encA);
        }
        encA0 = cc->EvalMult(encA0, makePlaintext(cc, isCKKS, nElem, [lado](int ell){return (unsigned)(ell-2*lado+1) < 1;}));
        #if ENABLE_DEBUG
        printf("encA0: 0/%d\n", 2*lado-1);
        #endif
        #if USE_OMP_TASKLOOP
        #pragma omp taskloop shared(encA0)
        #else
        #pragma omp parallel for shared(encA0) schedule(dynamic)
        #endif
        for (int k=-lado+2; k<lado; k++) {
            try {
                Ciphertext<Element> encA_temp;
                if (k != 0){
                    if (fastRot) {
                        encA_temp = cc->EvalFastRotation(encA, k, cc->GetCyclotomicOrder(), precomRot["encA"]);
                    } else {
                        encA_temp = cc->EvalRotate(encA, k);
                    }
                } else {
                    encA_temp = encA;
                }
                encA_temp = cc->EvalMult(encA_temp, makePlaintext(cc, isCKKS, nElem, [lado, k](int ell){return (k>=0) ? ((unsigned)(ell-lado*k) < (lado-k)) : ((unsigned)(ell-(lado+k)*lado+k) < (lado+k));}));
                #if ENABLE_DEBUG
                printf("encA0: %d/%d\n", k+lado-1, 2*lado-1);
                #endif
                #pragma omp critical
                encA0 += encA_temp;
            } catch (const std::exception& e) {
                newException(e, statsResults, "MatrixMultiplication_" + to_string_precise(reitActual) + ".A1[" + to_string_precise(k) + "]");
                reiteraciones = reitActual;
                error = true;
            }
        }
        #if ENABLE_DEBUG
        printf("encA0: END - Step 1-1\n");
        #endif
    #if USE_OMP_TASKLOOP
    } catch (const std::exception& e) {
    newException(e, statsResults, "MatrixMultiplication_" + to_string_precise(reitActual) + ".A0");
        reiteraciones = reitActual;
        error = true;
    }
    }
    #pragma omp task shared(encB0)
    {
    try {
    #else
        if (error)
            return;
    #endif
        // Step 1-2 encB0
        #if ENABLE_DEBUG
        printf("encB0: START - Step 1-2\n");
        #endif
        if (fastRot) {
            precomRot["encB"] = cc->EvalFastRotationPrecompute(encB);
        }
        encB0 = cc->EvalMult(encB, makePlaintext(cc, isCKKS, nElem, [lado](int ell){return mod(ell, lado) == 0;}));
        #if ENABLE_DEBUG
        printf("encB0: 0/%d\n", lado);
        #endif
        #if USE_OMP_TASKLOOP
        #pragma omp taskloop shared(encB0)
        #else
        #pragma omp parallel for shared(encB0) schedule(dynamic)
        #endif
        for (int k=1; k<lado; k++) {
            try {
                Ciphertext<Element> encB_temp;
                if (fastRot) {
                    encB_temp = cc->EvalFastRotation(encB, k*lado, cc->GetCyclotomicOrder(), precomRot["encB"]);
                } else {
                    encB_temp = cc->EvalRotate(encB, k*lado);
                }
                encB_temp = cc->EvalMult(encB_temp, makePlaintext(cc, isCKKS, nElem, [lado, k](int ell){return mod(ell, lado) == k;}));
                #if ENABLE_DEBUG
                printf("encB0: %d/%d\n", k, lado);
                #endif
                #pragma omp critical
                encB0 += encB_temp;
            } catch (const std::exception& e) {
                newException(e, statsResults, "MatrixMultiplication_" + to_string_precise(reitActual) + ".B1[" + to_string_precise(k) + "]");
                reiteraciones = reitActual;
                error = true;
            }
        }
        #if ENABLE_DEBUG
        printf("encB0: END - Step 1-2\n");
        #endif
    #if USE_OMP_TASKLOOP
    } catch (const std::exception& e) {
    newException(e, statsResults, "MatrixMultiplication_" + to_string_precise(reitActual) + ".B0");
        reiteraciones = reitActual;
        error = true;
    }
    }
    }
    // Step 2
    try {
    #else
        if (error)
            return;
    #endif
        #if ENABLE_DEBUG
        printf("encAB: START - Step 2\n");
        #endif
        encC = cc->EvalMult(encA0, encB0);
        if (fastRot){
            #if USE_OMP_TASKLOOP
            #pragma omp taskgroup
            {
            #pragma omp task shared(precomRot)
            {
            #endif
            precomRot["encA0"] = cc->EvalFastRotationPrecompute(encA0);
            #if USE_OMP_TASKLOOP
            }
            #pragma omp task shared(precomRot)
            {
            #endif
            precomRot["encB0"] = cc->EvalFastRotationPrecompute(encB0);
            #if USE_OMP_TASKLOOP
            }
            }
            #endif
        }
        #if ENABLE_DEBUG
        printf("encAB: 0/%d\n", lado);
        #endif
        #if USE_OMP_TASKLOOP
        #pragma omp taskloop shared(encC)
        #else
        #pragma omp parallel for shared(encC) schedule(dynamic)
        #endif
        for (int k=1; k<lado; k++) {
            try {
                Ciphertext<Element> encAB_temp, encA0_temp1, encA0_temp2, encB0_temp;
                #if USE_OMP_TASKLOOP
                #pragma omp taskgroup
                {
                #pragma omp task shared(encA0_temp1, encA0_temp2)
                {
                #pragma omp taskgroup
                {
                #pragma omp task shared(encA0_temp1)
                {
                #endif
                if (fastRot) {
                        encA0_temp1 = cc->EvalFastRotation(encA0, k, cc->GetCyclotomicOrder(), precomRot["encA0"]);
                } else {
                        encA0_temp1 = cc->EvalRotate(encA0, k);
                }
                encA0_temp1 = cc->EvalMult(encA0_temp1, makePlaintext(cc, isCKKS, nElem, [lado, k](int ell){return (unsigned)(mod(ell, lado)) < lado-k;}));
                #if USE_OMP_TASKLOOP
                }
                #pragma omp task shared(encA0_temp2)
                {
                #endif
                if (fastRot) {
                        encA0_temp2 = cc->EvalFastRotation(encA0, k-lado, cc->GetCyclotomicOrder(), precomRot["encA0"]);
                } else {
                        encA0_temp2 = cc->EvalRotate(encA0, k-lado);
                }
                encA0_temp2 = cc->EvalMult(encA0_temp2, makePlaintext(cc, isCKKS, nElem, [lado, k](int ell){return (unsigned)(mod(ell, lado)-lado+k) < k;}));
                #if USE_OMP_TASKLOOP
                }
                }
                #endif
                encA0_temp1 = cc->EvalAdd(encA0_temp1, encA0_temp2);
                #if USE_OMP_TASKLOOP
                }
                #pragma omp task shared(encB0_temp)
                {
                #endif
                if (fastRot) {
                        encB0_temp = cc->EvalFastRotation(encB0, lado*k, cc->GetCyclotomicOrder(), precomRot["encB0"]);
                } else {
                        encB0_temp = cc->EvalRotate(encB0, lado*k);
                }
                #if USE_OMP_TASKLOOP
                }
                }
                #endif
                encAB_temp = cc->EvalMult(encA0_temp1, encB0_temp);
                #if ENABLE_DEBUG
                printf("encAB: %d/%d\n", k, lado);
                #endif
                #pragma omp critical
                try {
                    encC += encAB_temp;
                } catch (const std::exception& e) {
                    newException(e, statsResults, "MatrixMultiplication_" + to_string_precise(reitActual) + ".AB1[" + to_string_precise(k) + "]");
                    reiteraciones = reitActual;
                    error = true;
                }
            } catch (const std::exception& e) {
                newException(e, statsResults, "MatrixMultiplication_" + to_string_precise(reitActual) + ".AB1[" + to_string_precise(k) + "]");
                reiteraciones = reitActual;
                error = true;
            }
        }
        #if ENABLE_DEBUG
        printf("ctAB: END - Step 2\n");
        #endif
    #if USE_OMP_TASKLOOP
    } catch (const std::exception& e) {
    newException(e, statsResults, "MatrixMultiplication_" + to_string_precise(reitActual) + ".AB");
        reiteraciones = reitActual;
        error = true;
    }
    #endif
}

template <class Element, typename T>
int MatMultHE(const CryptoContext<Element>& cc,
              const std::vector<std::vector<T>>& A,
              const std::vector<std::vector<T>>& B,
              std::vector<std::vector<std::vector<T>>>& result,
              size_t reiteraciones,
              std::map<std::string, std::string>& statsResults,
              const bool fastRot=false) {
    size_t lado = 1 << (int)ceil(log2(std::max({A.size(), A[0].size(), B.size(), B[0].size()})));
    #if ENABLE_DEBUG
    std::cout << "Starting Key Generation..." << std::endl;
    #endif
    KeyPair<Element> keys;
    std::set<int> rotations;
    measureBlock("KeyGeneration", statsResults, [&](){
        keys = cc->KeyGen();
        cc->EvalMultKeyGen(keys.secretKey);
        rotations = rotHE(lado);
        cc->EvalRotateKeyGen(keys.secretKey, std::vector<int>(rotations.begin(), rotations.end()));
    });

    #if ENABLE_DEBUG
    std::cout << "Starting Encryption..." << std::endl;
    #endif
    std::vector<T> flatA, flatB;
    measureBlock("Preprocess", statsResults, [&](){
        if (lado == A.size() && lado == A[0].size()) {
            flatA = flatten(A);
        } else {
            flatA = flatten(padMatrix(A, lado));
        }
        if (lado == B.size() && lado == B[0].size()) {
            flatB = flatten(B);
        } else {
            flatB = flatten(padMatrix(B, lado));
        }
    });
    Plaintext ptxA, ptxB;
    if constexpr (std::is_floating_point<T>::value) {
        measureBlock("Encode", statsResults, [&](){
            ptxA = cc->MakeCKKSPackedPlaintext(flatA);
            ptxB = cc->MakeCKKSPackedPlaintext(flatB);
        });
    } else {
        measureBlock("Encode", statsResults, [&](){
            ptxA = cc->MakePackedPlaintext(flatA);
            ptxB = cc->MakePackedPlaintext(flatB);
        });
    }
    Ciphertext<Element> encA, encB;
    measureBlock("Encrypt", statsResults, [&](){
        encA = cc->Encrypt(keys.publicKey, ptxA);
        encB = cc->Encrypt(keys.publicKey, ptxB);
    });
    addToMap(statsResults, "ctxA_size", to_string_precise(1));
    addToMap(statsResults, "ctxB_size", to_string_precise(1));

    std::vector<Ciphertext<Element>> encC(reiteraciones);
    std::vector<Plaintext> ptxC(reiteraciones);
    std::vector<std::vector<T>> decC(reiteraciones);
    #if USE_OMP_TASKLOOP
    #pragma omp parallel
    {
    #pragma omp single nowait
    {
    #endif
    for (size_t r = 0; r < reiteraciones; ++r) {
        #if ENABLE_DEBUG
        std::cout << "Starting Matrix Multiplication " << r+1 << "/" << reiteraciones << "..." << std::endl;
        #endif
        if (r > 0) {
            measureBlock("MatrixMultiplication", statsResults, [&](){
                matMultHE(encC[r-1], encB, encC[r], statsResults, reiteraciones, lado, r, fastRot);
            });
        } else {
            measureBlock("MatrixMultiplication", statsResults, [&](){
                matMultHE(encA, encB, encC[0], statsResults, reiteraciones, lado, r, fastRot);
            });
        }
        if(reiteraciones > r) {
            #if ENABLE_DEBUG
            std::cout << "Matrix Multiplication " << r+1 << "/" << reiteraciones << " completed." << std::endl;
            #endif
            calculateStats("ctx_level", statsResults, std::vector<size_t>{encC[r]->GetLevel()});
            addToMap(statsResults, "ctxC_size", to_string_precise(1));
            addToMap(statsResults, "ctx_hopLevel", to_string_precise(encC[r]->GetHopLevel()));
            addToMap(statsResults, "ctx_noiseScaleDeg", to_string_precise(encC[r]->GetNoiseScaleDeg()));
            addToMap(statsResults, "ctx_scalingFactor", to_string_precise(encC[r]->GetScalingFactor()));
            measureBlock("Decrypt", statsResults, [&](){
                try {
                    cc->Decrypt(keys.secretKey, encC[r], &ptxC[r]);
                } catch (const std::exception& e) {
                    newException(e, statsResults, "Decrypt_" + to_string_precise(r));
                    reiteraciones = r; // Stop further repetitions
                }
            });
            if(reiteraciones > r) {
                calculateStats("ptx_logError", statsResults, std::vector<double>{ptxC[r]->GetLogError()});
                calculateStats("ptx_logPrecision", statsResults, std::vector<double>{ptxC[r]->GetLogPrecision()});
                addToMap(statsResults, "ptx_length", to_string_precise(ptxC[r]->GetLength()));
                if constexpr (std::is_floating_point<T>::value) {
                    measureBlock("Decode", statsResults, [&](){
                        try {
                            decC[r] = ptxC[r]->GetRealPackedValue();
                        } catch (const std::exception& e) {
                            newException(e, statsResults, "Decode_" + to_string_precise(r));
                            reiteraciones = r; // Stop further repetitions
                        }
                    });
                } else {
                    measureBlock("Decode", statsResults, [&](){
                        try {
                            decC[r] = ptxC[r]->GetPackedValue();
                        } catch (const std::exception& e) {
                            newException(e, statsResults, "Decode_" + to_string_precise(r));
                            reiteraciones = r; // Stop further repetitions
                        }
                    });
                }
                if(reiteraciones > r)
                    measureBlock("Postprocess", statsResults, [&](){
                        #pragma omp parallel for shared(result) schedule(dynamic)
                        for (size_t i = 0; i < A.size(); ++i)
                            result[r][i] = std::vector<T>(decC[r].begin()+i*lado, decC[r].begin()+i*lado+B[0].size());
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