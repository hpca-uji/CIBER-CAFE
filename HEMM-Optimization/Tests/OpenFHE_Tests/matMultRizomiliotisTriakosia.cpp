#include <openfhe.h>
#include <omp.h>
#include <map>

using namespace lbcrypto;


std::set<int> rotRizomiliotisTriakosia(const usint lado,
                                       const bool sumRowsCols,
                                       const usint lvlStrassen=0) {
    #if ENABLE_DEBUG
    printf("RotationsRizomiliotisTriakosia(%d, %d): ", lado, lvlStrassen);
    #endif
    int ladoStr = lado/(1 << lvlStrassen),
        l=floor(log2(lado/ladoStr));
    std::set<int> ret;
    for (int i=1; i<ladoStr; i++) {
        ret.insert(i-(mod(i*ladoStr, lado)+floor(i*ladoStr/lado)*lado*ladoStr));
        ret.insert(i*lado-(mod(i*ladoStr, lado)+floor(i*ladoStr/lado)*lado*ladoStr));
    }
    for (int i=0; (1<<i)<ladoStr; i++) {
        if (!sumRowsCols)
            ret.insert(-(1 << i)); // OK
        ret.insert(-(1 << i)*lado);
    }
    for(int i=0; i < lvlStrassen && (1<<i) < ladoStr; i++)
        ret.insert(ladoStr << i);
    if (!sumRowsCols)
        for (int i=0; (1<<(i+lvlStrassen)) < ladoStr; i++)
            ret.insert((lado*ladoStr) << i);
    #if ENABLE_DEBUG
    for (const auto& r : ret)
        printf("%d ", r);
    printf("\n");
    #endif
    return ret;
}

template <class Element>
int matMultRizomiliotisTriakosia(const Ciphertext<Element>& encA,
               const Ciphertext<Element>& encB,
               Ciphertext<Element>& encC,
               std::map<std::string, std::string>& statsResults,
               size_t& reiteraciones,
               const int lado,
               const size_t reitActual=0,
               const std::shared_ptr<std::map<uint32_t, EvalKey<Element>>> sumRowsColsKeys = nullptr,
               const bool fastRot=false,
               const usint lvlStrassen=0) {
    bool error = false;
    CryptoContext<Element> cc = encA->GetCryptoContext();
    bool isCKKS = (encA->GetEncodingType() == CKKS_PACKED_ENCODING);
    int ladoStr = lado >> lvlStrassen,
        slotsMax = cc->GetRingDimension()/2;
    uint32_t nSlotsStr = std::max(ladoStr*ladoStr*ladoStr, lado*lado);
    if (nSlotsStr > slotsMax) {
        newException(statsResults, "Number of slots max. (" + std::to_string(slotsMax) + ") is less than required for Rizomiliotis-Triakosia method (" + std::to_string(nSlotsStr) + ").");
        reiteraciones = reitActual;
        return 3;
    }
    Ciphertext<Element> encA0 = encA->Clone(),
                        encAT,
                        encB0 = encB->Clone(),
                        encBT;
    #if USE_OMP_TASKLOOP
    #pragma omp taskgroup
    {
    #pragma omp task shared(encAT)
    {
    #endif
    try {
        #if ENABLE_DEBUG
        printf("encA: START\n");
        #endif
        if (isCKKS && nSlotsStr > encA->GetSlots()) {
            encA0->SetSlots(nSlotsStr);
            #if ENABLE_DEBUG
            printf("encA: Resize(%u -> %d)\n", encA->GetSlots(), nSlotsStr);
            #endif
        }
        encAT = cc->EvalMult(
            encA0,
            makePlaintext(cc, isCKKS, nSlotsStr, [ladoStr, lado](int j){return j<ladoStr*lado && mod(j, ladoStr)==0 && mod(floor(j/ladoStr), lado/ladoStr)==0;}, nSlotsStr)
        );
        #if ENABLE_DEBUG
        printf("encA: 0/%d\n", ladoStr);
        #endif
        std::shared_ptr<std::vector<Element>> precompEncA0;
        if (fastRot) {
            precompEncA0 = cc->EvalFastRotationPrecompute(encA0);
        }
        #if USE_OMP_TASKLOOP
        #pragma omp taskloop shared(encAT)
        #else
        #pragma omp parallel for shared(encAT) schedule(dynamic)
        #endif
        for (int i=1; i<ladoStr; i++) {
            try {
                Ciphertext<Element> encA_temp;
                int rot = i-(mod(i*ladoStr, lado)+floor(i*ladoStr/lado)*lado*ladoStr);
                if (rot != 0) {
                    if (fastRot) {
                        encA_temp = cc->EvalFastRotation(encA0, rot, cc->GetCyclotomicOrder(), precompEncA0);
                    } else {
                        encA_temp = cc->EvalRotate(encA0, rot);
                    }
                } else {
                    encA_temp = encA0;
                }
                encA_temp = cc->EvalMult(
                    encA_temp,
                    makePlaintext(cc, isCKKS, nSlotsStr, [ladoStr, lado, i](int j){return j>=ladoStr*lado*floor(i*ladoStr/lado) && j<ladoStr*lado*(floor(i*ladoStr/lado)+1) && mod(j, ladoStr)==0 && mod(floor(j/ladoStr), lado/ladoStr)==mod(i, lado/ladoStr);}, nSlotsStr)
                );
                
                #if ENABLE_DEBUG
                printf("encA: %d/%d -> rot(%d)\n", i, ladoStr, rot);
                #endif
                #pragma omp critical
                try {
                    encAT += encA_temp;
                } catch (const std::exception& e) {
                    newException(e, statsResults, "MatrixMultiplication_" + to_string_precise(reitActual) + ".A2[" + to_string_precise(i) + "]");
                    reiteraciones = reitActual;
                    error = true;
                }
            } catch (const std::exception& e) {
                newException(e, statsResults, "MatrixMultiplication_" + to_string_precise(reitActual) + ".A1[" + to_string_precise(i) + "]");
                reiteraciones = reitActual;
                error = true;
            }
        }
        if (sumRowsColsKeys != nullptr) {
            encAT = cc->EvalSumCols(encAT, ladoStr, *sumRowsColsKeys);
        } else {
            for (int i=0; (1 << i)<ladoStr; i++) {
                encAT += cc->EvalRotate(encAT, -(1 << i));
                #if ENABLE_DEBUG
                printf("encA: Expand(%d/%d)->%d\n", i, (int)log2(ladoStr), -(1 << i));
                #endif
            }

        }
        #if ENABLE_DEBUG
        printf("encA: END\n");
        #endif
    } catch (const std::exception& e) {
        newException(e, statsResults, "MatrixMultiplication_" + to_string_precise(reitActual) + ".A0");
        reiteraciones = reitActual;
        error = true;
    }
    #if USE_OMP_TASKLOOP
    }
    #pragma omp task shared(encBT)
    {
    #else
    if (error)
        return 1;
    #endif
    try {
        #if ENABLE_DEBUG
        printf("encB: START\n");
        #endif
        if (nSlotsStr > encB->GetSlots()) {
            encB0->SetSlots(nSlotsStr);
            #if ENABLE_DEBUG
            printf("encB: Resize(%u -> %d)\n", encB->GetSlots(), nSlotsStr);
            #endif
        }
        encBT = cc->EvalMult(
            encB0,
            makePlaintext(cc, isCKKS, nSlotsStr, [ladoStr, lado](int j){return j<ladoStr;}, nSlotsStr)
        );
        #if ENABLE_DEBUG
        printf("encB: 0/%d\n", ladoStr);
        #endif
        std::shared_ptr<std::vector<Element>> precompEncB0;
        if (fastRot) {
            precompEncB0 = cc->EvalFastRotationPrecompute(encB0);
        }
        #if USE_OMP_TASKLOOP
        #pragma omp taskloop shared(encBT)
        #else
        #pragma omp parallel for shared(encBT) schedule(dynamic)
        #endif
        for (int i=1; i<ladoStr; i++) {
            try {
                Ciphertext<Element> encB_temp;
                int rot = i*lado-(mod(i*ladoStr, lado)+floor(i*ladoStr/lado)*lado*ladoStr);
                //int rot = -i*lado*(ladoStr-1);
                if (rot != 0) {
                    if (fastRot) {
                        encB_temp = cc->EvalFastRotation(encB0, rot, cc->GetCyclotomicOrder(), precompEncB0);
                    } else {
                        encB_temp = cc->EvalRotate(encB0, rot);
                    }
                } else {
                    encB_temp = encB0;
                }
                encB_temp = cc->EvalMult(
                    encB_temp,
                    makePlaintext(cc, isCKKS, nSlotsStr, [ladoStr, lado, i](int j){return j>=(mod(i, lado/ladoStr)+lado*floor(i*ladoStr/lado))*ladoStr && j<(mod(i, lado/ladoStr)+lado*floor(i*ladoStr/lado)+1)*ladoStr;}, nSlotsStr)
                );
                #if ENABLE_DEBUG
                printf("encB: %d/%d -> rot(%d)\n", i, ladoStr, rot);
                #endif
                #pragma omp critical
                try {
                    encBT += encB_temp;
                } catch (const std::exception& e) {
                    newException(e, statsResults, "MatrixMultiplication_" + to_string_precise(reitActual) + ".B2[" + to_string_precise(i) + "]");
                    reiteraciones = reitActual;
                    error = true;
                }
            } catch (const std::exception& e) {
                newException(e, statsResults, "MatrixMultiplication_" + to_string_precise(reitActual) + ".B1[" + to_string_precise(i) + "]");
                reiteraciones = reitActual;
                error = true;
            }
        }
        for (int i=0; (1 << i)<ladoStr; i++) {
            encBT += cc->EvalRotate(encBT, -(1 << i)*lado);
            #if ENABLE_DEBUG
            printf("encB: Expand(%d/%d)->%d\n", i, (int)log2(ladoStr), -(1 << i));
            #endif
        }
        #if ENABLE_DEBUG
        printf("encB: END\n");
        #endif
    } catch (const std::exception& e) {
        newException(e, statsResults, "MatrixMultiplication_" + to_string_precise(reitActual) + ".B0");
        reiteraciones = reitActual;
        error = true;
    }
    #if USE_OMP_TASKLOOP
    }
    }
    #endif
    if (error)
        return 1;
    try {
        // Multiplica A & B
        #if ENABLE_DEBUG
        printf("encC: Prepared(encAT(lvl=%zu, noise=%zu), encBT(lvl=%zu, noise=%zu))\n", encAT->GetLevel(), encAT->GetNoiseScaleDeg(), encBT->GetLevel(), encBT->GetNoiseScaleDeg());
        #endif
        encC = cc->EvalMult(encAT, encBT);
        #if ENABLE_DEBUG
        printf("encC: Calculated lvl=%zu, noise=%zu\n", encC->GetLevel(), encC->GetNoiseScaleDeg());
        #endif
        for(int i=0; i < lvlStrassen && (1<<i) < ladoStr; i++) {
            encC += cc->EvalRotate(encC, ladoStr << i);
            #if ENABLE_DEBUG
            printf("encC: Rotate.Row(%d)\n", ladoStr << i);
            #endif
        }
        if (sumRowsColsKeys != nullptr) {
            encC = cc->EvalSumRows(encC, lado*ladoStr, *sumRowsColsKeys, nSlotsStr*4);
            #if ENABLE_DEBUG
            printf("encC: SumRows.Cols(%d; %d)\n", lado*ladoStr, nSlotsStr*4);
            #endif
        } else {
            for (int i=0; (1<<(i+lvlStrassen)) < ladoStr; i++) {
                encC += cc->EvalRotate(encC, (lado*ladoStr) << i);
                #if ENABLE_DEBUG
                printf("encC: Rotate.Col(%d)\n", (lado*ladoStr) << i);
                #endif
            }
        }
        if (nSlotsStr > encA->GetSlots()) {
            encC->SetSlots(encA->GetSlots());
            #if ENABLE_DEBUG
            printf("encC: Resize(%d -> %u)\n", nSlotsStr, encA->GetSlots());
            #endif
        }
        #if ENABLE_DEBUG
        printf("encC: END slots=%u lvl=%zu, noise=%zu\n", encC->GetSlots(), encC->GetLevel(), encC->GetNoiseScaleDeg());
        #endif
    } catch (const std::exception& e) {
        newException(e, statsResults, "MatrixMultiplication_" + to_string_precise(reitActual) + ".C");
        reiteraciones = reitActual;
        error = true;
    }
    if (error) {
        return 1;
    } else {
        return 0;
    }
}

template <class Element, typename T>
int MatMultRizomiliotisTriakosia(const CryptoContext<Element>& cc,
                                 const std::vector<std::vector<T>>& A,
                                 const std::vector<std::vector<T>>& B,
                                 std::vector<std::vector<std::vector<T>>>& result,
                                 size_t reiteraciones,
                                 std::map<std::string, std::string>& statsResults,
                                 const bool sumRowsCols=false,
                                 const bool fastRot=false) {
    size_t lado = 1 << (int)ceil(log2(std::max({A.size(), A[0].size(), B.size(), B[0].size()})));
    if (lado*lado*lado > cc->GetRingDimension()/2) {
        newException(statsResults, "Number of slots max. (" + std::to_string(cc->GetRingDimension()/2) + ") is less than required for Rizomiliotis-Triakosia method (" + std::to_string(lado*lado*lado) + ").");
        return 3;
    }
    #if ENABLE_DEBUG
    std::cout << "Starting Key Generation..." << std::endl;
    #endif
    KeyPair<Element> keys;
    std::set<int> rotations;
    std::shared_ptr<std::map<uint32_t, EvalKey<Element>>> sumRowsColsKeys;
    measureBlock("KeyGeneration", statsResults, [&](){
        keys = cc->KeyGen();
        cc->EvalMultKeyGen(keys.secretKey);
        rotations = rotRizomiliotisTriakosia(lado, sumRowsCols);
        cc->EvalRotateKeyGen(keys.secretKey, std::vector<int>(rotations.begin(), rotations.end()));
        if (sumRowsCols) {
            sumRowsColsKeys = cc->EvalSumColsKeyGen(keys.secretKey);
            auto rowsKeys = cc->EvalSumRowsKeyGen(keys.secretKey, nullptr, lado*lado, lado*lado*lado*4);
            sumRowsColsKeys->insert(rowsKeys->begin(), rowsKeys->end());
        }
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
    int ret = 0;
    #if USE_OMP_TASKLOOP
    #pragma omp parallel
    {
    #pragma omp single nowait
    {
    #endif
    for (size_t r = 0; r < reiteraciones && ret == 0; ++r) {
        #if ENABLE_DEBUG
        std::cout << "Starting Matrix Multiplication " << r+1 << "/" << reiteraciones << "..." << std::endl;
        #endif
        if (r > 0) {
            measureBlock("MatrixMultiplication", statsResults, [&](){
                ret = matMultRizomiliotisTriakosia(encC[r-1], encB, encC[r], statsResults, reiteraciones, lado, r, sumRowsColsKeys, fastRot);
            });
        } else {
            measureBlock("MatrixMultiplication", statsResults, [&](){
                ret = matMultRizomiliotisTriakosia(encA, encB, encC[0], statsResults, reiteraciones, lado, r, sumRowsColsKeys, fastRot);
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
    return ret == 0 ? (reiteraciones == 0 ? 2 : 0) : ret;
}