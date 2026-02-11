#include <helib/helib.h>
#include <omp.h>
#include <map>

using namespace helib;

int matMultRizomiliotisTriakosia(const Ctxt& encA,
                                 const Ctxt& encB,
                                 Ctxt& encC,
                                 std::map<std::string, std::string>& statsResults,
                                 size_t& reiteraciones,
                                 const int lado,
                                 const size_t reitActual=0,
                                 const uint lvlStrassen=0) {
    bool error = false;
    int ladoStr = lado >> lvlStrassen,
        slotsMax = encA.getContext().getNSlots();
    uint32_t nSlotsStr = std::max(ladoStr*ladoStr*ladoStr, lado*lado);
    if (nSlotsStr > slotsMax) {
        newException(statsResults, "Number of slots max. (" + std::to_string(slotsMax) + ") is less than required for Rizomiliotis-Triakosia method (" + std::to_string(nSlotsStr) + ").");
        reiteraciones = reitActual;
        return 3;
    }
    Ctxt encAT(encA),
         encBT(encB);
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
        encAT *= makePlaintext(encA.getContext(), nSlotsStr, [ladoStr, lado](int j){return j<ladoStr*lado && mod(j, ladoStr)==0 && mod(floor(j/ladoStr), lado/ladoStr)==0;});
        #if ENABLE_DEBUG
        printf("encA: 0/%d\n", ladoStr);
        #endif
        #if USE_OMP_TASKLOOP
        #pragma omp taskloop shared(encAT)
        #else
        #pragma omp parallel for shared(encAT) schedule(dynamic)
        #endif
        for (int i=1; i<ladoStr; i++) {
            try {
                Ctxt encA_temp(encA);
                int rot = (mod(i*ladoStr, lado)+floor(i*ladoStr/lado)*lado*ladoStr)-i;
                if (rot != 0) {
                    rotate(encA_temp, rot);
                }
                encA_temp *= makePlaintext(encA.getContext(), nSlotsStr, [ladoStr, lado, i](int j){return j>=ladoStr*lado*floor(i*ladoStr/lado) && j<ladoStr*lado*(floor(i*ladoStr/lado)+1) && mod(j, ladoStr)==0 && mod(floor(j/ladoStr), lado/ladoStr)==mod(i, lado/ladoStr);});
                
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
        for (int i=0; (1 << i)<ladoStr; i++) {
            Ctxt encA_temp(encAT);
            rotate(encA_temp, (1 << i));
            encAT += encA_temp;
            #if ENABLE_DEBUG
            printf("encA: Expand(%d/%d)->%d\n", i, (int)log2(ladoStr), (1 << i));
            #endif
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
        encBT *= makePlaintext(encB.getContext(), nSlotsStr, [ladoStr, lado](int j){return j<ladoStr;});
        #if ENABLE_DEBUG
        printf("encB: 0/%d\n", ladoStr);
        #endif
        #if USE_OMP_TASKLOOP
        #pragma omp taskloop shared(encBT)
        #else
        #pragma omp parallel for shared(encBT) schedule(dynamic)
        #endif
        for (int i=1; i<ladoStr; i++) {
            try {
                Ctxt encB_temp(encB);
                int rot = (mod(i*ladoStr, lado)+floor(i*ladoStr/lado)*lado*ladoStr)-i*lado;
                //int rot = -i*lado*(ladoStr-1);
                if (rot != 0) {
                    rotate(encB_temp, rot);
                }
                encB_temp *= makePlaintext(encB.getContext(), nSlotsStr, [ladoStr, lado, i](int j){return j>=(mod(i, lado/ladoStr)+lado*floor(i*ladoStr/lado))*ladoStr && j<(mod(i, lado/ladoStr)+lado*floor(i*ladoStr/lado)+1)*ladoStr;});
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
            Ctxt encB_temp(encBT);
            rotate(encB_temp, (1 << i)*lado);
            encBT += encB_temp;
            #if ENABLE_DEBUG
            printf("encB: Expand(%d/%d)->%d\n", i, (int)log2(ladoStr), (1 << i));
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
        printf("encC: Prepared(encAT(capacity=%f, errorBound=%f), encBT(capacity=%f, errorBound=%f))\n", encAT.capacity(), encAT.errorBound(), encBT.capacity(), encBT.errorBound());
        #endif
        encC = encAT;
        encC *= encBT;
        #if ENABLE_DEBUG
        printf("encC: Calculated capacity=%f, errorBound=%f\n", encC.capacity(), encC.errorBound());
        #endif
        for(int i=0; i < lvlStrassen && (1<<i) < ladoStr; i++) {
            Ctxt encC_temp(encC);
            rotate(encC_temp, -(ladoStr << i));
            encC += encC_temp;
            #if ENABLE_DEBUG
            printf("encC: Rotate.Row(%d)\n", -(ladoStr << i));
            #endif
        }
        for (int i=0; (1<<(i+lvlStrassen)) < ladoStr; i++) {
            Ctxt encC_temp(encC);
            rotate(encC_temp, -((lado*ladoStr) << i));
            encC += encC_temp;
            #if ENABLE_DEBUG
            printf("encC: Rotate.Col(%d)\n", -((lado*ladoStr) << i));
            #endif
        }
        #if ENABLE_DEBUG
        printf("encC: END capacity=%f, errorBound=%f\n", encC.capacity(), encC.errorBound());
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

template <typename T>
int MatMultRizomiliotisTriakosia(const Context& cc,
                                 const std::vector<std::vector<T>>& A,
                                 const std::vector<std::vector<T>>& B,
                                 std::vector<std::vector<std::vector<T>>>& result,
                                 size_t reiteraciones,
                                 std::map<std::string, std::string>& statsResults) {
    size_t lado = 1 << (int)ceil(log2(std::max({A.size(), A[0].size(), B.size(), B[0].size()})));
    if (lado*lado*lado > cc.getNSlots()) {
        newException(statsResults, "Number of slots max. (" + std::to_string(cc.getNSlots()) + ") is less than required for Rizomiliotis-Triakosia method (" + std::to_string(lado*lado*lado) + ").");
        return 3;
    }
    #if ENABLE_DEBUG
    std::cout << "Starting Key Generation..." << std::endl;
    #endif
    SecKey secretKey(cc);
    measureBlock("KeyGeneration", statsResults, [&](){
        secretKey.GenSecKey();
        addSome1DMatrices(secretKey);
    });
    const PubKey& publicKey = secretKey;

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
        /*
        for (size_t i = lado*lado; i < cc.getNSlots(); i*=2){
            flatA.insert(flatA.end(), flatA.begin(), flatA.end());
            flatB.insert(flatB.end(), flatB.begin(), flatB.end());
        }
        */
    });
    
    PtxtArray ptxA(cc), ptxB(cc);
    measureBlock("Encode", statsResults, [&](){
        ptxA = PtxtArray(cc, flatA);
        ptxB = PtxtArray(cc, flatB);
    });

    Ctxt encA(publicKey), encB(publicKey);
    measureBlock("Encrypt", statsResults, [&](){
        ptxA.encrypt(encA);
        ptxB.encrypt(encB);
    });
    addToMap(statsResults, "ctxA_size", to_string_precise(1));
    calculateStats("ctxA_capacity", statsResults, std::vector<double>{static_cast<double>(encA.capacity())});
    calculateStats("ctxA_errorBound", statsResults, std::vector<double>{static_cast<double>(encA.errorBound())});
    addToMap(statsResults, "ctxB_size", to_string_precise(1));
    calculateStats("ctxB_capacity", statsResults, std::vector<double>{static_cast<double>(encB.capacity())});
    calculateStats("ctxB_errorBound", statsResults, std::vector<double>{static_cast<double>(encB.errorBound())});
    
    std::vector<Ctxt> encC(reiteraciones, Ctxt(publicKey));
    std::vector<PtxtArray> ptxC(reiteraciones, PtxtArray(cc));
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
                ret = matMultRizomiliotisTriakosia(encC[r-1], encB, encC[r], statsResults, reiteraciones, lado, r);
            });
        } else {
            measureBlock("MatrixMultiplication", statsResults, [&](){
                ret = matMultRizomiliotisTriakosia(encA, encB, encC[0], statsResults, reiteraciones, lado, r);
            });
        }
        if(reiteraciones > r) {
            addToMap(statsResults, "ctxC_size", to_string_precise(1));
            calculateStats("ctxC_capacity", statsResults, std::vector<double>{static_cast<double>(encC[r].capacity())});
            calculateStats("ctxC_errorBound", statsResults, std::vector<double>{static_cast<double>(encC[r].errorBound())});
            measureBlock("Decrypt", statsResults, [&](){
                try {
                    ptxC[r].decrypt(encC[r], secretKey);
                } catch (const std::exception& e) {
                    newException(e, statsResults, "Decrypt_" + to_string_precise(r));
                    reiteraciones = r; // Stop further repetitions
                }
            });
            if(reiteraciones > r) {
                measureBlock("Decode", statsResults, [&](){
                    try {
                        ptxC[r].store(decC[r]);
                    } catch (const std::exception& e) {
                        newException(e, statsResults, "Decode_" + to_string_precise(r));
                        reiteraciones = r; // Stop further repetitions
                    }
                });
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