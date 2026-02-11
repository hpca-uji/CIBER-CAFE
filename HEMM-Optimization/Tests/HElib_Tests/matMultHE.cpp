#include <helib/helib.h>
#include <omp.h>
#include <map>

using namespace helib;

void matMultHE(const Ctxt& encA,
               const Ctxt& encB,
               Ctxt& encC,
               std::map<std::string, std::string>& statsResults,
               size_t& reiteraciones,
               const int lado,
               const size_t reitActual=0,
               const uint lvlStrassen=0) {
    if (lvlStrassen > 0) {
        std::cerr << "HEMatMult does not support Strassen's algorithm" << std::endl;
        addToMap(statsResults, "exceptions", "HEMatMult does not support Strassen's algorithm");
        reiteraciones = reitActual;
        return;
    }
    bool error = false;
    uint nElem = lado*lado;
    uint nSlots = encA.getContext().getNSlots();
    Ctxt encA0(encA), encB0(encB);

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
            rotate(encA0, lado-1);
        }
        encA0 *= makePlaintext(encA.getContext(), nElem, [lado](int ell){return (unsigned)(ell-2*lado+1) < 1;}, nSlots);
        #if ENABLE_DEBUG
        printf("encA0: 0/%d\n", 2*lado-1);
        #endif
        #if USE_OMP_TASKLOOP
        #pragma omp taskloop shared(encA0)
        #else
        #pragma omp parallel for shared(encA0) schedule(dynamic)
        #endif
        for (int k=-lado+2; k<lado; k++)
            try {
                Ctxt encA_temp(encA);
                if (k != 0) {
                    rotate(encA_temp, -k);
                }
                encA_temp *= makePlaintext(encA.getContext(), nElem, [lado, k](int ell){return (k>=0) ? ((unsigned)(ell-lado*k) < (lado-k)) : ((unsigned)(ell-(lado+k)*lado+k) < (lado+k));}, nSlots);
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
        encB0 *= makePlaintext(encB.getContext(), nElem, [lado](int ell){return mod(ell, lado) == 0;}, nSlots);
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
                Ctxt encB_temp(encB);
                rotate(encB_temp, -k*lado);
                encB_temp *= makePlaintext(encB.getContext(), nElem, [lado, k](int ell){return mod(ell, lado) == k;}, nSlots);
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
        encC = encA0;
        encC *= encB0;
        #if ENABLE_DEBUG
        printf("encAB: 0/%d\n", lado);
        #endif
        #if USE_OMP_TASKLOOP
        #pragma omp taskloop shared(encC)
        #else
        #pragma omp parallel for shared(encC) schedule(dynamic)
        #endif
        for (int k=1; k<lado; k++)
            try {
                Ctxt encAB_temp(encA.getPubKey()), encA0_temp1(encA0), encA0_temp2(encA0), encB0_temp(encB0);
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
                rotate(encA0_temp1, -k);
                encA0_temp1 *= makePlaintext(encA.getContext(), nElem, [lado, k](int ell){return (unsigned)(mod(ell, lado)) < lado-k;}, nSlots);
                #if USE_OMP_TASKLOOP
                }
                #pragma omp task shared(encA0_temp2)
                {
                #endif
                rotate(encA0_temp2, -k+lado);
                encA0_temp2 *= makePlaintext(encA.getContext(), nElem, [lado, k](int ell){return (unsigned)(mod(ell, lado)-lado+k) < k;}, nSlots);
                #if USE_OMP_TASKLOOP
                }
                }
                #endif
                encA0_temp1 += encA0_temp2;
                #if USE_OMP_TASKLOOP
                }
                #pragma omp task shared(encB0_temp)
                {
                #endif
                rotate(encB0_temp, -lado*k);
                #if USE_OMP_TASKLOOP
                }
                }
                #endif
                encAB_temp = encA0_temp1;
                encAB_temp *= encB0_temp;
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

template <typename T>
int MatMultHE(const Context& cc,
              const std::vector<std::vector<T>>& A,
              const std::vector<std::vector<T>>& B,
              std::vector<std::vector<std::vector<T>>>& result,
              size_t reiteraciones,
              std::map<std::string, std::string>& statsResults) {
    size_t lado = 1 << (int)ceil(log2(std::max({A.size(), A[0].size(), B.size(), B[0].size()})));
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
        for (size_t i = lado*lado; i < cc.getNSlots(); i*=2){
            flatA.insert(flatA.end(), flatA.begin(), flatA.end());
            flatB.insert(flatB.end(), flatB.begin(), flatB.end());
        }
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
                matMultHE(encC[r-1], encB, encC[r], statsResults, reiteraciones, lado, r);
            });
        } else {
            measureBlock("MatrixMultiplication", statsResults, [&](){
                matMultHE(encA, encB, encC[0], statsResults, reiteraciones, lado, r);
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
    return reiteraciones == 0 ? 2 : 0;
}