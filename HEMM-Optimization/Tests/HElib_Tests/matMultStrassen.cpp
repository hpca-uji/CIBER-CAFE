#include <helib/helib.h>
#include <omp.h>
#include <map>

using namespace helib;

int rotQuadrantStrassen(const uint quadrantOrig, const uint quadrantDest, const uint lado, const uint lvlStrassen=0) {
    if (quadrantOrig == quadrantDest) {
        return 0;
    } else {
        uint quadMin = std::min(quadrantOrig, quadrantDest),
             quadMax = std::max(quadrantOrig, quadrantDest);
        int ret=(lado >> (lvlStrassen+1));
        switch (quadMin^quadMax) {
            case 2:
                ret*=lado;
                break;
            case 3:
                ret*=lado+(quadMin==1 ? -1 : 1);
                break;
        }
        if (quadrantDest<quadrantOrig)
            ret = -ret;
        return ret;
    }
}

void matMultStrassen(const Ctxt& encA,
                     const Ctxt& encB,
                     Ctxt& encC,
                     std::map<std::string, std::string>& statsResults,
                     size_t& reiteraciones,
                     const int lado,
                     const uint strassenAtSize=1,
                     const Algorithm subAlgorithm=Algorithm::Traditional,
                     const size_t reitActual=0,
                     const uint lvlStrassen=0) {
    #if ENABLE_DEBUG
    std::cout << "matMultStrassen(lado=" << lado << ", strassenAtSize=" << strassenAtSize << ", lvlStrassen=" << lvlStrassen << ", subAlgorithm=" << AlgorithmToString.at(subAlgorithm) << ", reiteraciones=" << reiteraciones << ", reitActual=" << reitActual << ")" << std::endl;
    #endif
    int ladoStrassen = lado/(1 << lvlStrassen);
    if (ladoStrassen <= strassenAtSize || strassenAtSize == 0) {
        if (ladoStrassen == 1) {
            #if ENABLE_DEBUG
            std::cout << "matMultStrassen(lado=" << lado << ", strassenAtSize=" << strassenAtSize << ", lvlStrassen=" << lvlStrassen << "): Base case 1x1x1 multiplication." << std::endl;
            #endif
            encC = encA;
            encC *= encB;
        } else {
            #if ENABLE_DEBUG
            std::cout << "matMultStrassen(lado=" << lado << ", strassenAtSize=" << strassenAtSize << ", lvlStrassen=" << lvlStrassen << "): Base case multiplication with subAlgorithm=" << AlgorithmToString.at(subAlgorithm) << "." << std::endl;
            #endif
            switch (subAlgorithm) {
                case HEMatMult:
                    matMultHE(encA, encB, encC, statsResults, reiteraciones, lado, reitActual, lvlStrassen);
                    break;
                case RizomiliotisTriakosia:
                    matMultRizomiliotisTriakosia(encA, encB, encC, statsResults, reiteraciones, lado, reitActual, lvlStrassen);
                    break;
                default:
                    newException(statsResults, "matMultStrassen: Algorithm not implemented for base case multiplication: " + AlgorithmToString.at(subAlgorithm) + ".");
                    reiteraciones = reitActual;
                    return;
            }
        }
        if (lvlStrassen > 0) {
            #if ENABLE_DEBUG
            std::cout << "encC: Clean " << ladoStrassen << " firsts columns and rows of " << lado << " columns, slots: " << encC.getContext().getNSlots() << std::endl;
            //std::cout << "\t" << vec_from_pred<double>(encC.getContext().getNSlots(), [lado, ladoStrassen](int j){return j < lado*ladoStrassen && mod(j, lado) < ladoStrassen;}) << std::endl;
            #endif
            encC *= makePlaintext(encC.getContext(), encC.getContext().getNSlots(), [lado, ladoStrassen](int j){return j < lado*ladoStrassen && mod(j, lado) < ladoStrassen;});
        }
    } else {
        std::vector<Ctxt> encA_rot = std::vector<Ctxt>(3, encA),
                          encB_rot = std::vector<Ctxt>(3, encB),
                          enc_prep_m = std::vector<Ctxt>(10, Ctxt(encA.getPubKey())),
                          enc_m = std::vector<Ctxt>(7, Ctxt(encA.getPubKey())),
                          enc_c = std::vector<Ctxt>(4, Ctxt(encA.getPubKey()));
        // Cuadrantes
        #if USE_OMP_TASKLOOP
        bool error = false;
        #pragma omp taskgroup
        {
        try {
        #else
        #pragma omp parallel for shared(encA_rot, encB_rot) schedule(dynamic)
        #endif
        for (int i=0; i<3; i++) {
            int rot = rotQuadrantStrassen(i+1, 0, lado, lvlStrassen);
            #if USE_OMP_TASKLOOP
            #pragma omp task shared(encA_rot)
            {
            #endif
            #if ENABLE_DEBUG
            printf("lvl%d -> RotA(%d->0): %d\n", lvlStrassen, i+1, rot);
            #endif
            rotate(encA_rot[i], rot);
            #if USE_OMP_TASKLOOP
            }
            #pragma omp task shared(encB_rot)
            {
            #endif
            #if ENABLE_DEBUG
            printf("lvl%d -> RotB(%d->0): %d\n", lvlStrassen, i+1, rot);
            #endif
            rotate(encB_rot[i], rot);
            #if USE_OMP_TASKLOOP
            }
            #endif
        }
        #if USE_OMP_TASKLOOP
        } catch (const std::exception& e) {
            newException(e, statsResults, "MatrixMultiplication_" + to_string_precise(reitActual) + "_"+to_string_precise(lvlStrassen)+".Rot(encA & encB)");
            reiteraciones = reitActual;
            error = true;
        }
        }
        if (!error) {
        // Prepara M1-M7
        #pragma omp taskgroup
        {
        #pragma omp task shared(enc_prep_m)
        {
        #endif
        #if ENABLE_DEBUG
        printf("lvl%d -> PM0=M1.1=A1+A4\n", lvlStrassen);
        #endif
        enc_prep_m[0] = encA;
        enc_prep_m[0] += encA_rot[2];  // PM0
        #if USE_OMP_TASKLOOP
        }
        #pragma omp task shared(enc_prep_m)
        {
        #endif
        #if ENABLE_DEBUG
        printf("lvl%d -> PM1=M1.2=B1+B4\n", lvlStrassen);
        #endif
        enc_prep_m[1] = encB;
        enc_prep_m[1] += encB_rot[2];  // PM1
        #if USE_OMP_TASKLOOP
        }
        #pragma omp task shared(enc_prep_m)
        {
        #endif
        #if ENABLE_DEBUG
        printf("lvl%d -> PM2=M2.1=A3+A4\n", lvlStrassen);
        #endif
        enc_prep_m[2] = encA_rot[1];
        enc_prep_m[2] += encA_rot[2]; // PM2
        #if USE_OMP_TASKLOOP
        }
        #pragma omp task shared(enc_prep_m)
        {
        #endif
        #if ENABLE_DEBUG
        printf("lvl%d -> PM3=M3.2=B2-B4\n", lvlStrassen);
        #endif
        enc_prep_m[3] = encB_rot[0];
        enc_prep_m[3] -= encB_rot[2]; // PM3
        #if USE_OMP_TASKLOOP
        }
        #pragma omp task shared(enc_prep_m)
        {
        #endif
        #if ENABLE_DEBUG
        printf("lvl%d -> PM4=M4.2=B3-B1\n", lvlStrassen);
        #endif
        enc_prep_m[4] = encB_rot[1];
        enc_prep_m[4] -= encB;  // PM4
        #if USE_OMP_TASKLOOP
        }
        #pragma omp task shared(enc_prep_m)
        {
        #endif
        #if ENABLE_DEBUG
        printf("lvl%d -> PM5=M5.1=A1+A2\n", lvlStrassen);
        #endif
        enc_prep_m[5] = encA;
        enc_prep_m[5] += encA_rot[0];  // PM5
        #if USE_OMP_TASKLOOP
        }
        #pragma omp task shared(enc_prep_m)
        {
        #endif
        #if ENABLE_DEBUG
        printf("lvl%d -> PM6=M6.1=A3-A1\n", lvlStrassen);
        #endif
        enc_prep_m[6] = encA_rot[1];
        enc_prep_m[6] -= encA;  // PM6
        #if USE_OMP_TASKLOOP
        }
        #pragma omp task shared(enc_prep_m)
        {
        #endif
        #if ENABLE_DEBUG
        printf("lvl%d -> PM7=M6.2=B1+B2\n", lvlStrassen);
        #endif
        enc_prep_m[7] = encB;
        enc_prep_m[7] += encB_rot[0];  // PM7
        #if USE_OMP_TASKLOOP
        }
        #pragma omp task shared(enc_prep_m)
        {
        #endif
        #if ENABLE_DEBUG
        printf("lvl%d -> PM8=M7.1=A2-A4\n", lvlStrassen);
        #endif
        enc_prep_m[8] = encA_rot[0];
        enc_prep_m[8] -= encA_rot[2]; // PM8
        #if USE_OMP_TASKLOOP
        }
        #pragma omp task shared(enc_prep_m)
        {
        #endif
        #if ENABLE_DEBUG
        printf("lvl%d -> PM9=M7.2=B3+B4\n", lvlStrassen);
        #endif
        enc_prep_m[9] = encB_rot[1];
        enc_prep_m[9] += encB_rot[2];  // PM9
        #if USE_OMP_TASKLOOP
        }
        }
        // M1-M7
        #pragma omp taskgroup
        {
        #pragma omp task shared(enc_m)
        try {
        #endif
        #if ENABLE_DEBUG
        printf("lvl%d -> M1=(PM0=M1.1=A1+A4)*(PM1=M1.2=B1+B4)\n", lvlStrassen);
        #endif
        matMultStrassen(enc_prep_m[0], enc_prep_m[1], enc_m[0], statsResults, reiteraciones, lado, strassenAtSize, subAlgorithm, reitActual, lvlStrassen+1); // M1
        #if USE_OMP_TASKLOOP
        } catch (const std::exception& e) {
            newException(e, statsResults, "MatrixMultiplication_" + to_string_precise(reitActual) + "_"+to_string_precise(lvlStrassen)+".M1");
            reiteraciones = reitActual;
            error = true;
        }
        #pragma omp task shared(enc_m)
        try {
        #endif
        #if ENABLE_DEBUG
        printf("lvl%d -> M2=(PM2=M2.1=A3+A4)*B1\n", lvlStrassen);
        #endif
        matMultStrassen(enc_prep_m[2], encB, enc_m[1], statsResults, reiteraciones, lado, strassenAtSize, subAlgorithm, reitActual, lvlStrassen+1); // M2
        #if USE_OMP_TASKLOOP
        } catch (const std::exception& e) {
            newException(e, statsResults, "MatrixMultiplication_" + to_string_precise(reitActual) + "_"+to_string_precise(lvlStrassen)+".M2");
            reiteraciones = reitActual;
            error = true;
        }
        #pragma omp task shared(enc_m)
        try {
        #endif
        #if ENABLE_DEBUG
        printf("lvl%d -> M3=A1*(PM3=M3.2=B2-B4)\n", lvlStrassen);
        #endif
        matMultStrassen(encA, enc_prep_m[3], enc_m[2], statsResults, reiteraciones, lado, strassenAtSize, subAlgorithm, reitActual, lvlStrassen+1); // M3
        #if USE_OMP_TASKLOOP
        } catch (const std::exception& e) {
            newException(e, statsResults, "MatrixMultiplication_" + to_string_precise(reitActual) + "_"+to_string_precise(lvlStrassen)+".M3");
            reiteraciones = reitActual;
            error = true;
        }
        #pragma omp task shared(enc_m)
        try {
        #endif
        #if ENABLE_DEBUG
        printf("lvl%d -> M4=A4*(PM4=M4.2=B3-B1)\n", lvlStrassen);
        #endif
        matMultStrassen(encA_rot[2], enc_prep_m[4], enc_m[3], statsResults, reiteraciones, lado, strassenAtSize, subAlgorithm, reitActual, lvlStrassen+1); // M4
        #if USE_OMP_TASKLOOP
        } catch (const std::exception& e) {
            newException(e, statsResults, "MatrixMultiplication_" + to_string_precise(reitActual) + "_"+to_string_precise(lvlStrassen)+".M4");
            reiteraciones = reitActual;
            error = true;
        }
        #pragma omp task shared(enc_m)
        try {
        #endif
        #if ENABLE_DEBUG
        printf("lvl%d -> M5=(PM5=M5.1=A1+A2)*B4\n", lvlStrassen);
        #endif
        matMultStrassen(enc_prep_m[5], encB_rot[2], enc_m[4], statsResults, reiteraciones, lado, strassenAtSize, subAlgorithm, reitActual, lvlStrassen+1); // M5
        #if USE_OMP_TASKLOOP
        } catch (const std::exception& e) {
            newException(e, statsResults, "MatrixMultiplication_" + to_string_precise(reitActual) + "_"+to_string_precise(lvlStrassen)+".M5");
            reiteraciones = reitActual;
            error = true;
        }
        #pragma omp task shared(enc_m)
        try {
        #endif
        #if ENABLE_DEBUG
        printf("lvl%d -> M6=(PM6=M6.1=A3-A1)*(PM7=M6.2=B1+B2)\n", lvlStrassen);
        #endif
        matMultStrassen(enc_prep_m[6], enc_prep_m[7], enc_m[5], statsResults, reiteraciones, lado, strassenAtSize, subAlgorithm, reitActual, lvlStrassen+1); // M6
        #if USE_OMP_TASKLOOP
        } catch (const std::exception& e) {
            newException(e, statsResults, "MatrixMultiplication_" + to_string_precise(reitActual) + "_"+to_string_precise(lvlStrassen)+".M6");
            reiteraciones = reitActual;
            error = true;
        }
        #pragma omp task shared(enc_m)
        try {
        #endif
        #if ENABLE_DEBUG
        printf("lvl%d -> M7=(PM8=M7.1=A2-A4)*(PM9=M7.2=B3+B4)\n", lvlStrassen);
        #endif
        matMultStrassen(enc_prep_m[8], enc_prep_m[9], enc_m[6], statsResults, reiteraciones, lado, strassenAtSize, subAlgorithm, reitActual, lvlStrassen+1); // M7
        #if USE_OMP_TASKLOOP
        } catch (const std::exception& e) {
            newException(e, statsResults, "MatrixMultiplication_" + to_string_precise(reitActual) + "_"+to_string_precise(lvlStrassen)+".M7");
            reiteraciones = reitActual;
            error = true;
        }
        }
        #endif
        encA_rot.clear();
        encB_rot.clear();
        enc_prep_m.clear();
        // C1-C4
        #if USE_OMP_TASKLOOP
        if (!error) {
        #pragma omp taskgroup
        {
        #pragma omp task shared(enc_c)
        {
        #endif
        #if ENABLE_DEBUG
        printf("lvl%d -> C1=M1+M4+M7-M5\n", lvlStrassen);
        #endif
        enc_c[0] = enc_m[0];
        enc_c[0] += enc_m[3];
        enc_c[0] += enc_m[6];
        enc_c[0] -= enc_m[4];  // C1
        #if USE_OMP_TASKLOOP
        }
        #pragma omp task shared(enc_c)
        {
        #endif
        #if ENABLE_DEBUG
        printf("lvl%d -> C2=Rot(M3+M5, %d)\n", lvlStrassen, rotQuadrantStrassen(0, 1, lado, lvlStrassen));
        #endif
        enc_c[1] = enc_m[2];
        enc_c[1] += enc_m[4];
        rotate(enc_c[1], rotQuadrantStrassen(0, 1, lado, lvlStrassen));   // C2
        #if USE_OMP_TASKLOOP
        }
        #pragma omp task shared(enc_c)
        {
        #endif
        #if ENABLE_DEBUG
        printf("lvl%d -> C3=Rot(M2+M4, %d)\n", lvlStrassen, rotQuadrantStrassen(0, 2, lado, lvlStrassen));
        #endif
        enc_c[2] = enc_m[1];
        enc_c[2] += enc_m[3];
        rotate(enc_c[2], rotQuadrantStrassen(0, 2, lado, lvlStrassen));   // C3
        #if USE_OMP_TASKLOOP
        }
        #pragma omp task shared(enc_c)
        {
        #endif
        #if ENABLE_DEBUG
        printf("lvl%d -> C4=Rot(M1+M3+M6-M2, %d)\n", lvlStrassen, rotQuadrantStrassen(0, 3, lado, lvlStrassen));
        #endif
        enc_c[3] = enc_m[0];
        enc_c[3] += enc_m[2];
        enc_c[3] += enc_m[5];
        enc_c[3] -= enc_m[1];
        rotate(enc_c[3], rotQuadrantStrassen(0, 3, lado, lvlStrassen));    // C4
        #if USE_OMP_TASKLOOP
        }
        }
        }
        }
        #endif
        enc_m.clear();
        encC = enc_c[0];
        encC += enc_c[1];
        encC += enc_c[2];
        encC += enc_c[3];
        enc_c.clear();
        if (lvlStrassen != 0)
            try {
                encC *= makePlaintext(encC.getContext(), encC.getContext().getNSlots(), [lado, ladoStrassen](int i){return i < lado * ladoStrassen && mod(i, lado) < ladoStrassen;});
            } catch (const std::exception& e) {
                newException(e, statsResults, "MatrixMultiplication_" + to_string_precise(reitActual) + "_"+to_string_precise(lvlStrassen)+".Clean");
                reiteraciones = reitActual;
                return;
            }
    }
}

template <typename T>
int MatMultStrassen(const Context& cc,
                    const std::vector<std::vector<T>>& A,
                    const std::vector<std::vector<T>>& B,
                    std::vector<std::vector<std::vector<T>>>& result,
                    size_t reiteraciones,
                    std::map<std::string, std::string>& statsResults,
                    size_t strassenAtSize=0,
                    const Algorithm subAlgorithm=Algorithm::Traditional) {
    size_t lado = 1 << (int)ceil(log2(std::max({A.size(), A[0].size(), B.size(), B[0].size()})));
    if (strassenAtSize == 0) {
        switch (subAlgorithm)
        {
        case Algorithm::RizomiliotisTriakosia:
            if (cc.getNSlots() < lado*lado*lado)
                strassenAtSize = (1 << (int)(floor(log2(cbrt(cc.getNSlots())))));
            break;
        default:
            strassenAtSize = 1;
            break;
        }
    }

    #if ENABLE_DEBUG
    std::cout << "strassenAtSize=" << strassenAtSize << std::endl;
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
                matMultStrassen(encC[r-1], encB, encC[r], statsResults, reiteraciones, lado, strassenAtSize, subAlgorithm, r);
            });
        } else {
            measureBlock("MatrixMultiplication", statsResults, [&](){
                matMultStrassen(encA, encB, encC[0], statsResults, reiteraciones, lado, strassenAtSize, subAlgorithm, r);
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
                if(reiteraciones > r) {
                    measureBlock("Postprocess", statsResults, [&](){
                        #pragma omp parallel for shared(result) schedule(dynamic)
                        for (size_t i = 0; i < A.size(); ++i)
                            result[r][i] = std::vector<T>(decC[r].begin() + i * lado, decC[r].begin() + i * lado + B[0].size());
                    });
                }
            }
        }
    }
    #if USE_OMP_TASKLOOP
    }
    }
    #endif
    return reiteraciones == 0 ? 2 : 0;
}