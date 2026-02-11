#include <openfhe.h>
#include <omp.h>
#include <map>

using namespace lbcrypto;

int rotQuadrantStrassen(const usint quadrantOrig, const usint quadrantDest, const usint lado, const usint lvlStrassen=0) {
    if (quadrantOrig == quadrantDest) {
        return 0;
    } else {
        usint quadMin = std::min(quadrantOrig, quadrantDest),
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
        if (quadrantDest>quadrantOrig)
            ret = -ret;
        return ret;
    }
}

std::set<int> rotStrassen(std::map<std::string, std::string>& statsResults,
                          const usint lado,
                          const usint strassenAtSize=1,
                          const Algorithm subAlg=Algorithm::Traditional,
                          const usint lvlStrassen=0) {
    std::set<int> ret;
    if (strassenAtSize == 0 || (lado >> lvlStrassen) <= strassenAtSize) {
        switch (subAlg) {
            case Algorithm::RizomiliotisTriakosia:
            case Algorithm::RizomiliotisTriakosiaFastRot:
                ret = rotRizomiliotisTriakosia(lado, false, lvlStrassen);
                break;
            case Algorithm::RizomiliotisTriakosiaSums:
            case Algorithm::RizomiliotisTriakosiaSumsFastRot:
                ret = rotRizomiliotisTriakosia(lado, true, lvlStrassen);
                break;
            case Algorithm::Traditional:
                break;
            default:
                newException(statsResults, "rotStrassen("+std::to_string(lado)+", "+std::to_string(strassenAtSize)+"): Algorithm invalid (" + AlgorithmToString.at(subAlg).c_str() + ") for Strassen rotations.");
                break;
        }
    } else {
        ret = {
            rotQuadrantStrassen(1, 0, lado, lvlStrassen),
            rotQuadrantStrassen(2, 0, lado, lvlStrassen),
            rotQuadrantStrassen(3, 0, lado, lvlStrassen),
            rotQuadrantStrassen(0, 1, lado, lvlStrassen),
            rotQuadrantStrassen(0, 2, lado, lvlStrassen),
            rotQuadrantStrassen(0, 3, lado, lvlStrassen),
        };
        for (int i : rotStrassen(statsResults, lado, strassenAtSize, subAlg, lvlStrassen+1))
            ret.insert(i);
    }
    #if ENABLE_DEBUG
    printf("rotStrassen(%u, %u, %s, %u)(%lu): ", (unsigned)lado, (unsigned)strassenAtSize, AlgorithmToString.at(subAlg).c_str(), (unsigned)lvlStrassen, ret.size());
    for (const auto& r : ret)
        printf("%d ", r);
    printf("\n");
    #endif
    return ret;
}

template <class Element>
void matMultStrassen(const Ciphertext<Element>& encA,
                     const Ciphertext<Element>& encB,
                     Ciphertext<Element>& encC,
                     std::map<std::string, std::string>& statsResults,
                     size_t& reiteraciones,
                     const int lado,
                     const uint strassenAtSize=1,
                     const Algorithm subAlgorithm=Algorithm::Traditional,
                     const size_t reitActual=0,
                     const std::shared_ptr<std::map<uint32_t, EvalKey<Element>>> sumRowsColsKeys = nullptr,
                     const bool fastRot=false,
                     const usint lvlStrassen=0) {
    #if ENABLE_DEBUG
    std::cout << "matMultStrassen(lado=" << lado << ", strassenAtSize=" << strassenAtSize << ", lvlStrassen=" << lvlStrassen << ", subAlgorithm=" << AlgorithmToString.at(subAlgorithm) << ", reiteraciones=" << reiteraciones << ", reitActual=" << reitActual << ", sumRowsCols=" << (sumRowsColsKeys ? "true" : "false") << ", fastRot=" << (fastRot ? "true" : "false") << ")" << std::endl;
    #endif
    CryptoContext<Element> cc = encA->GetCryptoContext();
    int ladoStrassen = lado/(1 << lvlStrassen);
    bool isCKKS = (encA->GetEncodingType() == CKKS_PACKED_ENCODING);
    if (ladoStrassen <= strassenAtSize || strassenAtSize == 0) {
        if (ladoStrassen == 1) {
            #if ENABLE_DEBUG
            std::cout << "matMultStrassen(lado=" << lado << ", strassenAtSize=" << strassenAtSize << ", lvlStrassen=" << lvlStrassen << "): Base case 1x1x1 multiplication." << std::endl;
            #endif
            encC = cc->EvalMult(encA, encB);
        } else {
            #if ENABLE_DEBUG
            std::cout << "matMultStrassen(lado=" << lado << ", strassenAtSize=" << strassenAtSize << ", lvlStrassen=" << lvlStrassen << "): Base case multiplication with subAlgorithm=" << AlgorithmToString.at(subAlgorithm) << "." << std::endl;
            #endif
            switch (subAlgorithm) {
                case HEMatMult:
                case HEMatMultFastRot:
                    matMultHE(encA, encB, encC, statsResults, reiteraciones, lado, reitActual, fastRot, lvlStrassen);
                    break;
                case RizomiliotisTriakosia:
                case RizomiliotisTriakosiaFastRot:
                case RizomiliotisTriakosiaSums:
                case RizomiliotisTriakosiaSumsFastRot:
                    matMultRizomiliotisTriakosia(encA, encB, encC, statsResults, reiteraciones, lado, reitActual, sumRowsColsKeys, fastRot, lvlStrassen);
                    break;
                default:
                    newException(statsResults, "matMultStrassen: Algorithm not implemented for base case multiplication: " + AlgorithmToString.at(subAlgorithm) + ".");
                    reiteraciones = reitActual;
                    return;
            }

        }
        if (lvlStrassen > 0) {
            #if ENABLE_DEBUG
            std::cout << "encC: Clean " << ladoStrassen << " firsts columns and rows of " << lado << " columns and " << encC->GetSlots() << " slots\n";
            //std::cout << "\t" << vec_from_pred<double>(encC->GetSlots(), [lado, ladoStrassen](int j){return j < lado*ladoStrassen && mod(j, lado) < ladoStrassen;}) << std::endl;
            #endif
            encC = cc->EvalMult(encC,
                makePlaintext(cc, isCKKS, encC->GetSlots(), [lado, ladoStrassen](int j){return j < lado*ladoStrassen && mod(j, lado) < ladoStrassen;})
            );
        }
    } else {
        std::vector<Ciphertext<Element>> encA_rot = std::vector<Ciphertext<Element>>(3),
                                         encB_rot = std::vector<Ciphertext<Element>>(3),
                                         enc_prep_m = std::vector<Ciphertext<Element>>(10),
                                         enc_m = std::vector<Ciphertext<Element>>(7),
                                         enc_c = std::vector<Ciphertext<Element>>(4);
        // Cuadrantes
        std::shared_ptr<std::vector<Element>> precompEncA, precompEncB;
        if (fastRot) {
            #if USE_OMP_TASKLOOP
            #pragma omp taskgroup
            {
            #pragma omp task shared(precompEncA)
            #endif
            precompEncA = cc->EvalFastRotationPrecompute(encA);
            #if USE_OMP_TASKLOOP
            #pragma omp task shared(precompEncB)
            #endif
            precompEncB = cc->EvalFastRotationPrecompute(encB);
            #if USE_OMP_TASKLOOP
            }
            #endif
        }
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
            if (fastRot) {
                encA_rot[i] = cc->EvalFastRotation(encA, rot, cc->GetCyclotomicOrder(), precompEncA);
            } else {
                encA_rot[i] = cc->EvalRotate(encA, rot);
            }
            #if USE_OMP_TASKLOOP
            }
            #pragma omp task shared(encB_rot)
            {
            #endif
            #if ENABLE_DEBUG
            printf("lvl%d -> RotB(%d->0): %d\n", lvlStrassen, i+1, rot);
            #endif
            if (fastRot) {
                encB_rot[i] = cc->EvalFastRotation(encB, rot, cc->GetCyclotomicOrder(), precompEncB);
            } else {
                encB_rot[i] = cc->EvalRotate(encB, rot);
            }
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
        enc_prep_m[0] = cc->EvalAdd(encA, encA_rot[2]);  // PM0
        #if USE_OMP_TASKLOOP
        }
        #pragma omp task shared(enc_prep_m)
        {
        #endif
        #if ENABLE_DEBUG
        printf("lvl%d -> PM1=M1.2=B1+B4\n", lvlStrassen);
        #endif
        enc_prep_m[1] = cc->EvalAdd(encB, encB_rot[2]);  // PM1
        #if USE_OMP_TASKLOOP
        }
        #pragma omp task shared(enc_prep_m)
        {
        #endif
        #if ENABLE_DEBUG
        printf("lvl%d -> PM2=M2.1=A3+A4\n", lvlStrassen);
        #endif
        enc_prep_m[2] = cc->EvalAdd(encA_rot[1], encA_rot[2]); // PM2
        #if USE_OMP_TASKLOOP
        }
        #pragma omp task shared(enc_prep_m)
        {
        #endif
        #if ENABLE_DEBUG
        printf("lvl%d -> PM3=M3.2=B2-B4\n", lvlStrassen);
        #endif
        enc_prep_m[3] = cc->EvalSub(encB_rot[0], encB_rot[2]); // PM3
        #if USE_OMP_TASKLOOP
        }
        #pragma omp task shared(enc_prep_m)
        {
        #endif
        #if ENABLE_DEBUG
        printf("lvl%d -> PM4=M4.2=B3-B1\n", lvlStrassen);
        #endif
        enc_prep_m[4] = cc->EvalSub(encB_rot[1], encB);  // PM4
        #if USE_OMP_TASKLOOP
        }
        #pragma omp task shared(enc_prep_m)
        {
        #endif
        #if ENABLE_DEBUG
        printf("lvl%d -> PM5=M5.1=A1+A2\n", lvlStrassen);
        #endif
        enc_prep_m[5] = cc->EvalAdd(encA, encA_rot[0]);  // PM5
        #if USE_OMP_TASKLOOP
        }
        #pragma omp task shared(enc_prep_m)
        {
        #endif
        #if ENABLE_DEBUG
        printf("lvl%d -> PM6=M6.1=A3-A1\n", lvlStrassen);
        #endif
        enc_prep_m[6] = cc->EvalSub(encA_rot[1], encA);  // PM6
        #if USE_OMP_TASKLOOP
        }
        #pragma omp task shared(enc_prep_m)
        {
        #endif
        #if ENABLE_DEBUG
        printf("lvl%d -> PM7=M6.2=B1+B2\n", lvlStrassen);
        #endif
        enc_prep_m[7] = cc->EvalAdd(encB, encB_rot[0]);  // PM7
        #if USE_OMP_TASKLOOP
        }
        #pragma omp task shared(enc_prep_m)
        {
        #endif
        #if ENABLE_DEBUG
        printf("lvl%d -> PM8=M7.1=A2-A4\n", lvlStrassen);
        #endif
        enc_prep_m[8] = cc->EvalSub(encA_rot[0], encA_rot[2]); // PM8
        #if USE_OMP_TASKLOOP
        }
        #pragma omp task shared(enc_prep_m)
        {
        #endif
        #if ENABLE_DEBUG
        printf("lvl%d -> PM9=M7.2=B3+B4\n", lvlStrassen);
        #endif
        enc_prep_m[9] = cc->EvalAdd(encB_rot[1], encB_rot[2]);  // PM9
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
        matMultStrassen(enc_prep_m[0], enc_prep_m[1], enc_m[0], statsResults, reiteraciones, lado, strassenAtSize, subAlgorithm, reitActual, sumRowsColsKeys, fastRot, lvlStrassen+1); // M1
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
        matMultStrassen(enc_prep_m[2], encB, enc_m[1], statsResults, reiteraciones, lado, strassenAtSize, subAlgorithm, reitActual, sumRowsColsKeys, fastRot, lvlStrassen+1); // M2
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
        matMultStrassen(encA, enc_prep_m[3], enc_m[2], statsResults, reiteraciones, lado, strassenAtSize, subAlgorithm, reitActual, sumRowsColsKeys, fastRot, lvlStrassen+1); // M3
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
        matMultStrassen(encA_rot[2], enc_prep_m[4], enc_m[3], statsResults, reiteraciones, lado, strassenAtSize, subAlgorithm, reitActual, sumRowsColsKeys, fastRot, lvlStrassen+1); // M4
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
        matMultStrassen(enc_prep_m[5], encB_rot[2], enc_m[4], statsResults, reiteraciones, lado, strassenAtSize, subAlgorithm, reitActual, sumRowsColsKeys, fastRot, lvlStrassen+1); // M5
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
        matMultStrassen(enc_prep_m[6], enc_prep_m[7], enc_m[5], statsResults, reiteraciones, lado, strassenAtSize, subAlgorithm, reitActual, sumRowsColsKeys, fastRot, lvlStrassen+1); // M6
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
        matMultStrassen(enc_prep_m[8], enc_prep_m[9], enc_m[6], statsResults, reiteraciones, lado, strassenAtSize, subAlgorithm, reitActual, sumRowsColsKeys, fastRot, lvlStrassen+1); // M7
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
        enc_c[0] = cc->EvalSub(cc->EvalAddMany({enc_m[0], enc_m[3], enc_m[6]}), enc_m[4]);  // C1
        #if USE_OMP_TASKLOOP
        }
        #pragma omp task shared(enc_c)
        {
        #endif
        #if ENABLE_DEBUG
        printf("lvl%d -> C2=Rot(M3+M5, %d)\n", lvlStrassen, rotQuadrantStrassen(0, 1, lado, lvlStrassen));
        #endif
        enc_c[1] = cc->EvalRotate(cc->EvalAdd(enc_m[2], enc_m[4]), rotQuadrantStrassen(0, 1, lado, lvlStrassen));   // C2
        #if USE_OMP_TASKLOOP
        }
        #pragma omp task shared(enc_c)
        {
        #endif
        #if ENABLE_DEBUG
        printf("lvl%d -> C3=Rot(M2+M4, %d)\n", lvlStrassen, rotQuadrantStrassen(0, 2, lado, lvlStrassen));
        #endif
        enc_c[2] = cc->EvalRotate(cc->EvalAdd(enc_m[1], enc_m[3]), rotQuadrantStrassen(0, 2, lado, lvlStrassen));   // C3
        #if USE_OMP_TASKLOOP
        }
        #pragma omp task shared(enc_c)
        {
        #endif
        #if ENABLE_DEBUG
        printf("lvl%d -> C4=Rot(M1+M3+M6-M2, %d)\n", lvlStrassen, rotQuadrantStrassen(0, 3, lado, lvlStrassen));
        #endif
        enc_c[3] = cc->EvalRotate(cc->EvalSub(cc->EvalAddMany({enc_m[0], enc_m[2], enc_m[5]}), enc_m[1]), rotQuadrantStrassen(0, 3, lado, lvlStrassen));    // C4
        #if USE_OMP_TASKLOOP
        }
        }
        }
        }
        #endif
        enc_m.clear();
        encC = cc->EvalAddMany(enc_c);
        enc_c.clear();
        if (lvlStrassen != 0)
            try {
                encC = cc->EvalMult(encC, makePlaintext(cc, isCKKS, encC->GetSlots(), [lado, ladoStrassen](int i){return i < lado * ladoStrassen && mod(i, lado) < ladoStrassen;}));
            } catch (const std::exception& e) {
                newException(e, statsResults, "MatrixMultiplication_" + to_string_precise(reitActual) + "_"+to_string_precise(lvlStrassen)+".Clean");
                reiteraciones = reitActual;
                return;
            }
    }
}

template <class Element, typename T>
int MatMultStrassen(const CryptoContext<Element>& cc,
                    const std::vector<std::vector<T>>& A,
                    const std::vector<std::vector<T>>& B,
                    std::vector<std::vector<std::vector<T>>>& result,
                    size_t reiteraciones,
                    std::map<std::string, std::string>& statsResults,
                    size_t strassenAtSize=0,
                    const Algorithm subAlgorithm=Algorithm::Traditional,
                    const bool fastRot=false,
                    const bool sumRowsCols=false) {
    size_t lado = 1 << (int)ceil(log2(std::max({A.size(), A[0].size(), B.size(), B[0].size()})));
    if (strassenAtSize == 0) {
        switch (subAlgorithm)
        {
        case Algorithm::RizomiliotisTriakosia:
        case Algorithm::RizomiliotisTriakosiaFastRot:
        case Algorithm::RizomiliotisTriakosiaSums:
        case Algorithm::RizomiliotisTriakosiaSumsFastRot:
            if (cc->GetRingDimension()/2 < lado*lado*lado)
                strassenAtSize = (1 << (int)(floor(log2(cbrt(cc->GetRingDimension()/2)))));
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
    KeyPair<Element> keys;
    std::set<int> rotations;
    std::shared_ptr<std::map<uint32_t, EvalKey<Element>>> sumRowsColsKeys;
    measureBlock("KeyGeneration", statsResults, [&](){
        keys = cc->KeyGen();
        cc->EvalMultKeyGen(keys.secretKey);
        rotations = rotStrassen(statsResults, lado, strassenAtSize, subAlgorithm);
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
                matMultStrassen(encC[r-1], encB, encC[r], statsResults, reiteraciones, lado, strassenAtSize, subAlgorithm, r, sumRowsColsKeys, fastRot);
            });
        } else {
            measureBlock("MatrixMultiplication", statsResults, [&](){
                matMultStrassen(encA, encB, encC[0], statsResults, reiteraciones, lado, strassenAtSize, subAlgorithm, r, sumRowsColsKeys, fastRot);
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