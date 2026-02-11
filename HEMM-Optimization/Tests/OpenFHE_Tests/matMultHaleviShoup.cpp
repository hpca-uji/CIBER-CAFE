#include <openfhe.h>
#include <omp.h>
#include <map>

using namespace lbcrypto;

std::set<int> rotHaleviShoup(int rowsM1, int colsM1rowsM2, int colsM2) {
    std::set<int> ret;
    for (int j = 0; j < colsM1rowsM2; j++) {
        if (j != 0)
            ret.insert(j);
        int lastPos = std::min(rowsM1, colsM1rowsM2-j);
        for (int k = 0; k*colsM1rowsM2+lastPos < rowsM1; k++) {
            if (k*colsM1rowsM2+lastPos != 0)
                ret.insert(-(k*colsM1rowsM2+lastPos));
        }
    }
    return ret;
}

template<typename T>
std::vector<std::vector<T>> encodeHaleviShoup(const std::vector<std::vector<T>>& x) {
    if (x.empty() || x[0].empty()) return {};
    size_t rows = x.size(), cols = x[0].size();
    std::vector<std::vector<T>> ret(cols, std::vector<T>(rows));
    #if USE_OMP_TASKLOOP
    #pragma omp taskloop collapse(2) shared(ret)
    #else
    #pragma omp parallel for collapse(2) shared(ret) schedule(dynamic)
    #endif
    for (size_t i = 0; i < cols; i++)
        for (size_t j = 0; j < rows; j++)
            ret[i][j] = x[j][mod(j+i, cols)];
    return ret;
}

template<typename T>
std::vector<std::vector<T>> decodeHaleviShoup(const std::vector<std::vector<T>>& x, size_t rows) {
    if (x.empty() || x[0].empty()) return {};
    size_t cols = x.size();
    std::vector<std::vector<T>> ret(rows, std::vector<T>(cols));
    #pragma omp parallel for collapse(2) shared(ret) schedule(dynamic)
    for (size_t i = 0; i < cols; i++)
        for (size_t j = 0; j < rows; j++)
            ret[j][mod(i+j, cols)] = x[i][j];
    return ret;
}

template <class Element>
void matMultHaleviShoup(const std::vector<Ciphertext<Element>>& encA,
             const std::vector<Ciphertext<Element>>& encB,
             std::vector<Ciphertext<Element>>& encC,
             std::map<std::string, std::string>& statsResults,
             size_t& reiteraciones, int rowsM1 = 0,
             const size_t reitActual=0, bool fastRot=false) {
    bool error = false;
    CryptoContext<Element> cc = encA[0]->GetCryptoContext();
    bool isCKKS = (encA[0]->GetEncodingType() == CKKS_PACKED_ENCODING);
    size_t slots = encA[0]->GetSlots();
    if (rowsM1 == 0)
        rowsM1 = slots;
    int colsM1rowsM2 = encA.size(),
        colsM2 = encB.size();
    encC = std::vector<Ciphertext<Element>>(colsM2);
    
    std::vector<std::shared_ptr<std::vector<Element>>> precomRot;
    if (fastRot) {
        precomRot = std::vector<std::shared_ptr<std::vector<Element>>>(encB.size());
        #pragma omp parallel for shared(precomRot)
        for (int i = 0; i < encB.size(); i++)
            precomRot[i] = cc->EvalFastRotationPrecompute(encB[i]);
    }

    #if USE_OMP_TASKLOOP
    #pragma omp taskloop shared(encC)
    for (int i = 0; i < colsM2; i++)
        try {
            std::vector<Ciphertext<Element>> temp = std::vector<Ciphertext<Element>>(colsM1rowsM2);
            #pragma omp taskloop shared(temp)
    #else
    std::vector<std::vector<Ciphertext<Element>>> temp = std::vector<std::vector<Ciphertext<Element>>>(colsM2, std::vector<Ciphertext<Element>>(colsM1rowsM2));
    #pragma omp parallel for collapse(2) shared(temp) schedule(dynamic)
    for (int i = 0; i < colsM2; i++)
    #endif
            for (int j = 0; j < colsM1rowsM2; j++)
                try {
                    int t_idx = mod(colsM2-j+i, colsM2);
                    Ciphertext<Element> t = encB[t_idx];
                    if (j != 0) {
                        if(fastRot) {
                            t = cc->EvalFastRotation(t, j, cc->GetCyclotomicOrder(), precomRot[t_idx]);
                        } else {
                            t = cc->EvalRotate(t, j);
                        }
                    }
                    if (slots != colsM1rowsM2 || mod(colsM2-colsM1rowsM2, colsM2) != 0) {
                        int lastPos = std::min(rowsM1, colsM1rowsM2-j);
                        if (colsM1rowsM2+rowsM1-slots-lastPos > 0) {
                            t = cc->EvalMult(t, makePlaintext(cc, isCKKS, slots, [lastPos](int ell){return (ell<lastPos);}));
                        }
                        for (int k = 0; k*colsM1rowsM2+lastPos < rowsM1; k++) {
                            int tt_idx = mod(colsM2-j+i-(k+1)*(colsM2-colsM1rowsM2), colsM2);
                            Ciphertext<Element> tt = encB[tt_idx];
                            int exceso = (k+1)*colsM1rowsM2+lastPos-slots;
                            if (exceso > 0)
                                tt = cc->EvalMult(tt, makePlaintext(cc, isCKKS, slots, [colsM1rowsM2, exceso](int ell){return ell < colsM1rowsM2-exceso;}));
                            if (k*colsM1rowsM2+lastPos != 0){
                                if(fastRot) {
                                    tt = cc->EvalFastRotation(tt, -(k*colsM1rowsM2+lastPos), cc->GetCyclotomicOrder(), precomRot[t_idx]);
                                } else {
                                    tt = cc->EvalRotate(tt, -(k*colsM1rowsM2+lastPos));
                                }
                            }
                            t = cc->EvalAdd(t, tt);
                        }
                    }
                    #if USE_OMP_TASKLOOP
                    temp[j] = cc->EvalMult(encA[j], t);
                    #else
                    temp[i][j] = cc->EvalMult(encA[j], t);
                    #endif
                } catch (const std::exception& e) {
                    newException(e, statsResults, "MatrixMultiplication_" + to_string_precise(reitActual) + ".1[" + to_string_precise(i) + "][" + to_string_precise(j) + "]");
                    reiteraciones = reitActual;
                    error = true;
                }
    #if USE_OMP_TASKLOOP
            encC[i] = cc->EvalAddMany(temp);
        } catch (const std::exception& e) {
            newException(e, statsResults, "MatrixMultiplication_" + to_string_precise(reitActual) + ".0[" + to_string_precise(i) + "]");
            reiteraciones = reitActual;
            error = true;
        }
    #else
    if (error)
        return;
    #pragma omp parallel for shared(encC) schedule(dynamic)
    for (int i = 0; i < colsM2; i++)
        encC[i] = cc->EvalAddMany(temp[i]);
    #endif

}

template <class Element, typename T>
int MatMultHaleviShoup(const CryptoContext<Element>& cc,
                       const std::vector<std::vector<T>>& A,
                       const std::vector<std::vector<T>>& B,
                       std::vector<std::vector<std::vector<T>>>& result,
                       size_t reiteraciones,
                       std::map<std::string, std::string>& statsResults,
                       bool fastRot=false) {
    #if ENABLE_DEBUG
    std::cout << "Starting Key Generation..." << std::endl;
    #endif
    KeyPair<Element> keys;
    std::set<int> rotations;
    measureBlock("KeyGeneration", statsResults, [&](){
        keys = cc->KeyGen();
        cc->EvalMultKeyGen(keys.secretKey);
        rotations = rotHaleviShoup(A.size(), A[0].size(), B[0].size());
        cc->EvalRotateKeyGen(keys.secretKey, std::vector<int>(rotations.begin(), rotations.end()));
    });
    #if ENABLE_DEBUG
    std::cout << "Rotations(HaleviShoup): ";
    for (int r : rotations)
        std::cout << r << ", ";
    std::cout << std::endl;
    #endif

    #if ENABLE_DEBUG
    std::cout << "Starting Encryption..." << std::endl;
    #endif
    std::vector<std::vector<T>> diagA, diagB;
    measureBlock("Preprocess", statsResults, [&](){
        diagA = encodeHaleviShoup(A);
        diagB = encodeHaleviShoup(B);
    });
    std::vector<Plaintext> ptxA(diagA.size()), ptxB(diagB.size());
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
            for (size_t i = 0; i < diagA.size(); ++i)
                ptxA[i] = cc->MakeCKKSPackedPlaintext(diagA[i]);
            #if USE_OMP_TASKLOOP
            }
            #pragma omp task
            {
            #pragma omp taskloop shared(ptxB)
            #else
            #pragma omp parallel for shared(ptxB) schedule(dynamic)
            #endif
            for (size_t i = 0; i < diagB.size(); ++i)
                ptxB[i] = cc->MakeCKKSPackedPlaintext(diagB[i]);
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
            for (size_t i = 0; i < diagA.size(); ++i)
                ptxA[i] = cc->MakePackedPlaintext(diagA[i]);
            #if USE_OMP_TASKLOOP
            }
            #pragma omp task
            {
            #pragma omp taskloop shared(ptxB)
            #else
            #pragma omp parallel for shared(ptxB) schedule(dynamic)
            #endif
            for (size_t i = 0; i < diagB.size(); ++i)
                ptxB[i] = cc->MakePackedPlaintext(diagB[i]);
            #if USE_OMP_TASKLOOP
            }
            }
            #endif
        });
    }
    std::vector<Ciphertext<Element>> encA(diagA.size()), encB(diagB.size());
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
        for (size_t i = 0; i < diagA.size(); ++i)
            encA[i] = cc->Encrypt(keys.publicKey, ptxA[i]);
        #if USE_OMP_TASKLOOP
        }
        #pragma omp task
        {
        #pragma omp taskloop shared(encB)
        #else
        #pragma omp parallel for shared(encB) schedule(dynamic)
        #endif
        for (size_t i = 0; i < diagB.size(); ++i)
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
                matMultHaleviShoup(encC[r-1], encB, encC[r], statsResults, reiteraciones, A.size(), r, fastRot);
            });
        } else {
            measureBlock("MatrixMultiplication", statsResults, [&](){
                matMultHaleviShoup(encA, encB, encC[0], statsResults, reiteraciones, A.size(), r, fastRot);
            });
        }
        if(reiteraciones > r) {
            #if ENABLE_DEBUG
            std::cout << "Matrix Multiplication " << r+1 << "/" << reiteraciones << " completed." << std::endl;
            #endif
            std::vector<size_t> level = std::vector<size_t>(encC[0].size());
            #pragma omp parallel for shared(level) schedule(dynamic)
            for (size_t i = 0; i < encC[0].size(); ++i)
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
                        //#pragma omp parallel for shared(decC) schedule(dynamic)
                        //for (size_t i = 0; i < decC[0].size(); ++i)
                        //    decC[r][i] = std::vector<T>(decC[r][i].begin(), decC[r][i].begin() + A.size());
                        result[r] = decodeHaleviShoup(decC[r], A.size());
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