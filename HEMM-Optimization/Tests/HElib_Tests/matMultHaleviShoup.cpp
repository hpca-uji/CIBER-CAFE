#include <helib/helib.h>
#include <omp.h>
#include <map>

using namespace helib;

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

void matMultHaleviShoup(const std::vector<Ctxt>& encA,
                        const std::vector<Ctxt>& encB,
                        std::vector<Ctxt>& encC,
                        std::map<std::string, std::string>& statsResults,
                        size_t& reiteraciones,
                        int rowsM1 = 0,
                        const size_t reitActual=0,
                        bool useShift=false) {
    bool error = false;
    size_t slots = encA[0].getContext().getNSlots();
    if (rowsM1 == 0)
        rowsM1 = slots;
    int colsM1rowsM2 = encA.size(),
        colsM2 = encB.size();
    encC = std::vector<Ctxt>(colsM2, Ctxt(encB[0].getPubKey()));
    
    #if USE_OMP_TASKLOOP
    #pragma omp taskloop shared(encC)
    for (int i = 0; i < colsM2; i++)
        try {
            std::vector<Ctxt> temp = std::vector<Ctxt>(colsM1rowsM2, Ctxt(encB[0].getPubKey()));
            #pragma omp taskloop shared(temp)
    #else
    std::vector<std::vector<Ctxt>> temp = std::vector<std::vector<Ctxt>>(colsM2, std::vector<Ctxt>(colsM1rowsM2, Ctxt(encB[0].getPubKey())));
    #pragma omp parallel for collapse(2) shared(temp) schedule(dynamic)
    for (int i = 0; i < colsM2; i++)
    #endif
            for (int j = 0; j < colsM1rowsM2; j++)
                try {
                    int t_idx = mod(colsM2-j+i, colsM2);
                    Ctxt t = encB[t_idx];
                    if (j != 0) {
                        if(useShift) {
                            shift(t, -j);
                        } else {
                            rotate(t, -j);
                        }
                    }
                    if (slots != colsM1rowsM2 || mod(colsM2-colsM1rowsM2, colsM2) != 0) {
                        int lastPos = std::min(rowsM1, colsM1rowsM2-j);
                        if (!useShift && colsM1rowsM2+rowsM1-slots-lastPos > 0) {
                            t *= makePlaintext(encB[0].getContext(), slots, [lastPos](int ell){return (ell<lastPos);});
                        }
                        for (int k = 0; k*colsM1rowsM2+lastPos < rowsM1; k++) {
                            int tt_idx = mod(colsM2-j+i-(k+1)*(colsM2-colsM1rowsM2), colsM2);
                            Ctxt tt = encB[tt_idx];
                            int exceso = (k+1)*colsM1rowsM2+lastPos-slots;
                            if (!useShift && exceso > 0)
                                tt *= makePlaintext(encB[0].getContext(), slots, [colsM1rowsM2, exceso](int ell){return ell < colsM1rowsM2-exceso;});
                            if (k*colsM1rowsM2+lastPos != 0){
                                if(useShift) {
                                    shift(tt, (k*colsM1rowsM2+lastPos));
                                } else {
                                    rotate(tt, (k*colsM1rowsM2+lastPos));
                                }
                            }
                            t += tt;
                        }
                    }
                    #if USE_OMP_TASKLOOP
                    temp[j] = encA[j];
                    temp[j] *= t;
                    #else
                    temp[i][j] = encA[j];
                    temp[i][j] *= t;
                    #endif
                } catch (const std::exception& e) {
                    newException(e, statsResults, "MatrixMultiplication_" + to_string_precise(reitActual) + ".1[" + to_string_precise(i) + "][" + to_string_precise(j) + "]");
                    reiteraciones = reitActual;
                    error = true;
                }
    #if USE_OMP_TASKLOOP
            encC[i] = temp[0];
            for (int j = 1; j < temp.size(); j++)
                encC[i] += temp[j];
        } catch (const std::exception& e) {
            newException(e, statsResults, "MatrixMultiplication_" + to_string_precise(reitActual) + ".0[" + to_string_precise(i) + "]");
            reiteraciones = reitActual;
            error = true;
        }
    #else
    if (error)
        return;
    #pragma omp parallel for shared(encC) schedule(dynamic)
    for (int i = 0; i < colsM2; i++) {
        encC[i] = temp[i][0];
        for (int j = 1; j < temp[i].size(); j++)
            encC[i] += temp[i][j];
    }
    #endif

}

template <typename T>
int MatMultHaleviShoup(const Context& cc,
                       const std::vector<std::vector<T>>& A,
                       const std::vector<std::vector<T>>& B,
                       std::vector<std::vector<std::vector<T>>>& result,
                       size_t reiteraciones,
                       std::map<std::string, std::string>& statsResults,
                       bool useShift=false) {
    if (cc.getNSlots() < A.size() || cc.getNSlots() < B.size()) {
        newException(statsResults, "Not enough slots for the given matrices");
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
    std::vector<std::vector<T>> diagA, diagB;
    measureBlock("Preprocess", statsResults, [&](){
        diagA = encodeHaleviShoup(A);
        diagB = encodeHaleviShoup(B);
    });
    std::vector<PtxtArray> ptxA(diagA.size(), PtxtArray(cc)), ptxB(diagB.size(), PtxtArray(cc));
    #if USE_OMP_TASKLOOP
    #pragma omp parallel
    {
    #pragma omp single nowait
    {
    #endif
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
            ptxA[i] = PtxtArray(cc, diagA[i]);
        #if USE_OMP_TASKLOOP
        }
        #pragma omp task
        {
        #pragma omp taskloop shared(ptxB)
        #else
        #pragma omp parallel for shared(ptxB) schedule(dynamic)
        #endif
        for (size_t i = 0; i < diagB.size(); ++i)
            ptxB[i] = PtxtArray(cc, diagB[i]);
        #if USE_OMP_TASKLOOP
        }
        }
        #endif
    });
    std::vector<Ctxt> encA(ptxA.size(), Ctxt(publicKey)), encB(ptxB.size(), Ctxt(publicKey));
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
        for (size_t i = 0; i < ptxA.size(); ++i)
            ptxA[i].encrypt(encA[i]);
        #if USE_OMP_TASKLOOP
        }
        #pragma omp task
        {
        #pragma omp taskloop shared(encB)
        #else
        #pragma omp parallel for shared(encB) schedule(dynamic)
        #endif
        for (size_t i = 0; i < ptxB.size(); ++i)
            ptxB[i].encrypt(encB[i]);
        #if USE_OMP_TASKLOOP
        }
        }
        #endif
    });
    addToMap(statsResults, "ctxA_size", to_string_precise(encA.size()));
    std::vector<size_t> capacity = std::vector<size_t>(encA.size()),
                        errorBound = std::vector<size_t>(encA.size());
    #pragma omp parallel for shared(capacity, errorBound) schedule(dynamic)
    for (size_t i = 0; i < encA.size(); ++i) {
        capacity[i] = encA[i].capacity();
        errorBound[i] = encA[i].errorBound();
    }
    calculateStats("ctxA_capacity", statsResults, capacity);
    calculateStats("ctxA_errorBound", statsResults, errorBound);
    addToMap(statsResults, "ctxB_size", to_string_precise(encB.size()));
    capacity = std::vector<size_t>(encB.size());
    errorBound = std::vector<size_t>(encB.size());
    #pragma omp parallel for shared(capacity, errorBound) schedule(dynamic)
    for (size_t i = 0; i < encB.size(); ++i) {
        capacity[i] = encB[i].capacity();
        errorBound[i] = encB[i].errorBound();
    }
    calculateStats("ctxB_capacity", statsResults, capacity);
    calculateStats("ctxB_errorBound", statsResults, errorBound);

    std::vector<std::vector<Ctxt>> encC(reiteraciones);
    std::vector<std::vector<PtxtArray>> ptxC(reiteraciones);
    std::vector<std::vector<std::vector<T>>> decC(reiteraciones);
    for (size_t r = 0; r < reiteraciones; ++r) {
        #if ENABLE_DEBUG
        std::cout << "Starting Matrix Multiplication " << r+1 << "/" << reiteraciones << "..." << std::endl;
        #endif
        if (r > 0) {
            measureBlock("MatrixMultiplication", statsResults, [&](){
                matMultHaleviShoup(encC[r-1], encB, encC[r], statsResults, reiteraciones, A.size(), r, useShift);
            });
        } else {
            measureBlock("MatrixMultiplication", statsResults, [&](){
                matMultHaleviShoup(encA, encB, encC[0], statsResults, reiteraciones, A.size(), r, useShift);
            });
        }
        if(reiteraciones > r) {
            addToMap(statsResults, "ctxC_size", to_string_precise(encC[r].size()));
            capacity = std::vector<size_t>(encC[r].size());
            errorBound = std::vector<size_t>(encC[r].size());
            #pragma omp parallel for shared(capacity, errorBound) schedule(dynamic)
            for (size_t i = 0; i < encC[r].size(); ++i) {
                capacity[i] = encC[r][i].capacity();
                errorBound[i] = encC[r][i].errorBound();
            }
            calculateStats("ctxC_capacity", statsResults, capacity);
            calculateStats("ctxC_errorBound", statsResults, errorBound);
            measureBlock("Decrypt", statsResults, [&](){
                ptxC[r] = std::vector<PtxtArray>(encC[0].size(), PtxtArray(cc));
                #pragma omp parallel for shared(ptxC, statsResults, reiteraciones) schedule(dynamic)
                for (size_t i = 0; i < encC[0].size(); ++i)
                    try {
                        ptxC[r][i].decrypt(encC[r][i], secretKey);
                    } catch (const std::exception& e) {
                        newException(e, statsResults, "Decrypt_" + to_string_precise(r) + "[" + to_string_precise(i) + "]");
                        reiteraciones = r; // Stop further repetitions
                    }
            });
            if(reiteraciones > r) {
                measureBlock("Decode", statsResults, [&](){
                    decC[r] = std::vector<std::vector<T>>(encC[0].size());
                    #pragma omp parallel for shared(decC, statsResults, reiteraciones) schedule(dynamic)
                    for (size_t i = 0; i < encC[0].size(); ++i)
                        try {
                            ptxC[r][i].store(decC[r][i]);
                        } catch (const std::exception& e) {
                            newException(e, statsResults, "Decode_" + to_string_precise(r) + "[" + to_string_precise(i) + "]");
                            reiteraciones = r; // Stop further repetitions
                        }
                });
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