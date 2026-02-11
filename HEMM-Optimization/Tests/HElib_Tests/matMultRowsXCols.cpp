#include <helib/helib.h>
#include <omp.h>
#include <map>

using namespace helib;

void matMultRowsXCols(const std::vector<Ctxt>& encA,
                      const std::vector<Ctxt>& encB,
                      std::vector<Ctxt>& encC,
                      std::map<std::string, std::string>& statsResults,
                      size_t& reiteraciones,
                      size_t k = 0,
                      const size_t reitActual=0,
                      const bool sum = false) {
    bool error = false;
    const PubKey& publicKey = encA[0].getPubKey();
    if (k==0)
        k = encA[0].getContext().getNSlots();
    std::vector<Ctxt> encTemp = std::vector<Ctxt>(encA.size()*encB.size(), Ctxt(publicKey));
    encC = std::vector<Ctxt>(encA.size(), Ctxt(publicKey));
    std::vector<PtxtArray> ptxtMasks = std::vector<PtxtArray>(encB.size(), PtxtArray(encA[0].getContext()));
    #if USE_OMP_TASKLOOP
    #pragma omp taskloop shared(ptxtMasks)
    #else
    #pragma omp parallel for shared(ptxtMasks) schedule(dynamic)
    #endif
    for (size_t j = 0; j < encB.size(); j++)
        ptxtMasks[j] = makePlaintext(encA[0].getContext(), encB.size(), [j](int ell){return ell == j;});
    #if USE_OMP_TASKLOOP
    #pragma omp taskloop collapse(2) shared(encTemp, reiteraciones, statsResults, error)
    #else
    #pragma omp parallel for collapse(2) shared(encTemp, reiteraciones, statsResults, error) schedule(dynamic)
    #endif
    for (size_t i = 0; i < encA.size(); i++)
        for (size_t j = 0; j < encB.size(); j++)
            try {
                encTemp[i*encB.size()+j] = encA[i];
                encTemp[i*encB.size()+j] *= encB[j];
                if (sum) {
                    totalSums(encTemp[i*encB.size()+j]);
                } else {
                    for (size_t kk = 0; (1<<kk) < k; ++kk) {
                        Ctxt temp = encTemp[i*encB.size()+j];
                        rotate(temp, (1<<kk)*(j&(1<<kk) ? 1 : -1));
                        encTemp[i*encB.size()+j] += temp;
                    }
                }
                encTemp[i*encB.size()+j] *= ptxtMasks[j];
            } catch (const std::exception& e) {
                newException(e, statsResults, "MatrixMultiplication_" + to_string_precise(reitActual) + ".1[" + to_string_precise(i) + "][" + to_string_precise(j) + "]");
                reiteraciones = reitActual;
                error = true;
            }
    if (error)
        return;
    #if USE_OMP_TASKLOOP
    #pragma omp taskloop shared(encC, reiteraciones, statsResults)
    #else
    #pragma omp parallel for shared(encC, reiteraciones, statsResults) schedule(dynamic)
    #endif
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

template <typename T = double>
int MatMultRowsXCols(const Context& cc,
                    const std::vector<std::vector<T>>& A,
                    const std::vector<std::vector<T>>& B,
                    std::vector<std::vector<std::vector<T>>>& result,
                    size_t reiteraciones,
                    std::map<std::string, std::string>& statsResults,
                    const bool sum = false) {
    if (cc.getNSlots() < A[0].size() || cc.getNSlots() < B.size()) {
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
    std::vector<std::vector<T>> transB;
    measureBlock("Preprocess", statsResults, [&](){
        transB = std::vector<std::vector<T>>(B[0].size(), std::vector<T>(B.size()));
        #pragma omp parallel for collapse(2) shared(transB) schedule(dynamic)
        for (size_t i = 0; i < B.size(); ++i)
            for (size_t j = 0; j < B[0].size(); ++j)
                transB[j][i] = B[i][j];
    });
    std::vector<PtxtArray> ptxA(A.size(), PtxtArray(cc)), ptxB(transB.size(), PtxtArray(cc));
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
        for (size_t i = 0; i < A.size(); ++i)
            ptxA[i] = PtxtArray(cc, A[i]);
        #if USE_OMP_TASKLOOP
        }
        #pragma omp task
        {
        #pragma omp taskloop shared(ptxB)
        #else
        #pragma omp parallel for shared(ptxB) schedule(dynamic)
        #endif
        for (size_t i = 0; i < transB.size(); ++i)
            ptxB[i] = PtxtArray(cc, transB[i]);
        #if USE_OMP_TASKLOOP
        }
        }
        #endif
    });
    std::vector<Ctxt> encA(A.size(), Ctxt(publicKey)), encB(transB.size(), Ctxt(publicKey));
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
            ptxA[i].encrypt(encA[i]);
        #if USE_OMP_TASKLOOP
        }
        #pragma omp task
        {
        #pragma omp taskloop shared(encB)
        #else
        #pragma omp parallel for shared(encB) schedule(dynamic)
        #endif
        for (size_t i = 0; i < transB.size(); ++i)
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
                matMultRowsXCols(encC[r-1], encB, encC[r], statsResults, reiteraciones, A[0].size(), r, sum);
            });
        } else {
            measureBlock("MatrixMultiplication", statsResults, [&](){
                matMultRowsXCols(encA, encB, encC[0], statsResults, reiteraciones, A[0].size(), 0, sum);
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
                ptxC[r] = std::vector<PtxtArray>(encC[r].size(), PtxtArray(cc));
                #pragma omp parallel for shared(ptxC, statsResults, reiteraciones) schedule(dynamic)
                for (size_t i = 0; i < encC[r].size(); ++i)
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
