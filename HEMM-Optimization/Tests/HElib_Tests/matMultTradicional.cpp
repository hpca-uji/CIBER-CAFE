#include <helib/helib.h>
#include <omp.h>
#include <map>

using namespace helib;

void matMultTradicional(const std::vector<std::vector<Ctxt>>& encA,
                        const std::vector<std::vector<Ctxt>>& encB,
                        std::vector<std::vector<Ctxt>>& encC,
                        std::map<std::string, std::string>& statsResults,
                        size_t& reiteraciones,
                        const size_t reitActual=0) {
    encC = std::vector<std::vector<Ctxt>>(encA.size(), std::vector<Ctxt>(encB[0].size(), Ctxt(encA[0][0].getPubKey())));
    #pragma omp parallel for collapse(2) shared(encC, statsResults, reiteraciones) schedule(dynamic)
    for (size_t i = 0; i < encA.size(); ++i)
        for (size_t j = 0; j < encB[0].size(); ++j)
            try {
                encC[i][j] = encA[i][0];
                encC[i][j] *= encB[0][j];
                for (size_t k = 1; k < encA[0].size(); ++k) {
                    Ctxt tempCtxt = encA[i][k];
                    tempCtxt *= encB[k][j];
                    encC[i][j] += tempCtxt;
                }
            } catch (const std::exception& e) {
                newException(e, statsResults, "MatrixMultiplication_" + to_string_precise(reitActual) + "[" + to_string_precise(i) + "][" + to_string_precise(j) + "]");
                reiteraciones = reitActual; // Stop further repetitions
            }
}

template <typename T = double>
int MatMultTradicional(const Context& cc,
                       const std::vector<std::vector<T>>& A,
                       const std::vector<std::vector<T>>& B,
                       std::vector<std::vector<std::vector<T>>>& result,
                       size_t reiteraciones,
                       std::map<std::string, std::string>& statsResults) {
    #if ENABLE_DEBUG
    std::cout << "Starting Key Generation..." << std::endl;
    #endif
    SecKey secretKey(cc);
    measureBlock("KeyGeneration", statsResults, [&](){
        secretKey.GenSecKey();
    });
    const PubKey& publicKey = secretKey;

    #if ENABLE_DEBUG
    std::cout << "Starting Encryption..." << std::endl;
    #endif
    std::vector<std::vector<Ctxt>> encA(A.size(), std::vector<Ctxt>(A[0].size(), Ctxt(publicKey))), encB(B.size(), std::vector<Ctxt>(B[0].size(), Ctxt(publicKey)));
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
    std::vector<std::vector<PtxtArray>> ptxA(A.size(), std::vector<PtxtArray>(A[0].size(), PtxtArray(cc))), ptxB(B.size(), std::vector<PtxtArray>(B[0].size(), PtxtArray(cc)));
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
                ptxA[i][j] = PtxtArray(cc, preA[i][j]);
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
                ptxB[i][j] = PtxtArray(cc, preB[i][j]);
        #if USE_OMP_TASKLOOP
        }
        }
        #endif
    });
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
                ptxA[i][j].encrypt(encA[i][j]);
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
                ptxB[i][j].encrypt(encB[i][j]);
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
    std::vector<size_t> capacity = std::vector<size_t>(encA.size() * encA[0].size()),
                        errorBound = std::vector<size_t>(encA.size() * encA[0].size());
    #pragma omp parallel for collapse(2) shared(capacity, errorBound) schedule(dynamic)
    for (size_t i = 0; i < encA.size(); ++i)
        for (size_t j = 0; j < encA[0].size(); ++j) {
            capacity[i * encA[0].size() + j] = encA[i][j].capacity();
            errorBound[i * encA[0].size() + j] = encA[i][j].errorBound();
        }
    calculateStats("ctxA_capacity", statsResults, capacity);
    calculateStats("ctxA_errorBound", statsResults, errorBound);
    addToMap(statsResults, "ctxB_size", to_string_precise(encB.size() * encB[0].size()));
    capacity = std::vector<size_t>(encB.size() * encB[0].size());
    errorBound = std::vector<size_t>(encB.size() * encB[0].size());
    #pragma omp parallel for collapse(2) shared(capacity, errorBound) schedule(dynamic)
    for (size_t i = 0; i < encB.size(); ++i)
        for (size_t j = 0; j < encB[0].size(); ++j) {
            capacity[i * encB[0].size() + j] = encB[i][j].capacity();
            errorBound[i * encB[0].size() + j] = encB[i][j].errorBound();
        }
    calculateStats("ctxB_capacity", statsResults, capacity);
    calculateStats("ctxB_errorBound", statsResults, errorBound);

    std::vector<std::vector<std::vector<Ctxt>>> encC(reiteraciones);
    std::vector<std::vector<std::vector<PtxtArray>>> ptxC(reiteraciones);
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
            capacity = std::vector<size_t>(encA.size() * encB[0].size());
            errorBound = std::vector<size_t>(encA.size() * encB[0].size());
            addToMap(statsResults, "ctxC_size", to_string_precise(encC[r].size() * encC[r][0].size()));
            #pragma omp parallel for collapse(2) shared(capacity, errorBound) schedule(dynamic)
            for (size_t i = 0; i < encA.size(); ++i)
                for (size_t j = 0; j < encB[0].size(); ++j) {
                    capacity[i * encB[0].size() + j] = encC[r][i][j].capacity();
                    errorBound[i * encB[0].size() + j] = encC[r][i][j].errorBound();
                }
            calculateStats("ctxC_capacity", statsResults, capacity);
            calculateStats("ctxC_errorBound", statsResults, errorBound);
            measureBlock("Decrypt", statsResults, [&](){
                ptxC[r] = std::vector<std::vector<PtxtArray>>(encC[0].size(), std::vector<PtxtArray>(encC[0][0].size(), PtxtArray(cc)));
                #pragma omp parallel for collapse(2) shared(ptxC, statsResults, reiteraciones) schedule(dynamic)
                for (size_t i = 0; i < encC[0].size(); ++i)
                    for (size_t j = 0; j < encC[0][0].size(); ++j)
                        try {
                            ptxC[r][i][j].decrypt(encC[r][i][j], secretKey);
                        } catch (const std::exception& e) {
                            newException(e, statsResults, "Decrypt_" + to_string_precise(r) + "[" + to_string_precise(i) + "][" + to_string_precise(j) + "]");
                            reiteraciones = r; // Stop further repetitions
                        }
            });
            if(reiteraciones > r) {
                measureBlock("Decode", statsResults, [&](){
                    decC[r] = std::vector<std::vector<std::vector<T>>>(encC[0].size(), std::vector<std::vector<T>>(encC[0][0].size()));
                    #pragma omp parallel for collapse(2) shared(decC, statsResults, reiteraciones) schedule(dynamic)
                    for (size_t i = 0; i < encC[0].size(); ++i)
                        for (size_t j = 0; j < encC[0][0].size(); ++j)
                            try {
                                ptxC[r][i][j].store(decC[r][i][j]);
                            } catch (const std::exception& e) {
                                newException(e, statsResults, "Decode_" + to_string_precise(r) + "[" + to_string_precise(i) + "][" + to_string_precise(j) + "]");
                                reiteraciones = r; // Stop further repetitions
                            }
                });
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
