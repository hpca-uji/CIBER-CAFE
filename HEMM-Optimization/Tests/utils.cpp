#include <map>
#include <chrono>
#include <iostream>
#include <fstream>
#include <sys/resource.h>
#include <unistd.h>
#include <algorithm>
#include <numeric>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <regex>
#include <sstream>
#include <iomanip>


enum Algorithm {
    Traditional,
    RowsXCols,
    RowsXColsSum,
    HaleviShoup,
    HaleviShoupFastRot,
    HaleviShoupShift,
    RHEMatMult,
    RHEMatMultFastRot,
    HEMatMult,
    HEMatMultFastRot,
    RizomiliotisTriakosia,
    RizomiliotisTriakosiaFastRot,
    RizomiliotisTriakosiaSums,
    RizomiliotisTriakosiaSumsFastRot,
    Strassen1x1x1,
    //StrassenHEMatMult,
    StrassenRizomiliotisTriakosia,
    StrassenRizomiliotisTriakosiaFastRot,
    StrassenRizomiliotisTriakosiaSums,
    StrassenRizomiliotisTriakosiaSumsFastRot
};

static const bool isRectangular(Algorithm alg) {
    return (alg <= RHEMatMultFastRot);
}

static const std::map<Algorithm, std::string> AlgorithmToString = {
    {Traditional, "Traditional"},
    {RowsXCols, "RowsXCols"},
    {RowsXColsSum, "RowsXColsSum"},
    {HaleviShoup, "HaleviShoup"},
    {HaleviShoupFastRot, "HaleviShoupFastRot"},
    {HaleviShoupShift, "HaleviShoupShift"},
    {RHEMatMult, "RHEMatMult"},
    {RHEMatMultFastRot, "RHEMatMultFastRot"},
    {HEMatMult, "HEMatMult"},
    {HEMatMultFastRot, "HEMatMultFastRot"},
    {RizomiliotisTriakosia, "RizomiliotisTriakosia"},
    {RizomiliotisTriakosiaFastRot, "RizomiliotisTriakosiaFastRot"},
    {RizomiliotisTriakosiaSums, "RizomiliotisTriakosiaSums"},
    {RizomiliotisTriakosiaSumsFastRot, "RizomiliotisTriakosiaSumsFastRot"},
    {Strassen1x1x1, "Strassen1x1x1"},
    //{StrassenHEMatMult, "StrassenHEMatMult"},
    {StrassenRizomiliotisTriakosia, "StrassenRizomiliotisTriakosia"},
    {StrassenRizomiliotisTriakosiaFastRot, "StrassenRizomiliotisTriakosiaFastRot"},
    {StrassenRizomiliotisTriakosiaSums, "StrassenRizomiliotisTriakosiaSums"},
    {StrassenRizomiliotisTriakosiaSumsFastRot, "StrassenRizomiliotisTriakosiaSumsFastRot"}
};

static const std::map<Algorithm, Algorithm> StrassenSubAlgorithm = {
    {Strassen1x1x1, Traditional},
    //{StrassenHEMatMult, HEMatMult},
    {StrassenRizomiliotisTriakosia, RizomiliotisTriakosia},
    {StrassenRizomiliotisTriakosiaFastRot, RizomiliotisTriakosiaFastRot},
    {StrassenRizomiliotisTriakosiaSums, RizomiliotisTriakosiaSums},
    {StrassenRizomiliotisTriakosiaSumsFastRot, RizomiliotisTriakosiaSumsFastRot}
};


// 1) Leer VmRSS desde /proc/self/status -> bytes
size_t getRSS_via_status_bytes() {
    std::ifstream f("/proc/self/status");
    std::string line;
    while (std::getline(f, line)) {
        if (line.rfind("VmRSS:", 0) == 0) {
            std::istringstream iss(line);
            std::string key, value, unit;
            iss >> key >> value >> unit; // key=="VmRSS:", value in kB
            try {
                size_t kB = std::stoull(value);
                return kB * 1024;
            } catch (...) {
                return 0;
            }
        }
    }
    return 0;
}

// 3) Sumar Pss desde /proc/self/smaps (proportional set size) -> bytes
// Nota: más lento (analiza todo smaps) pero útil si quieres medir coste real
size_t getPss_via_smaps_bytes() {
    std::ifstream f("/proc/self/smaps");
    if (!f) return 0;
    std::string line;
    size_t total_kB = 0;
    while (std::getline(f, line)) {
        if (line.rfind("Pss:", 0) == 0) {
            std::istringstream iss(line);
            std::string key;
            size_t value;
            std::string unit;
            iss >> key >> value >> unit; // value en kB
            total_kB += value;
        }
    }
    return total_kB * 1024;
}

int mod(int a, int b) {
    if (a >= 0) {
        return a % b;
    } else {
        return (a - b * (a / b) + b) % b;
    }
}

template <typename T>
std::string to_string_precise(const T value,
                              const int precision = 18) {
    std::ostringstream oss;
    if (std::is_integral<T>::value || value == static_cast<unsigned long long>(value)) {
        oss << value;
    } else {
        // Convert to double for magnitude checks to avoid ambiguous abs overloads
        double v = static_cast<double>(value);
        double mag = std::fabs(v);
        // Use scientific notation for very small/large numbers, fixed otherwise
        if (mag > 0 && (mag < 1e-4 || mag > 1e6) && mag < pow(10, precision)) {
            oss << std::scientific << std::setprecision(precision) << v;
        } else {
            oss << std::fixed << std::setprecision(0) << v;
        }
    }
    return oss.str();
}

void addToMap(std::map<std::string, std::string>& statsResults,
              const std::string& key,
              const std::string& value) {
    if (statsResults.find(key) == statsResults.end()) {
        statsResults[key] = value;
    } else {
        statsResults[key] += ";" + value;
    }
}

void newException(const std::exception& e,
                  std::map<std::string, std::string>& statsResults,
                  const std::string& context = "") {
    if (context.empty()){
        std::cerr << e.what() << std::endl;
        addToMap(statsResults, "exceptions", "\"" + std::string(e.what()) + "\"");
    } else {
        std::cerr << context << ": " << e.what() << std::endl;
        addToMap(statsResults, "exceptions", "\"" + context + ": " + std::string(e.what()) + "\"");
    }
}

void newException(std::map<std::string, std::string>& statsResults,
                  const std::string& context = "") {
    std::cerr << context << std::endl;
    addToMap(statsResults, "exceptions", "\"" + context + "\"");
}

template <typename F>
void measureBlock(const std::string label,
                  std::map<std::string, std::string>& statsResults,
                  const F func) {
    size_t page_size = static_cast<size_t>(sysconf(_SC_PAGESIZE));
    std::chrono::high_resolution_clock::time_point start;
    int duration;
    if (page_size <= 0) page_size = 4096;
    size_t ramBefore = 0, ramAfter = 0;
    ramBefore = getPss_via_smaps_bytes();
    start = std::chrono::high_resolution_clock::now();
    func();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count();
    ramAfter = getPss_via_smaps_bytes();
    addToMap(statsResults, "ram_" + label, to_string_precise(ramAfter - ramBefore));
    addToMap(statsResults, "dur_" + label, to_string_precise(duration));
    #if ENABLE_DEBUG
    std::cout << label << "\tTiempo:\t" << duration << " ms\n" << label << "\tRAM:\t" << (ramAfter - ramBefore) << " B" << std::endl;
    #endif
}

// Helper to detect if a type is a std::vector
template <typename T>
struct is_vector : std::false_type {};

template <typename T, typename Alloc>
struct is_vector<std::vector<T, Alloc>> : std::true_type {};

// Recursive flatten for nested vectors
template <typename T>
auto flatten(const std::vector<T>& vec) -> std::enable_if_t<!is_vector<T>::value, std::vector<T>> {
    // Base case: T is not a vector
    return vec;
}

template <typename T>
auto flatten(const std::vector<T>& vec) -> std::enable_if_t<is_vector<T>::value, std::vector<typename T::value_type>> {
    // Recursive case: T is a vector
    using InnerType = typename T::value_type;
    std::vector<InnerType> result;
    for (const auto& subvec : vec) {
        auto flat = flatten(subvec);
        result.insert(result.end(), flat.begin(), flat.end());
    }
    return result;
}

template <typename T>
std::vector<std::vector<T>> padMatrix(const std::vector<std::vector<T>>& matrix,
                                        const size_t newSize) {
    size_t oldRows = matrix.size();
    size_t oldCols = matrix[0].size();
    std::vector<std::vector<T>> paddedMatrix(newSize, std::vector<T>(newSize, 0));
    for (size_t i = 0; i < oldRows; ++i) {
        for (size_t j = 0; j < oldCols; ++j) {
            paddedMatrix[i][j] = matrix[i][j];
        }
    }
    return paddedMatrix;
}
// Compute percentile (linear interpolation) on a copy of the data
template <typename E>
double percentile(std::vector<E> data,
                  const double percent) {
    if (data.empty()) return 0.0;
    std::sort(data.begin(), data.end());
    if (percent <= 0) return data.front();
    if (percent >= 100) return data.back();
    double idx = (percent / 100.0) * (static_cast<double>(data.size()) - 1.0);
    size_t lo = static_cast<size_t>(std::floor(idx));
    size_t hi = static_cast<size_t>(std::ceil(idx));
    if (lo == hi || hi >= data.size()) return data[lo];
    double frac = idx - static_cast<double>(lo);
    return (1.0 - frac) * static_cast<double>(data[lo]) + frac * static_cast<double>(data[hi]);
}

template <typename E>
double calcStd(const std::vector<E>& data) {
    if (data.size() <= 1) return 0.0;
    double mean = std::accumulate(data.begin(), data.end(), 0.0) / static_cast<double>(data.size());
    double accum = 0.0;
    for (const auto& d : data)
        accum += pow(static_cast<double>(d) - mean, 2);
    return std::sqrt(accum / static_cast<double>(data.size()));
}

template <typename E>
void calculateStats(const std::string label,
                    std::map<std::string, std::string>& statsResults,
                    const std::vector<E>& data) {
    addToMap(statsResults, label + "_max", to_string_precise(*std::max_element(data.begin(), data.end())));
    addToMap(statsResults, label + "_min", to_string_precise(*std::min_element(data.begin(), data.end())));
    addToMap(statsResults, label + "_mean", to_string_precise(std::accumulate(data.begin(), data.end(), E(0)) / static_cast<E>(data.size())));
    addToMap(statsResults, label + "_median", to_string_precise(percentile(data, 50)));
    addToMap(statsResults, label + "_std", to_string_precise(calcStd(data)));
    addToMap(statsResults, label + "_p90", to_string_precise(percentile(data, 90)));
    addToMap(statsResults, label + "_p75", to_string_precise(percentile(data, 75)));
    addToMap(statsResults, label + "_p25", to_string_precise(percentile(data, 25)));
    addToMap(statsResults, label + "_p10", to_string_precise(percentile(data, 10)));
}

std::vector<std::vector<double>> readMatrix(const std::string& path,
                                            const size_t rows,
                                            const size_t cols) {
    std::vector<std::vector<double>> matrix(rows, std::vector<double>(cols));
    std::ifstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + path);
    }
    std::string line;
    size_t row = 0;
    while (row < rows && std::getline(file, line)) {
        std::istringstream ss(line);
        std::string val;
        size_t col = 0;
        while (col < cols && std::getline(ss, val, ',')) {
            matrix[row][col] = std::stod(val);
            ++col;
        }
        ++row;
    }
    if (row < rows) {
        throw std::runtime_error("Not enough rows in file: " + path);
    }
    return matrix;
}

void writeStats(const std::string& path,
                const std::vector<std::string>& headerOrig,
                const std::map<std::string, std::string>& stats) {
    std::vector<std::string> header;
    bool fileExists = false;
    std::ifstream infile(path);
    if (infile.good()) {
        fileExists = true;
        std::string firstLine;
        if (std::getline(infile, firstLine)) {
            std::istringstream ss(firstLine);
            std::string col;
            while (std::getline(ss, col, ',')) {
                header.push_back(col);
            }
        }
        for (const auto& col : headerOrig) {
            if (std::find(header.begin(), header.end(), col) == header.end()) {
                header.push_back(col);
            }
        }
    } else {
        header = headerOrig;
    }
    infile.close();

    std::ofstream file;
    if (fileExists) {
        file.open(path, std::ios::app);
    } else {
        file.open(path, std::ios::out);
    }
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file for writing: " + path);
    }
    if (!fileExists) {
        for (size_t i = 0; i < header.size(); ++i) {
            file << header[i];
            if (i < header.size() - 1) file << ",";
        }
        file << "\n";
    }
    for (size_t i = 0; i < header.size(); ++i) {
        auto it = stats.find(header[i]);
        if (it != stats.end()) {
            file << it->second; // write the value
        } else {
            file << ""; // write empty if key not found
        }
        if (i < header.size() - 1) file << ",";
    }
    file << "\n";
    file.close();
}

template <typename T>
std::vector<T> vec_from_pred(const uint n,
                             const std::function<bool(uint)> pred,
                             const uint nSlots = 0,
                             const T val = 1) {
    std::vector<T> vec(n, T(0));
    #if USE_OMP_TASKLOOP
    #pragma omp taskloop shared(vec)
    #endif
    for (uint i = 0; i < n; i++)
        if (pred(i))
            vec[i] = val;
    for (uint i = n; i < nSlots; i*=2)
        vec.insert(vec.end(), vec.begin(), vec.end());
    return vec;
}

int calcError(const std::vector<std::vector<double>>& A,
               const std::vector<std::vector<double>>& B,
               const std::vector<std::vector<std::vector<double>>>& result,
               size_t& reiteraciones,
               std::map<std::string, std::string>& statsResults){
    #if ENABLE_DEBUG
    std::cout << "Result:" << std::endl;
    for (size_t r = 0; r < result.size(); ++r) {
        std::cout << "Iteration " << r+1 << ":" << std::endl;
        for (size_t i = 0; i < result[r].size(); ++i) {
            for (size_t j = 0; j < result[r][i].size(); ++j) {
                std::cout << result[r][i][j] << " ";
            }
            std::cout << std::endl;
        }
    }
    std::cout << "Starting Check Error..." << std::endl;
    #endif
    reiteraciones = result.size(); // En caso de error, result puede tener menos iteraciones
    std::vector<std::vector<std::vector<double>>> expected(reiteraciones, std::vector<std::vector<double>>(A.size(), std::vector<double>(B[0].size(), 0.0))),
                                                errorAbs(reiteraciones, std::vector<std::vector<double>>(A.size(), std::vector<double>(B[0].size(), 0.0))),
                                                errorRel(reiteraciones, std::vector<std::vector<double>>(A.size(), std::vector<double>(B[0].size(), 0.0))),
                                                errorSmape(reiteraciones, std::vector<std::vector<double>>(A.size(), std::vector<double>(B[0].size(), 0.0)));
    #pragma omp parallel for collapse(2) shared(expected, errorAbs, errorRel, errorSmape) schedule(dynamic)
    for (size_t i = 0; i < A.size(); ++i)
        for (size_t j = 0; j < B[0].size(); ++j) {
            for (size_t k = 0; k < A[0].size(); ++k)
                expected[0][i][j] += A[i][k] * B[k][j];
            errorAbs[0][i][j] = std::abs(expected[0][i][j] - result[0][i][j]);
            errorRel[0][i][j] = std::abs(errorAbs[0][i][j] / expected[0][i][j]);
            errorSmape[0][i][j] = errorAbs[0][i][j] / ((std::abs(expected[0][i][j]) + std::abs(result[0][i][j])) / 2.0);
        }
    calculateStats("errorAbs", statsResults, flatten(errorAbs[0]));
    calculateStats("errorRel", statsResults, flatten(errorRel[0]));
    calculateStats("errorSmape", statsResults, flatten(errorSmape[0]));
    for (size_t r = 1; r < reiteraciones; ++r) {
        #pragma omp parallel for collapse(2) shared(expected, errorAbs, errorRel, errorSmape) schedule(dynamic)
        for (size_t i = 0; i < A.size(); ++i)
            for (size_t j = 0; j < B[0].size(); ++j) {
                for (size_t k = 0; k < A[0].size(); ++k)
                    expected[r][i][j] += expected[r-1][i][k] * B[k][j];
                errorAbs[r][i][j] = std::abs(expected[r][i][j] - result[r][i][j]);
                errorRel[r][i][j] = std::abs(errorAbs[r][i][j] / expected[r][i][j]);
                errorSmape[r][i][j] = errorAbs[r][i][j] / ((std::abs(expected[r][i][j]) + std::abs(result[r][i][j])) / 2.0);
            }
        calculateStats("errorAbs", statsResults, flatten(errorAbs[r]));
        calculateStats("errorRel", statsResults, flatten(errorRel[r]));
        calculateStats("errorSmape", statsResults, flatten(errorSmape[r]));
    }
    #if ENABLE_DEBUG
    for (size_t r = 0; r < errorAbs.size(); ++r) {
        std::cout << "errorAbs " << r+1 << ":" << std::endl;
        for (size_t i = 0; i < errorAbs[r].size(); ++i) {
            for (size_t j = 0; j < errorAbs[r][i].size(); ++j)
                std::cout << errorAbs[r][i][j] << " ";
            std::cout << std::endl;
        }
    }
    for (size_t r = 0; r < errorRel.size(); ++r) {
        std::cout << "errorRel " << r+1 << ":" << std::endl;
        for (size_t i = 0; i < errorRel[r].size(); ++i) {
            for (size_t j = 0; j < errorRel[r][i].size(); ++j)
                std::cout << errorRel[r][i][j] << " ";
            std::cout << std::endl;
        }
    }
    for (size_t r = 0; r < errorSmape.size(); ++r) {
        std::cout << "errorSmape " << r+1 << ":" << std::endl;
        for (size_t i = 0; i < errorSmape[r].size(); ++i) {
            for (size_t j = 0; j < errorSmape[r][i].size(); ++j)
                std::cout << errorSmape[r][i][j] << " ";
            std::cout << std::endl;
        }
    }
    std::cout << "errorAbs: " << std::endl;
    for (auto &e :flatten(errorAbs[0]))
        std::cout << e << " ";
    std::cout << std::endl << "errorRel: " << std::endl;
    for (auto &e :flatten(errorRel[0]))
        std::cout << e << " ";
    std::cout << std::endl << "errorSmape: " << std::endl;
    for (auto &e :flatten(errorSmape[0]))
        std::cout << e << " ";
    std::cout << std::endl;
    auto data = flatten(errorSmape[0]);
    std::cout << "MaxError " << *std::max_element(data.begin(), data.end()) << std::endl;
    #endif
    for (const std::vector<std::vector<double>>& mtrx : errorRel)
        for (const std::vector<double>& row : mtrx)
            for (double val : row)
                if (val >= 1.0)
                    return 1;
    return 0;
}