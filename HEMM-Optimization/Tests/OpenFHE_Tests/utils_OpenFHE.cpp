#include <openfhe.h>

using namespace lbcrypto;

template <class Element>
Plaintext makePlaintext(const CryptoContext<Element>& cc, const bool isCKKS, const uint n, std::function<bool(uint)> pred, uint32_t slots=0) {
    if (isCKKS)
        return cc->MakeCKKSPackedPlaintext(vec_from_pred<double>(n, pred), 1, 0, nullptr, slots);
    else
        return cc->MakePackedPlaintext(vec_from_pred<int64_t>(n, pred));
}
