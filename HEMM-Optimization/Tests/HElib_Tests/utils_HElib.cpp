#include <helib/helib.h>

using namespace helib;

PtxtArray makePlaintext(const Context& cc, const uint n, std::function<bool(uint)> pred, uint32_t slots=0) {
    return PtxtArray(cc, vec_from_pred<double>(n, pred, slots));
}
