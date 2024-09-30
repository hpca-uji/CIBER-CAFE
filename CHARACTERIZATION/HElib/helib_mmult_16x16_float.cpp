#include <iostream>
#include <chrono>

#include <iomanip>

#include <helib/helib.h>
#include <helib/intraSlot.h>
#include "../functional_units/functional_units.hpp"

using namespace helib;
using namespace std;

int main(void) {
  Context context = helib::ContextBuilder<helib::CKKS>()
    .m(16384).precision(20).bits(120).c(3).build();
  SecKey secret_key(context);
  cout << "Security Level: " << context.securityLevel() << endl;
  secret_key.GenSecKey();
  addSome1DMatrices(secret_key);
  const PubKey& public_key = secret_key;
  long slots = context.getNSlots();
  vector<std::complex<double>> tmp(slots);
  PtxtArray ptxt(context);
  Ctxt tmp_(public_key);

  vector<Ctxt> m1_, m2_, res_;
  Ctxt t3mp_(public_key);
  int r1, r2, c2, i, j, k, i_r1, k_c2, i_c1;
  int m;
  r1 = 16;
  r2 = 16;
  c2 = 16;
  m = r1*r1;
  m1_.resize(m, tmp_);
  m2_.resize(m, tmp_);
  res_.resize(m, tmp_);
  for (i = 0; i<m; i += 1.0) {
    for (int tmp_ij = 0; tmp_ij < slots; tmp_ij++) {
      tmp[tmp_ij] = std::complex<double>(1.0+i, 0);
    }
    ptxt = tmp; // encode
    ptxt.encrypt(m1_[i]);
;
    for (int tmp_ij = 0; tmp_ij < slots; tmp_ij++) {
      tmp[tmp_ij] = std::complex<double>(1.0+i, 0);
    }
    ptxt = tmp; // encode
    ptxt.encrypt(m2_[i]);
    for (int tmp_ij = 0; tmp_ij < slots; tmp_ij++) {
      tmp[tmp_ij] = std::complex<double>(0.0, 0);
    }
    ptxt = tmp; // encode
    ptxt.encrypt(res_[i]);
  }
  auto start_timer = chrono::high_resolution_clock::now();
  for (i = 0; i<r1; i += 1.0) {
    i_r1 = i*r1;
    i_c1 = i*r2;
    for (j = 0; j<c2; j += 1.0) {
      for (k = 0; k<r2; k += 1.0) {
        k_c2 = k*c2;
        Ctxt tmp_0_(public_key);
        tmp_0_ = (m1_[(i_c1+k)]);
        tmp_0_ *= (m2_[(k_c2+j)]);
        t3mp_ = (tmp_0_);
        res_[(i_r1+j)] += t3mp_;
      }
    }
  }
  auto stop_timer = chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<chrono::milliseconds>(stop_timer-start_timer);
  cout << "Time: " << duration.count() << "ms" << endl;
  for (i = 0; i<(r1*c2); i += 1.0) {
    ptxt.decryptComplex(res_[i], secret_key);
    ptxt.store(tmp);
    cout << "dec(res_[i]) = " << fixed << setprecision(1) << real(tmp[0]) << endl;
  }
  return 0;
}