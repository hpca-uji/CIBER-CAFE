#include <iostream>
#include <chrono>

#include <helib/helib.h>
#include <helib/intraSlot.h>
#include "functional_units/functional_units.hpp"

using namespace helib;
using namespace std;

int main(void) {
  unsigned long p = 137, m = 3589, r = 1, bits = 30;
  unsigned long c = 2;
  Context context = helib::ContextBuilder<helib::BGV>()
    .m(m).p(p).r(r).bits(bits).c(c).build();
  SecKey secret_key(context);
  cout << "Security Level: " << context.securityLevel() << endl;
  secret_key.GenSecKey();
  addSome1DMatrices(secret_key);
  const PubKey& public_key = secret_key;
  const EncryptedArray& ea = context.getEA();
  long slots = ea.size();
  Ptxt<helib::BGV> tmp(context);
  Ctxt tmp_(public_key);

  vector<Ctxt> m1_, m2_, res_;
  int r1, r2, c2, i, j, k, i_r1, k_c2, i_c1;
  int dim;
  r1 = 4;
  r2 = 4;
  c2 = 4;
  dim = r1*r2;
  m1_.resize(dim, scratch);
  m2_.resize(dim, scratch);
  res_.resize(dim, scratch);
  for (i = 0; i<dim; i++) {
    for (int tmp_ij = 0; tmp_ij < slots; tmp_ij++) {
      tmp[tmp_ij] = 1+i;
    }
    public_key.Encrypt(m1_[i], tmp);
    for (int tmp_ij = 0; tmp_ij < slots; tmp_ij++) {
      tmp[tmp_ij] = 1+i;
    }
    public_key.Encrypt(m2_[i], tmp);
    for (int tmp_ij = 0; tmp_ij < slots; tmp_ij++) {
      tmp[tmp_ij] = 0;
    }
    public_key.Encrypt(res_[i], tmp);
  }
  auto start_timer = chrono::high_resolution_clock::now();
  for (i = 0; i<r1; i++) {
    i_r1 = i*r1;
    i_c1 = i*r2;
    for (j = 0; j<c2; j++) {
      for (k = 0; k<r2; k++) {
        k_c2 = k*c2;
        Ctxt tmp_0_(public_key);
        tmp_0_ = (m1_[(i_c1+k)]);
        tmp_0_ *= (m2_[(k_c2+j)]);
        res_[(i_r1+j)] += (tmp_0_);
      }
    }
  }
  auto stop_timer = chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<chrono::milliseconds>(stop_timer-start_timer);
  cout << "Time: " << duration.count() << "ms" << endl;
  for (i = 0; i<(r1*c2); i++) {
    secret_key.Decrypt(tmp, res_[i]);
    cout << "dec(res_[i]) = " << static_cast<long>(tmp[0]) << endl;
  }
  return 0;
}