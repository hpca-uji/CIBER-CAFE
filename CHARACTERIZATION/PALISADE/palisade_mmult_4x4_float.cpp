#include <iostream>
#include <chrono>

#include <iomanip>

#include "palisade.h"
#include "../functional_units/functional_units.hpp"

using namespace lbcrypto;
using namespace std;

int main(void) {
  uint32_t multDepth = 1;
  uint32_t scaleFactorBits = 50;
  uint32_t slots = 8;
  SecurityLevel securityLevel = HEStd_128_classic;
  uint32_t ringDimension = 0;
  CryptoContext<DCRTPoly> cc = CryptoContextFactory<
    DCRTPoly>::genCryptoContextCKKS(multDepth,
      scaleFactorBits, slots, securityLevel, ringDimension, EXACTRESCALE);
  cc->Enable(ENCRYPTION);
  cc->Enable(SHE);
  cc->Enable(LEVELEDSHE);
  auto keyPair = cc->KeyGen();
  cc->EvalMultKeyGen(keyPair.secretKey);
  vector<complex<double>> tmp_vec_(slots);
  Plaintext tmp;
  Ciphertext<DCRTPoly> tmp_;

  vector<Ciphertext<DCRTPoly>> m1_, m2_, res_;
  Ciphertext<DCRTPoly> t3mp_;
  int r1, r2, c2, i, j, k, i_r1, k_c2, i_c1;
  int m;
  r1 = 4;
  r2 = 4;
  c2 = 4;
  m = r1*r1;
  m1_.resize(m);
  m2_.resize(m);
  res_.resize(m);
  for (i = 0; i<m; i++) {
    vector<double> tmp_vec_1(slots, 1.0+i);
    tmp = cc->MakeCKKSPackedPlaintext(tmp_vec_1);
    m1_[i] = cc->Encrypt(keyPair.publicKey, tmp);
;
    vector<double> tmp_vec_2(slots, 1.0+i);
    tmp = cc->MakeCKKSPackedPlaintext(tmp_vec_2);
    m2_[i] = cc->Encrypt(keyPair.publicKey, tmp);
    vector<double> tmp_vec_3(slots, 0.0);
    tmp = cc->MakeCKKSPackedPlaintext(tmp_vec_3);
    res_[i] = cc->Encrypt(keyPair.publicKey, tmp);
  }
  auto start_timer = chrono::high_resolution_clock::now();
  for (i = 0; i<r1; i++) {
    i_r1 = i*r1;
    i_c1 = i*r2;
    for (j = 0; j<c2; j++) {
      for (k = 0; k<r2; k++) {
        k_c2 = k*c2;
        Ciphertext<DCRTPoly> tmp_4_;
        tmp_4_ = cc->EvalMultAndRelinearize((m1_[(i_c1+k)]), (m2_[(k_c2+j)]));
        t3mp_ = (tmp_4_);
        cc->Rescale(t3mp_);
        res_[(i_r1+j)] = cc->EvalAdd(res_[(i_r1+j)], t3mp_);
      }
    }
  }
  auto stop_timer = chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<chrono::milliseconds>(stop_timer-start_timer);
  cout << "Time: " << duration.count() << "ms" << endl;
  for (i = 0; i<(r1*c2); i++) {
    cc->Decrypt(keyPair.secretKey,res_[i], &tmp);
    tmp->SetLength(1);
    tmp_vec_ = tmp->GetCKKSPackedValue();
    cout << "dec(res_[i]) = " << fixed << setprecision(1) << real(tmp_vec_[0]) << endl;
  }
  return 0;
}