#include <iostream>
#include <chrono>

#include "palisade.h"
#include "../functional_units/functional_units.hpp"

using namespace lbcrypto;
using namespace std;

int main(void) {
  uint32_t depth = 1;
  double sigma = 3.2;
  SecurityLevel securityLevel = HEStd_128_classic;
  size_t plaintext_modulus = 40961;
  CryptoContext<DCRTPoly> cc = CryptoContextFactory<
    DCRTPoly>::genCryptoContextBFVrns(plaintext_modulus,
      securityLevel, sigma, 0, depth, 0, OPTIMIZED, depth, 0, 60, 4096);
  cc->Enable(ENCRYPTION);
  cc->Enable(SHE);
  auto keyPair = cc->KeyGen();
  cc->EvalMultKeyGen(keyPair.secretKey);
  size_t slots(cc->GetRingDimension());
  vector<int64_t> tmp_vec_(slots);
  Plaintext tmp;
  Ciphertext<DCRTPoly> tmp_;

  vector<Ciphertext<DCRTPoly>> m1_, m2_, res_;
  int r1, r2, c2, i, j, k, i_r1, k_c2, i_c1;
  int m;
  r1 = 4;
  r2 = 4;
  c2 = 4;
  m = r1*r2;
  m1_.resize(m);
  m2_.resize(m);
  res_.resize(m);
  for (i = 0; i<m; i++) {
    fill(tmp_vec_.begin(), tmp_vec_.end(), 1+i);
    tmp = cc->MakePackedPlaintext(tmp_vec_);
    m1_[i] = cc->Encrypt(keyPair.publicKey, tmp);
    fill(tmp_vec_.begin(), tmp_vec_.end(), 1+i);
    tmp = cc->MakePackedPlaintext(tmp_vec_);
    m2_[i] = cc->Encrypt(keyPair.publicKey, tmp);
    fill(tmp_vec_.begin(), tmp_vec_.end(), 0);
    tmp = cc->MakePackedPlaintext(tmp_vec_);
    res_[i] = cc->Encrypt(keyPair.publicKey, tmp);
  }
  auto start_timer = chrono::high_resolution_clock::now();
  for (i = 0; i<r1; i++) {
    i_r1 = i*r1;
    i_c1 = i*r2;
    for (j = 0; j<c2; j++) {
      for (k = 0; k<r2; k++) {
        k_c2 = k*c2;
        Ciphertext<DCRTPoly> tmp_1_;
        tmp_1_ = cc->EvalMultAndRelinearize((m1_[(i_c1+k)]), (m2_[(k_c2+j)]));
        res_[(i_r1+j)] = cc->EvalAdd(res_[(i_r1+j)], (tmp_1_));
      }
    }
  }
  auto stop_timer = chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<chrono::milliseconds>(stop_timer-start_timer);
  cout << "Time: " << duration.count() << "ms" << endl;
  for (i = 0; i<(r1*c2); i++) {
    cc->Decrypt(keyPair.secretKey,res_[i], &tmp);
    tmp->SetLength(1);
    tmp_vec_ = tmp->GetPackedValue();
    cout << "dec(res_[i]) = " << tmp_vec_[0] << endl;
  }
  return 0;
}