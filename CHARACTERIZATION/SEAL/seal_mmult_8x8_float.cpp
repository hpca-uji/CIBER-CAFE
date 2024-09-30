#include <iostream>
#include <chrono>

#include <iomanip>

#include "seal/seal.h"
#include "../functional_units/functional_units.hpp"

using namespace seal;
using namespace std;

int main(void) {
  size_t word_sz = 0;
  EncryptionParameters parms(scheme_type::ckks);
  size_t poly_modulus_degree = 8192;
  parms.set_poly_modulus_degree(poly_modulus_degree);
  parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { 60, 40, 60 }));
  double scale = pow(2.0, 40);
  SEALContext context(parms);
  KeyGenerator keygen(context);
  SecretKey secret_key = keygen.secret_key();
  PublicKey public_key;
  RelinKeys relin_keys;
  GaloisKeys gal_keys;
  keygen.create_public_key(public_key);
  keygen.create_relin_keys(relin_keys);
  keygen.create_galois_keys(gal_keys);
  Encryptor encryptor(context, public_key);
  Evaluator evaluator(context);
  Decryptor decryptor(context, secret_key);
  CKKSEncoder encoder(context);
  size_t slot_count = encoder.slot_count();
  Plaintext tmp;
  Ciphertext tmp_;

  vector<Ciphertext> m1_, m2_, res_;
  Ciphertext t3mp_;
  int r1, r2, c2, i, j, k, i_r1, k_c2, i_c1;
  int m;
  r1 = 8;
  r2 = 8;
  c2 = 8;
  m = r1*r1;
  m1_.resize(m);
  m2_.resize(m);
  res_.resize(m);
  for (i = 0; i<m; i++) {
    encoder.encode(1.0+i, scale, tmp);
    encryptor.encrypt(tmp, m1_[i]);
;
    encoder.encode(1.0+i, scale, tmp);
    encryptor.encrypt(tmp, m2_[i]);
    encoder.encode(0.0, scale, tmp);
    encryptor.encrypt(tmp, res_[i]);
  }
  auto start_timer = chrono::high_resolution_clock::now();
  for (i = 0; i<r1; i++) {
    i_r1 = i*r1;
    i_c1 = i*r2;
    for (j = 0; j<c2; j++) {
      for (k = 0; k<r2; k++) {
        k_c2 = k*c2;
        Ciphertext tmp_1_;
        evaluator.multiply((m1_[(i_c1+k)]), (m2_[(k_c2+j)]), tmp_1_);
        evaluator.relinearize_inplace(tmp_1_, relin_keys);
        t3mp_ = (tmp_1_);
        evaluator.rescale_to_next_inplace(t3mp_);
        evaluator.mod_switch_to_inplace(res_[(i_r1+j)], t3mp_.parms_id());
        res_[(i_r1+j)].scale() = scale;
        t3mp_.scale() = scale;
        evaluator.add(res_[(i_r1+j)], t3mp_, res_[(i_r1+j)]);
      }
    }
  }
  auto stop_timer = chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<chrono::milliseconds>(stop_timer-start_timer);
  cout << "Time: " << duration.count() << "ms" << endl;
  for (i = 0; i<(r1*c2); i++) {
    decryptor.decrypt(res_[i], tmp);
    vector<double> tmp_vec_2;
    encoder.decode(tmp, tmp_vec_2);
    cout << "dec(res_[i]) = " << fixed << setprecision(1) << tmp_vec_2[0] << endl;
  }
  return 0;
}