#include <iostream>
#include <chrono>

#include "seal/seal.h"
#include "../functional_units/functional_units.hpp"

using namespace seal;
using namespace std;

int main(void) {
  size_t word_sz = 0;
  size_t poly_modulus_degree = 4096;
  EncryptionParameters parms(scheme_type::bfv);
  parms.set_poly_modulus_degree(poly_modulus_degree);
  parms.set_coeff_modulus(CoeffModulus::BFVDefault(poly_modulus_degree));
  parms.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, 20));
  SEALContext context(parms);
  KeyGenerator keygen(context);
  SecretKey secret_key = keygen.secret_key();
  PublicKey public_key;
  RelinKeys relin_keys;
  keygen.create_public_key(public_key);
  keygen.create_relin_keys(relin_keys);
  Encryptor encryptor(context, public_key);
  Evaluator evaluator(context);
  Decryptor decryptor(context, secret_key);
  BatchEncoder batch_encoder(context);
  size_t slots = poly_modulus_degree/2;
  Plaintext tmp;
  Ciphertext tmp_;

  vector<Ciphertext> m1_, m2_, res_;
  int r1, r2, c2, i, j, k, i_r1, k_c2, i_c1;
  int m;
  r1 = 8;
  r2 = 8;
  c2 = 8;
  m = r1*r2;
  m1_.resize(m);
  m2_.resize(m);
  res_.resize(m);
  for (i = 0; i<m; i++) {
    tmp = uint64_to_hex_string(1+i);
    encryptor.encrypt(tmp, m1_[i]);
    tmp = uint64_to_hex_string(1+i);
    encryptor.encrypt(tmp, m2_[i]);
    tmp = uint64_to_hex_string(0);
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
        evaluator.add(res_[(i_r1+j)], (tmp_1_), res_[(i_r1+j)]);
      }
    }
  }
  auto stop_timer = chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<chrono::milliseconds>(stop_timer-start_timer);
  cout << "Time: " << duration.count() << "ms" << endl;
  for (i = 0; i<(r1*c2); i++) {
    decryptor.decrypt(res_[i], tmp);
    cout << "dec(res_[i]) = " << tmp << endl;
  }
  return 0;
}