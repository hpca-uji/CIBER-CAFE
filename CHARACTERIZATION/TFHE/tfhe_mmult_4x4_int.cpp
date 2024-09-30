#include <iostream>
#include <bitset>
#include <chrono>

#include <tfhe/tfhe.h>
#include <tfhe/tfhe_io.h>
#include <tfhe/tfhe_generic_streams.h>

#include "functional_units/functional_units.hpp"

using namespace std;

int main(void) {
  const size_t word_sz = 8;
  const size_t minimum_lambda = 128;
  TFheGateBootstrappingParameterSet* params =
    new_default_gate_bootstrapping_parameters(minimum_lambda);
  uint32_t seed[] = { 314, 1592, 657 };
  tfhe_random_generator_setSeed(seed, 3);
  TFheGateBootstrappingSecretKeySet* key =
    new_random_gate_bootstrapping_secret_keyset(params);
  LweSample *tmp, *tmp_;

  vector<vector<LweSample*>> m1_, m2_, res_;
  int r1, r2, c2, i, j, k, i_r1, k_c2, i_c1;
  int dim;
  r1 = 4;
  r2 = 4;
  c2 = 4;
  dim = r1*r2;
  m1_.resize(dim);
  m2_.resize(dim);
  res_.resize(dim);
  for (i = 0; i<dim; i++) {
    m1_[i] = e_cloud(1+i, word_sz, &key->cloud);
    m2_[i] = e_cloud(1+i, word_sz, &key->cloud);
    res_[i] = e_cloud(0, word_sz, &key->cloud);
  }
  auto start_timer = chrono::high_resolution_clock::now();
  for (i = 0; i<r1; i++) {
    i_r1 = i*r1;
    i_c1 = i*r2;
    for (j = 0; j<c2; j++) {
      for (k = 0; k<r2; k++) {
        k_c2 = k*c2;
        vector<LweSample*> tmp_1_(8);
        tmp_1_ = e_cloud(0, word_sz, &key->cloud);
        mult(tmp_1_, (m1_[(i_c1+k)]), (m2_[(k_c2+j)]), word_sz, &key->cloud);
        add(res_[(i_r1+j)], res_[(i_r1+j)], (tmp_1_), word_sz, &key->cloud);
      }
    }
  }
  auto stop_timer = chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<chrono::milliseconds>(stop_timer-start_timer);
  cout << "Time: " << duration.count() << "ms" << endl;
  for (i = 0; i<(r1*c2); i++) {
    vector<uint32_t> tmp_vec_2 = d_client(word_sz, res_[i], key);
    for (auto v : tmp_vec_2) {
      cout << "dec(res_[i]) = " << v << " ";
    }
    cout << endl;
  }
  return 0;
}