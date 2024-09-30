#################################################
# https://eprint.iacr.org/2018/1041.pdf
# https://github.com/mark-schultz/EVA-matrix-multiplication

import math
import numpy

def sum_list(vect1, vect2):
    ret = []
    for e1, e2 in zip(vect1, vect2):
        ret.append(e1 + e2)
    return ret

def mult_list(vect1, vect2):
    ret = []
    for e1, e2 in zip(vect1, vect2):
        ret.append(e1 * e2)
    return ret

def rotar_list(vect, k):
    return vect[k:]+vect[:k]

def vec_from_pred(n, pred):
    # Returns a vector v in {0,1}^n s.t.
    # v[i] = pred(i)
    return [1 if pred(ell) else 0 for ell in range(n)]

class HEMatMult:
    def __init__(self, tam_vect, scale, min_scale, coeff_modulus_array, evaluator, encoder, relin_keys, gal_keys, decryptor = None):
        self.tam_vect = tam_vect
        self.scale = scale
        self.min_scale = min_scale
        self.coeff_modulus_array = coeff_modulus_array
        # La escala maxima para cada coeff_modulus_size
        self.max_scales = {i: 2**(sum(coeff_modulus_array[:i])-2) for i in range(1, len(coeff_modulus_array))}
        self.evaluator = evaluator
        self.encoder = encoder
        self.relin_keys = relin_keys
        self.gal_keys = gal_keys
        self.decryptor = decryptor

    def reescalar(self, crypt):
        if crypt.scale() / (2**self.coeff_modulus_array[crypt.coeff_modulus_size()-1]) < self.min_scale:
            cod = self.encoder.encode(1.0, self.min_scale * (2**self.coeff_modulus_array[crypt.coeff_modulus_size()-1]) / crypt.scale())
            if crypt.parms_id() != cod.parms_id():
                self.evaluator.mod_switch_to_inplace(cod, crypt.parms_id())
            crypt = self.evaluator.multiply_plain(crypt, cod)
        self.evaluator.rescale_to_next_inplace(crypt)
        crypt.scale(2**int(math.log2(crypt.scale())))
        return crypt
    
    def multiply_plain(self, crypt, cod):
        # Si el crypt no tiene el coeff_modulus_size original cambia el cod para adaptarse a este
        if crypt.parms_id() != cod.parms_id():
            self.evaluator.mod_switch_to_inplace(cod, crypt.parms_id())

        # Comprueba que la multiplicación no exceda la escala maxima y si lo hace realiza un rescalado
        while crypt.scale()*cod.scale() > self.max_scales[crypt.coeff_modulus_size()]:
            crypt = self.reescalar(crypt)
            self.evaluator.mod_switch_to_next_inplace(cod)

        ret = self.evaluator.multiply_plain(crypt, cod)
        self.evaluator.relinearize_inplace(ret, self.relin_keys)
        return ret
    
    def multiply_vect(self, crypt, vect):
        # Si es una lista la multiplica para que ocupe todo el tamaño maximo
        # Sino se asegura que el numero sea un decimal
        return self.multiply_plain(crypt, self.encoder.encode((vect * int(self.tam_vect/len(vect))) if type(vect) == list else float(vect), self.scale))
    
    def multiply(self, crypt1, crypt2):
        # Comprueba que los dos vectores codificados tengan el mismo coeff_modulus_size
        # y sino lo tienen los adapta en consecuencia
        while crypt1.coeff_modulus_size() > crypt2.coeff_modulus_size():
            if crypt1.scale() != self.min_scale:
                crypt1 = self.reescalar(crypt1)
            else:
                self.evaluator.mod_switch_to_next_inplace(crypt1)

        while crypt1.coeff_modulus_size() < crypt2.coeff_modulus_size():
            if crypt2.scale() != self.min_scale:
                crypt2 = self.reescalar(crypt2)
            else:
                self.evaluator.mod_switch_to_next_inplace(crypt2)

        # Comprueba que la multiplicación no exceda la escala maxima y si lo hace realiza un rescalado
        while crypt1.scale()*crypt2.scale() > self.max_scales[crypt1.coeff_modulus_size()]:
            if crypt1.scale() > crypt2.scale():
                crypt1 = self.reescalar(crypt1)
                if crypt2.scale() != self.min_scale:
                    crypt2 = self.reescalar(crypt2)
                else:
                    self.evaluator.mod_switch_to_next_inplace(crypt2)
            elif crypt1.scale() < crypt2.scale():
                crypt2 = self.reescalar(crypt2)
                if crypt1.scale() != self.min_scale:
                    crypt1 = self.reescalar(crypt1)
                else:
                    self.evaluator.mod_switch_to_next_inplace(crypt1)
            elif crypt1.scale() == crypt2.scale() != self.min_scale:
                crypt1 = self.reescalar(crypt1)
                crypt2 = self.reescalar(crypt2)
            else:
                raise Exception(f"No se puede continuar el calculo con la escala minima establecida.")
            
        ret = self.evaluator.multiply(crypt1, crypt2)
        self.evaluator.relinearize_inplace(ret, self.relin_keys)
        return ret

    def add(self, crypt1, crypt2):
        # Se asegura que ambos esten en el mismo coeff_modulus_size
        while crypt1.coeff_modulus_size() > crypt2.coeff_modulus_size():
            if crypt1.scale() != self.min_scale:
                crypt1 = self.reescalar(crypt1)
            else:
                self.evaluator.mod_switch_to_next_inplace(crypt1)

        while crypt1.coeff_modulus_size() < crypt2.coeff_modulus_size():
            if crypt2.scale() != self.min_scale:
                crypt2 = self.reescalar(crypt2)
            else:
                self.evaluator.mod_switch_to_next_inplace(crypt2)
        
        # Se asegura que ambos tengan la misma escala
        if crypt1.scale() != crypt2.scale():
            if crypt1.scale() > crypt2.scale():
                crypt2 = self.multiply_plain(crypt2, self.encoder.encode(1.0, crypt1.scale() / crypt2.scale()))
            else:
                crypt1 = self.multiply_plain(crypt1, self.encoder.encode(1.0, crypt2.scale() / crypt1.scale()))
        
        ret = self.evaluator.add(crypt1, crypt2)
        return ret

    def multiplica(self, a, b, d):
        n = d**2
        if int(self.tam_vect/n) != self.tam_vect/n:
            raise Exception(f"El lado (d: {d}) no es valido, debe ser una potencia de 2 cuyo cuadrado quepa en el vector codificado (d²: {d**2} < {self.tam_vect}")
        k = -d+1

        # Step 1-1 - ctA0
        #print("Step 1-1 - ctA0")
        if k >= 0:
            uk = vec_from_pred(n, lambda ell: 0 <= ell - d * k < (d-k))
        else:
            uk = vec_from_pred(n, lambda ell: -k <= ell - (d+k) * d < d)
        ctA0 = self.multiply_vect(self.evaluator.rotate_vector(a, k, self.gal_keys), uk)
        for k in range(-d+2, d):
            if k >= 0:
                uk = vec_from_pred(n, lambda ell: 0 <= ell - d * k < (d-k))
            else:
                uk = vec_from_pred(n, lambda ell: -k <= ell - (d+k) * d < d)
            ctA0 = self.add(ctA0, self.multiply_vect(self.evaluator.rotate_vector(a, k, self.gal_keys) if k != 0 else a, uk))
        
        # Step 1-2 - ctB0
        #print("Step 1-2 - ctB0")
        k=0
        ctB0 = self.multiply_vect(b, vec_from_pred(n, lambda ell: (ell % d) == k))
        for k in range(1, d):
            ctB0 = self.add(ctB0, self.multiply_vect(self.evaluator.rotate_vector(b, d * k, self.gal_keys), vec_from_pred(n, lambda ell: ell % d == k)))

        # Step 2&3 - ctAB
        #print("Step 2&3 - ctAB")
        ctAB = self.multiply(ctA0, ctB0)
        for k in range(1, d):
            vk = vec_from_pred(n, lambda ell: 0 <= ell % d < d - k)
            vk_minus_d = vec_from_pred(n, lambda ell: (d-k) <= ell % d < d)
            ctAk = self.evaluator.add(self.multiply_vect(self.evaluator.rotate_vector(ctA0, k, self.gal_keys), vk),
                                      self.multiply_vect(self.evaluator.rotate_vector(ctA0, k-d, self.gal_keys), vk_minus_d))
            ctBk = self.evaluator.rotate_vector(ctB0, d * k, self.gal_keys)
            ctAB = self.add(ctAB, self.multiply(ctAk, ctBk))

        return ctAB
        
##################################################################################################
    def reescalar_coment(self, crypt, vectDecrypt):
        if crypt.scale() / (2**self.coeff_modulus_array[crypt.coeff_modulus_size()-1]) < self.min_scale:
            print(f"\t\t\tReescalando de {math.log2(crypt.scale())} a {math.log2(self.min_scale)} el rescalado restara {self.coeff_modulus_array[crypt.coeff_modulus_size()-1]}")
            cod = self.encoder.encode(1.0, self.min_scale * (2**self.coeff_modulus_array[crypt.coeff_modulus_size()-1]) / crypt.scale())
            if crypt.parms_id() != cod.parms_id():
                self.evaluator.mod_switch_to_inplace(cod, crypt.parms_id())
                print(f"\t\tcod - mod_switch_to: {math.log2(cod.scale())} | {cod.coeff_count()} | {cod.parms_id()}")
                print(f"\t\t\terror:{self.comparaDec(cod, [1]*self.tam_vect)}")
            crypt = self.evaluator.multiply_plain(crypt, cod)
        self.evaluator.rescale_to_next_inplace(crypt)
        crypt.scale(2**int(math.log2(crypt.scale())))
        return crypt
    
    def multiply_plain_coment(self, crypt, cod, vectDecrypt, vectDecod):
        print("\tmultiply_plain:")
        print(f"\t\tcrypt: {math.log2(crypt.scale())} | {crypt.size()}/{crypt.size_capacity()} | {crypt.coeff_modulus_size()} | {crypt.parms_id()}")
        print(f"\t\tcod: {math.log2(cod.scale())} | {cod.coeff_count()} | {cod.parms_id()}")
        # Si el crypt no tiene el coeff_modulus_size original cambia el cod para adaptarse a este
        if crypt.parms_id() != cod.parms_id():
            self.evaluator.mod_switch_to_inplace(cod, crypt.parms_id())
            print(f"\t\tcod - mod_switch_to: {math.log2(cod.scale())} | {cod.coeff_count()} | {cod.parms_id()}")
            print(f"\t\t\terror:{self.comparaDec(cod, vectDecod)}")
        # Comprueba que la multiplicación no exceda la escala maxima y si lo hace realiza un rescalado
        while crypt.scale()*cod.scale() > self.max_scales[crypt.coeff_modulus_size()]:
            crypt = self.reescalar_coment(crypt, vectDecrypt)
            print(f"\t\tcrypt - rescale_to: {math.log2(crypt.scale())} | {crypt.size()}/{crypt.size_capacity()} | {crypt.coeff_modulus_size()} | {crypt.parms_id()}")
            print(f"\t\t\terror:{self.compara(crypt, vectDecrypt)}")
            self.evaluator.mod_switch_to_next_inplace(cod)
            print(f"\t\tcod - mod_switch_to: {math.log2(cod.scale())} | {cod.coeff_count()} | {cod.parms_id()}")
            print(f"\t\t\terror:{self.comparaDec(cod, vectDecod)}")

        ret = self.evaluator.multiply_plain(crypt, cod)
        print(f"\t\tmultiply_plain: {math.log2(ret.scale())} | {ret.size()}/{ret.size_capacity()} | {ret.coeff_modulus_size()} | {ret.parms_id()}")
        print(f"\t\t\terror:{self.compara(ret, mult_list(vectDecrypt, vectDecod))}")
        self.evaluator.relinearize_inplace(ret, self.relin_keys)
        print(f"\t\trelinearize: {math.log2(ret.scale())} | {ret.size()}/{ret.size_capacity()} | {ret.coeff_modulus_size()} | {ret.parms_id()}")
        print(f"\t\t\terror:{self.compara(ret, mult_list(vectDecrypt, vectDecod))}")
        return ret
    
    def multiply_vect_coment(self, crypt, vect, vectDecrypt):
        # Si es una lista la multiplica para que ocupe todo el tamaño maximo
        # Sino se asegura que el numero sea un decimal
        vect = (vect * int(self.tam_vect/len(vect))) if type(vect) == list else float(vect)
        return self.multiply_plain_coment(crypt, self.encoder.encode(vect, self.scale), vect, vectDecrypt)
    
    def multiply_coment(self, crypt1, crypt2, vectDecrypt1, vectDecrypt2):
        print("\tmultiply:")
        print(f"\t\tcrypt1: {math.log2(crypt1.scale())} | {crypt1.size()}/{crypt1.size_capacity()} | {crypt1.coeff_modulus_size()} | {crypt1.parms_id()}")
        print(f"\t\tcrypt2: {math.log2(crypt2.scale())} | {crypt2.size()}/{crypt2.size_capacity()} | {crypt2.coeff_modulus_size()} | {crypt2.parms_id()}")

        # Comprueba que los dos vectores codificados tengan el mismo coeff_modulus_size
        # y sino lo tienen los adapta en consecuencia
        while crypt1.coeff_modulus_size() > crypt2.coeff_modulus_size():
            if crypt1.scale() != self.min_scale:
                crypt1 = self.reescalar_coment(crypt1, vectDecrypt1)
                print(f"\t\tcrypt1 - rescale_to: {math.log2(crypt1.scale())} | {crypt1.size()}/{crypt1.size_capacity()} | {crypt1.coeff_modulus_size()} | {crypt1.parms_id()}")
            else:
                self.evaluator.mod_switch_to_next_inplace(crypt1)
                print(f"\t\tcrypt1 - mod_switch_to: {math.log2(crypt1.scale())} | {crypt1.size()}/{crypt1.size_capacity()} | {crypt1.coeff_modulus_size()} | {crypt1.parms_id()}")
            
        while crypt1.coeff_modulus_size() < crypt2.coeff_modulus_size():
            if crypt2.scale() != self.min_scale:
                crypt2 = self.reescalar_coment(crypt2, vectDecrypt2)
                print(f"\t\tcrypt2 - rescale_to: {math.log2(crypt2.scale())} | {crypt2.size()}/{crypt2.size_capacity()} | {crypt2.coeff_modulus_size()} | {crypt2.parms_id()}")
            else:
                self.evaluator.mod_switch_to_next_inplace(crypt2)
                print(f"\t\tcrypt2 - mod_switch_to: {math.log2(crypt2.scale())} | {crypt2.size()}/{crypt2.size_capacity()} | {crypt2.coeff_modulus_size()} | {crypt2.parms_id()}")
        
        # Comprueba que la multiplicación no exceda la escala maxima y si lo hace realiza un rescalado
        while crypt1.scale()*crypt2.scale() > self.max_scales[crypt1.coeff_modulus_size()]:
            if crypt1.scale() > crypt2.scale():
                crypt1 = self.reescalar_coment(crypt1, vectDecrypt1)
                print(f"\t\tcrypt1 - rescale_to: {math.log2(crypt1.scale())} | {crypt1.size()}/{crypt1.size_capacity()} | {crypt1.coeff_modulus_size()} | {crypt1.parms_id()}")
                if crypt2.scale() != self.min_scale:
                    crypt2 = self.reescalar_coment(crypt2, vectDecrypt2)
                    print(f"\t\tcrypt2 - rescale_to: {math.log2(crypt2.scale())} | {crypt2.size()}/{crypt2.size_capacity()} | {crypt2.coeff_modulus_size()} | {crypt2.parms_id()}")
                else:
                    self.evaluator.mod_switch_to_next_inplace(crypt2)
                    print(f"\t\tcrypt2 - mod_switch_to: {math.log2(crypt2.scale())} | {crypt2.size()}/{crypt2.size_capacity()} | {crypt2.coeff_modulus_size()} | {crypt2.parms_id()}")
            elif crypt1.scale() < crypt2.scale():
                if crypt1.scale() != self.min_scale:
                    crypt1 = self.reescalar_coment(crypt1, vectDecrypt1)
                    print(f"\t\tcrypt1 - rescale_to: {math.log2(crypt1.scale())} | {crypt1.size()}/{crypt1.size_capacity()} | {crypt1.coeff_modulus_size()} | {crypt1.parms_id()}")
                else:
                    self.evaluator.mod_switch_to_next_inplace(crypt1)
                    print(f"\t\tcrypt1 - mod_switch_to: {math.log2(crypt1.scale())} | {crypt1.size()}/{crypt1.size_capacity()} | {crypt1.coeff_modulus_size()} | {crypt1.parms_id()}")
                crypt2 = self.reescalar_coment(crypt2, vectDecrypt2)
                print(f"\t\tcrypt2 - rescale_to: {math.log2(crypt2.scale())} | {crypt2.size()}/{crypt2.size_capacity()} | {crypt2.coeff_modulus_size()} | {crypt2.parms_id()}")
            elif crypt1.scale() == crypt2.scale() != self.min_scale:
                crypt1 = self.reescalar_coment(crypt1, vectDecrypt1)
                print(f"\t\tcrypt1 - rescale_to: {math.log2(crypt1.scale())} | {crypt1.size()}/{crypt1.size_capacity()} | {crypt1.coeff_modulus_size()} | {crypt1.parms_id()}")
                crypt2 = self.reescalar_coment(crypt2, vectDecrypt2)
                print(f"\t\tcrypt2 - rescale_to: {math.log2(crypt2.scale())} | {crypt2.size()}/{crypt2.size_capacity()} | {crypt2.coeff_modulus_size()} | {crypt2.parms_id()}")
            else:
                raise Exception(f"No se puede continuar el calculo con la escala minima establecida.")
            
        ret = self.evaluator.multiply(crypt1, crypt2)
        print(f"\t\tmultiply: {math.log2(ret.scale())} | {ret.size()}/{ret.size_capacity()} | {ret.coeff_modulus_size()} | {ret.parms_id()}")
        print(f"\t\t\terror:{self.compara(ret, mult_list(vectDecrypt1, vectDecrypt2))}")
        self.evaluator.relinearize_inplace(ret, self.relin_keys)
        print(f"\t\tmultiply - relinearize: {math.log2(ret.scale())} | {ret.size()}/{ret.size_capacity()} | {ret.coeff_modulus_size()} | {ret.parms_id()}")
        print(f"\t\t\terror:{self.compara(ret, mult_list(vectDecrypt1, vectDecrypt2))}")
        
        return ret

    def add_coment(self, crypt1, crypt2, vectDecrypt1, vectDecrypt2):
        print("\tadd:")
        print(f"\t\tcrypt1: {math.log2(crypt1.scale())} | {crypt1.size()}/{crypt1.size_capacity()} | {crypt1.coeff_modulus_size()} | {crypt1.parms_id()}")
        print(f"\t\tcrypt2: {math.log2(crypt2.scale())} | {crypt2.size()}/{crypt2.size_capacity()} | {crypt2.coeff_modulus_size()} | {crypt2.parms_id()}")
        # Se asegura que ambos esten en el mismo coeff_modulus_size
        while crypt1.coeff_modulus_size() > crypt2.coeff_modulus_size():
            if crypt1.scale() != self.min_scale:
                if crypt1.scale() / (2**self.coeff_modulus_array[crypt1.coeff_modulus_size()-1]) < self.min_scale:
                    crypt1 = self.multiply_plain_coment(crypt1, self.encoder.encode(1.0, crypt1.scale() / self.min_scale / (2**self.coeff_modulus_array[crypt1.coeff_modulus_size()-1])), vectDecrypt1, [1]*self.tam_vect)
                self.evaluator.rescale_to_next_inplace(crypt1)
                crypt1.scale(2**int(math.log2(crypt1.scale())))
                print(f"\t\tcrypt1 - rescale_to: {math.log2(crypt1.scale())} | {crypt1.size()}/{crypt1.size_capacity()} | {crypt1.coeff_modulus_size()} | {crypt1.parms_id()}")
            else:
                self.evaluator.mod_switch_to_next_inplace(crypt1)
                print(f"\t\tcrypt1 - mod_switch_to: {math.log2(crypt1.scale())} | {crypt1.size()}/{crypt1.size_capacity()} | {crypt1.coeff_modulus_size()} | {crypt1.parms_id()}")
        while crypt1.coeff_modulus_size() < crypt2.coeff_modulus_size():
            if crypt2.scale() - 2**self.coeff_modulus_array[crypt2.coeff_modulus_size()-1] >= self.min_scale:
                self.evaluator.rescale_to_next_inplace(crypt2)
                crypt2.scale(2**int(math.log2(crypt2.scale())))
                print(f"\t\tcrypt2 - rescale_to: {math.log2(crypt2.scale())} | {crypt2.size()}/{crypt2.size_capacity()} | {crypt2.coeff_modulus_size()} | {crypt2.parms_id()}")
            elif crypt2.scale() <= self.max_scales[crypt2.coeff_modulus_size()-1]:
                self.evaluator.mod_switch_to_next_inplace(crypt2)
                print(f"\t\tcrypt2 - mod_switch_to: {math.log2(crypt2.scale())} | {crypt2.size()}/{crypt2.size_capacity()} | {crypt2.coeff_modulus_size()} | {crypt2.parms_id()}")
            else:
                raise Exception("No se puede equiparar el 'coeff_modulus_size' sin descender de la escala deseada ni superar la escala restante")
        
        # Se asegura que ambos tengan la misma escala
        if crypt1.scale() != crypt2.scale():
            if crypt1.scale() > crypt2.scale():
                crypt2 = self.multiply_plain_coment(crypt2, self.encoder.encode(1.0, crypt1.scale() / crypt2.scale()), vectDecrypt2, [1]*self.tam_vect)
            else:
                crypt1 = self.multiply_plain_coment(crypt1, self.encoder.encode(1.0, crypt2.scale() / crypt1.scale()), vectDecrypt1, [1]*self.tam_vect)
        
        ret = self.evaluator.add(crypt1, crypt2)
        print(f"\t\tadd: {math.log2(ret.scale())} | {ret.size()}/{ret.size_capacity()} | {ret.coeff_modulus_size()} | {ret.parms_id()}")
        print(f"\t\t\terror:{self.compara(ret, sum_list(vectDecrypt1, vectDecrypt2))}")
        return ret

    def multiplica_coment(self, a, b, d, vA, vB):
        n = d**2
        vA *= int(self.tam_vect/len(vA))
        print(len(vA))
        vB *= int(self.tam_vect/len(vB))
        print(len(vB))
        if int(self.tam_vect/n) != self.tam_vect/n:
            raise Exception(f"El lado (d: {d}) no es valido, debe ser una potencia de 2 cuyo cuadrado quepa en el vector codificado (d²: {d**2} < {self.tam_vect}")
        k = -d+1

        # Step 1-1 - ctA0
        print("Step 1-1 - ctA0")
        if k >= 0:
            uk = vec_from_pred(n, lambda ell: 0 <= ell - d * k < (d-k))
        else:
            uk = vec_from_pred(n, lambda ell: -k <= ell - (d+k) * d < d)
        uk *= int(self.tam_vect/len(uk))
        ctA0 = self.multiply_vect_coment(self.evaluator.rotate_vector(a, k, self.gal_keys), uk, rotar_list(vA, k))
        vA0 = mult_list(rotar_list(vA, k), uk)

        print(f"ctA0[0]: {math.log2(ctA0.scale())} | {ctA0.size()}/{ctA0.size_capacity()} | {ctA0.coeff_modulus_size()} | {ctA0.parms_id()}")
        print(f"ctA0[0]: {self.compara(ctA0, vA0)}")
        for k in range(-d+2, d):
            if k >= 0:
                uk = vec_from_pred(n, lambda ell: 0 <= ell - d * k < (d-k))
            else:
                uk = vec_from_pred(n, lambda ell: -k <= ell - (d+k) * d < d)
            uk *= int(self.tam_vect/len(uk))
            ctA0 = self.add_coment(ctA0, self.multiply_vect_coment(self.evaluator.rotate_vector(a, k, self.gal_keys) if k != 0 else a, uk, rotar_list(vA, k)), vA0, mult_list(rotar_list(vA, k), uk))
            vA0 = sum_list(vA0, mult_list(rotar_list(vA, k), uk))
        print(f"ctA0: {math.log2(ctA0.scale())} | {ctA0.size()}/{ctA0.size_capacity()} | {ctA0.coeff_modulus_size()} | {ctA0.parms_id()}")
        print(f"ctA0: {self.compara(ctA0, vA0)}")
        
        # Step 1-2 - ctB0
        print("Step 1-2 - ctB0")
        k = 0
        uk = vec_from_pred(n, lambda ell: (ell % d) == k)
        uk *= int(self.tam_vect/len(uk))
        ctB0 = self.multiply_vect_coment(b, uk, vB)
        vB0 = mult_list(vB, uk)
        print(f"ctB0[0]: {math.log2(ctB0.scale())} | {ctB0.size()}/{ctB0.size_capacity()} | {ctB0.coeff_modulus_size()} | {ctB0.parms_id()}")
        print(f"ctB0[0]: {self.compara(ctB0, vB0)}")
        for k in range(1, d):
            uk = vec_from_pred(n, lambda ell: ell % d == k)
            uk *= int(self.tam_vect/len(uk))
            ctB0 = self.add_coment(ctB0, self.multiply_vect_coment(self.evaluator.rotate_vector(b, d * k, self.gal_keys), uk, rotar_list(vB, d*k)), vB0, mult_list(rotar_list(vB, d*k), uk))
            vB0 = sum_list(vB0, mult_list(rotar_list(vB, d*k), uk))
        print(f"ctB0: {math.log2(ctB0.scale())} | {ctB0.size()}/{ctB0.size_capacity()} | {ctB0.coeff_modulus_size()} | {ctB0.parms_id()}")
        print(f"ctB0: {self.compara(ctB0, vB0)}")
        # Step 2&3 - ctAB
        print("Step 2&3 - ctAB")
        ctAB = self.multiply_coment(ctA0, ctB0, vA0, vB0)
        vAB = mult_list(vA0, vB0)
        print(f"ctAB[0]: {math.log2(ctAB.scale())} | {ctAB.size()}/{ctAB.size_capacity()} | {ctAB.coeff_modulus_size()} | {ctAB.parms_id()}")
        print(f"ctAB[0]: {self.compara(ctAB, vAB)}")
        for k in range(1, d):
            vk = vec_from_pred(n, lambda ell: 0 <= ell % d < d - k)
            vk *= int(self.tam_vect/len(vk))
            vk_minus_d = vec_from_pred(n, lambda ell: (d-k) <= ell % d < d)
            vk_minus_d *= int(self.tam_vect/len(vk_minus_d))
            ctAk = self.add_coment(self.multiply_vect_coment(self.evaluator.rotate_vector(ctA0, k, self.gal_keys), vk, rotar_list(vA0, k)),
                                   self.multiply_vect_coment(self.evaluator.rotate_vector(ctA0, k-d, self.gal_keys), vk_minus_d, rotar_list(vA0, k-d)),
                                   mult_list(rotar_list(vA0, k), vk),
                                   mult_list(rotar_list(vA0, k-d), vk_minus_d))
            vAk = sum_list(mult_list(rotar_list(vA0, k), vk),
                           mult_list(rotar_list(vA0, k-d), vk_minus_d))
            ctBk = self.evaluator.rotate_vector(ctB0, d * k, self.gal_keys)
            vBk = rotar_list(vB0, d * k)
            ctAB = self.add_coment(ctAB, self.multiply_coment(ctAk, ctBk, vAk, vBk), vAB, mult_list(vAk, vBk))
            vAB = sum_list(vAB, mult_list(vAk, vBk))
        print(f"ctAB: {math.log2(ctAB.scale())} | {ctAB.size()}/{ctAB.size_capacity()} | {ctAB.coeff_modulus_size()} | {ctAB.parms_id()}")
        print(f"ctAB: {self.compara(ctAB, vAB)}")
        
        return ctAB
    
    def compara(self, ct, v):
        return self.comparaDec(self.decryptor.decrypt(ct), v)
    
    def comparaDec(self, cod, v):
        cod = self.encoder.decode(cod)
        t = abs(cod-v)
        return f"\t±{numpy.mean(t)}"
    
