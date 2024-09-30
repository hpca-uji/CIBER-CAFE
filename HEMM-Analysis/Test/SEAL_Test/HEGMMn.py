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

#   Input matrix A(m x l) y B(l x n)
#   Output matrix C(m x n) = A(m x l) x B(l x n)
class HEGMMEn:
    def __init__(self, fil_mA, col_mA_fil_mB, col_mB, tam_vect, scale, min_scale, coeff_modulus_array, evaluator, encoder, relin_keys, gal_keys, decryptor = None):
        self.m = fil_mA
        self.l = col_mA_fil_mB
        self.n = col_mB
        self.tam_vect = tam_vect
        self.chekTam()
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
    
    def chekTam(self):
        if max(self.m*self.l, self.l*self.n, self.m*self.n) > self.tam_vect:
            raise Exception(f"El tamaño del vector ({self.tam_vect}) es demasiado pequeño para almacenar alguna de las matrices: A({self.m}x{self.l})={self.m*self.l}, B({self.l}x{self.n})={self.l*self.n} o C({self.m}x{self.n})={self.m*self.n}")

    def setTamMatA(self, fil_mA, col_mA):
        self.m = fil_mA
        self.l = col_mA
        self.chekTam()

    def setTamMatB(self, fil_mB, col_mB):
        self.l = fil_mB
        self.n = col_mB
        self.chekTam()

    def getTamMatC(self):
        return self.m, self.n
    
    def p(self):
        return min(self.m, self.l, self.n)
    
    def t(self):
        return math.ceil(self.l/self.p())

    def M(self):
        return (self.m*self.t()) if self.p() == self.m else self.m

    def N(self):
        return (self.n*self.t()) if self.p() != self.m and self.p() == self.n else self.n
    
    def prepA(self, A):
        if self.p() == self.m:
            tA = A * self.t()
        else:
            tA = A
        oA = []
        for i in range(self.M()):
            oA += tA[i*self.l+(i%self.l):(i+1)*self.l] + tA[i*self.l:(i+1)*self.l-(self.l-i%self.l)]
            if self.N() > self.l:
                oA += oA[i*self.M():i*self.M()+self.N()-self.l]
        return oA
# ERROR: FALTA ARREGLAR EL CASO (m == n < l) - NO FUNCIONA
    def prepB(self, B):
        if self.p() != self.m and self.p() == self.n:
            tB = []
            for i in range(self.l):
                tB += B[i*self.n:(i+1)*self.n] * self.t()
        else:
            tB = B
        rB = []
        for i in range(self.l * self.N()):
            rB += [tB[(i+(i%self.N())*self.N())%(self.l*self.N())]]
        
        rB += rB[:(self.M()-self.l)*self.N()]
        return rB

    def arregloC(self, C):
        ret = []
        if self.N() != self.n:
            for i in range(self.M()):
                ti = C[i*self.N():(i+1)*self.N()]
                tj = ti[:self.n]
                for j in range(1, self.t()):
                    tj = sum_list(tj, ti[j*self.n:(j+1)*self.n])
                ret += tj
        elif self.M() != self.m:
            ret = C[:self.m*self.n]
            for i in range(1,self.t()):
                ret = sum_list(ret, C[i*self.m*self.n:(i+1)*self.m*self.n])
        else:
            ret = C
        return ret

##################################################################################################
    def multiply_plain(self, crypt, cod):
        # Si el crypt no tiene el coeff_modulus_size original cambia el cod para adaptarse a este
        if crypt.parms_id() != cod.parms_id():
            self.evaluator.mod_switch_to_inplace(cod, crypt.parms_id())

        # Comprueba que la multiplicación no exceda la escala maxima y si lo hace realiza un rescalado
        while crypt.scale()*cod.scale() > self.max_scales[crypt.coeff_modulus_size()]:
            if crypt.coeff_modulus_size() > 1 and crypt.scale() / (2**self.coeff_modulus_array[crypt.coeff_modulus_size()-1]) >= self.min_scale:
                self.evaluator.rescale_to_next_inplace(crypt)
                crypt.scale(2**int(math.log2(crypt.scale())))
                self.evaluator.mod_switch_to_next_inplace(cod)
            else:
                raise Exception("No se puede continuar la multiplicación sin descender de la escala mínima.")

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
            if crypt1.scale() - 2**self.coeff_modulus_array[crypt1.coeff_modulus_size()-1] >= self.min_scale:
                self.evaluator.rescale_to_next_inplace(crypt1)
                crypt1.scale(2**int(math.log2(crypt1.scale())))
            elif crypt1.scale() <= self.max_scales[crypt1.coeff_modulus_size()-1]:
                self.evaluator.mod_switch_to_next_inplace(crypt1)
            else:
                raise Exception("No se puede equiparar el 'coeff_modulus_size' sin descender de la escala deseada ni superar la escala restante")
        while crypt1.coeff_modulus_size() < crypt2.coeff_modulus_size():
            if crypt2.scale() - 2**self.coeff_modulus_array[crypt2.coeff_modulus_size()-1] >= self.min_scale:
                self.evaluator.rescale_to_next_inplace(crypt2)
                crypt2.scale(2**int(math.log2(crypt2.scale())))
            elif crypt2.scale() <= self.max_scales[crypt2.coeff_modulus_size()-1]:
                self.evaluator.mod_switch_to_next_inplace(crypt2)
            else:
                raise Exception("No se puede equiparar el 'coeff_modulus_size' sin descender de la escala deseada ni superar la escala restante")

        # Comprueba que la multiplicación no exceda la escala maxima y si lo hace realiza un rescalado
        while crypt1.scale()*crypt2.scale() > self.max_scales[crypt1.coeff_modulus_size()]:
            if crypt1.scale() > crypt2.scale():
                if crypt1.coeff_modulus_size() > 1 and crypt1.scale() / (2**self.coeff_modulus_array[crypt1.coeff_modulus_size()-1]) >= self.min_scale:
                    self.evaluator.rescale_to_next_inplace(crypt1)
                    crypt1.scale(2**int(math.log2(crypt1.scale())))
                else:
                    raise Exception(f"Precisión mínima necesaria ({math.ceil(math.log2(crypt1.scale() / (2**self.coeff_modulus_array[crypt1.coeff_modulus_size()-1])))})")
                if crypt2.coeff_modulus_size() > 1 and crypt2.scale() / (2**self.coeff_modulus_array[crypt2.coeff_modulus_size()-1]) >= self.min_scale:
                    self.evaluator.rescale_to_next_inplace(crypt2)
                    crypt2.scale(2**int(math.log2(crypt2.scale())))
                else:
                    self.evaluator.mod_switch_to_next_inplace(crypt2)
            elif crypt1.scale() < crypt2.scale():
                if crypt1.coeff_modulus_size() > 1 and crypt1.scale() / (2**self.coeff_modulus_array[crypt1.coeff_modulus_size()-1]) >= self.min_scale:
                    self.evaluator.rescale_to_next_inplace(crypt1)
                    crypt1.scale(2**int(math.log2(crypt1.scale())))
                else:
                    self.evaluator.mod_switch_to_next_inplace(crypt1)
                if crypt2.coeff_modulus_size() > 1 and crypt2.scale() / (2**self.coeff_modulus_array[crypt2.coeff_modulus_size()-1]) >= self.min_scale:
                    self.evaluator.rescale_to_next_inplace(crypt2)
                    crypt2.scale(2**int(math.log2(crypt2.scale())))
                else:
                    raise Exception(f"Precisión mínima necesaria ({math.ceil(math.log2(crypt2.scale() / (2**self.coeff_modulus_array[crypt2.coeff_modulus_size()-1])))})")
            else:
                if crypt1.coeff_modulus_size() > 1 and crypt1.scale() / (2**self.coeff_modulus_array[crypt1.coeff_modulus_size()-1]) >= self.min_scale:
                    self.evaluator.rescale_to_next_inplace(crypt1)
                    crypt1.scale(2**int(math.log2(crypt1.scale())))
                else:
                    raise Exception(f"Precisión mínima necesaria ({math.ceil(math.log2(crypt1.scale() / (2**self.coeff_modulus_array[crypt1.coeff_modulus_size()-1])))})")
                if crypt2.coeff_modulus_size() > 1 and crypt2.scale() / (2**self.coeff_modulus_array[crypt2.coeff_modulus_size()-1]) >= self.min_scale:
                    self.evaluator.rescale_to_next_inplace(crypt2)
                    crypt2.scale(2**int(math.log2(crypt2.scale())))
                else:
                    raise Exception(f"Precisión mínima necesaria ({math.ceil(math.log2(crypt2.scale() / (2**self.coeff_modulus_array[crypt2.coeff_modulus_size()-1])))})")
            
        ret = self.evaluator.multiply(crypt1, crypt2)
        self.evaluator.relinearize_inplace(ret, self.relin_keys)
        return ret

    def add(self, crypt1, crypt2):
        # Se asegura que ambos esten en el mismo coeff_modulus_size
        while crypt1.coeff_modulus_size() > crypt2.coeff_modulus_size():
            if crypt1.scale() - 2**self.coeff_modulus_array[crypt1.coeff_modulus_size()-1] >= self.min_scale:
                self.evaluator.rescale_to_next_inplace(crypt1)
                crypt1.scale(2**int(math.log2(crypt1.scale())))
            elif crypt1.scale() <= self.max_scales[crypt1.coeff_modulus_size()-1]:
                self.evaluator.mod_switch_to_next_inplace(crypt1)
            else:
                raise Exception("No se puede equiparar el 'coeff_modulus_size' sin descender de la escala deseada ni superar la escala restante")
        while crypt1.coeff_modulus_size() < crypt2.coeff_modulus_size():
            if crypt2.scale() - 2**self.coeff_modulus_array[crypt2.coeff_modulus_size()-1] >= self.min_scale:
                self.evaluator.rescale_to_next_inplace(crypt2)
                crypt2.scale(2**int(math.log2(crypt2.scale())))
            elif crypt2.scale() <= self.max_scales[crypt2.coeff_modulus_size()-1]:
                self.evaluator.mod_switch_to_next_inplace(crypt2)
            else:
                raise Exception("No se puede equiparar el 'coeff_modulus_size' sin descender de la escala deseada ni superar la escala restante")
        
        # Se asegura que ambos tengan la misma escala
        if crypt1.scale() != crypt2.scale():
            if crypt1.scale() > crypt2.scale():
                crypt2 = self.multiply_plain(crypt2, self.encoder.encode(1.0, crypt1.scale() / crypt2.scale()))
            else:
                crypt1 = self.multiply_plain(crypt1, self.encoder.encode(1.0, crypt2.scale() / crypt1.scale()))
        
        ret = self.evaluator.add(crypt1, crypt2)
        return ret

    def rotar_col(self, crypt, k):
        k %= self.N()
        if k != 0:
            eA1 = self.evaluator.rotate_vector(self.multiply_vect(crypt, [1 if i % self.N() < k + self.N() - self.l else 0 for i in range(self.M() * self.N())]), k - self.l, self.gal_keys)
            eA2 = self.evaluator.rotate_vector(self.multiply_vect(crypt, [1 if k <= i % self.N() < self.l else 0 for i in range(self.M() * self.N())]), k, self.gal_keys)
            return self.add(eA1, eA2)
        else:
            return crypt
# ERROR: FALTA ARREGLAR EL CASO (m == n < l) - NO FUNCIONA
    def rotar_fil(self, crypt, k):
        k %= self.M()
        if k != 0:
            wB1 = self.evaluator.rotate_vector(self.multiply_vect(crypt, [1 if i < (k + self.M() - self.l) * self.N() else 0 for i in range(self.M() * self.N())]), (k - self.l) * self.N(), self.gal_keys)
            wB2 = self.evaluator.rotate_vector(self.multiply_vect(crypt, [1 if k * self.N() <= i < self.l * self.N() else 0 for i in range(self.M() * self.N())]), k * self.N(), self.gal_keys)
            return self.add(wB1, wB2)
        else:
            return crypt
    
    def multiplica(self, a, b):
        c = self.multiply(a, b)
        js = set()
        for i in range(self.t()):
            js.add(i * self.p())
        for k in range(1, self.p()):
            ctA = self.rotar_col(a, k)
            ctB = self.rotar_fil(b, k)
            ctCtemp = self.multiply(ctA, ctB)
            mask = []
            for i in range(self.t()):
                j = (k + i * self.p())% self.l
                if j in js:
                    mask += [0] * self.m * self.n
                else:
                    mask += [1] * self.m * self.n
                    js.add(j)
            if any(map(lambda x: x == 0, mask)):
                ctCtemp = self.multiply_vect(ctCtemp, mask)
            c = self.add(c, ctCtemp)
        return c

##################################################################################################
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
            if crypt.coeff_modulus_size() > 1 and crypt.scale() / (2**self.coeff_modulus_array[crypt.coeff_modulus_size()-1]) >= self.min_scale:
                self.evaluator.rescale_to_next_inplace(crypt)
                crypt.scale(2**int(math.log2(crypt.scale())))
                self.evaluator.mod_switch_to_next_inplace(cod)
                print(f"\t\tcrypt - rescale_to: {math.log2(crypt.scale())} | {crypt.size()}/{crypt.size_capacity()} | {crypt.coeff_modulus_size()} | {crypt.parms_id()}")
                print(f"\t\t\terror:{self.compara(crypt, vectDecrypt)}")
                print(f"\t\tcod - mod_switch_to: {math.log2(cod.scale())} | {cod.coeff_count()} | {cod.parms_id()}")
                print(f"\t\t\terror:{self.comparaDec(cod, vectDecod)}")
            else:
                raise Exception("Se ha descendido de la escala deseada")

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
            if crypt1.scale() - 2**self.coeff_modulus_array[crypt1.coeff_modulus_size()-1] >= self.min_scale:
                self.evaluator.rescale_to_next_inplace(crypt1)
                crypt1.scale(2**int(math.log2(crypt1.scale())))
                print(f"\t\tcrypt1 - rescale_to: {math.log2(crypt1.scale())} | {crypt1.size()}/{crypt1.size_capacity()} | {crypt1.coeff_modulus_size()} | {crypt1.parms_id()}")
            elif crypt1.scale() <= self.max_scales[crypt1.coeff_modulus_size()-1]:
                self.evaluator.mod_switch_to_next_inplace(crypt1)
                print(f"\t\tcrypt1 - mod_switch_to: {math.log2(crypt1.scale())} | {crypt1.size()}/{crypt1.size_capacity()} | {crypt1.coeff_modulus_size()} | {crypt1.parms_id()}")
            else:
                raise Exception("No se puede equiparar el 'coeff_modulus_size' sin descender de la escala deseada ni superar la escala restante")
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
        
        # Comprueba que la multiplicación no exceda la escala maxima y si lo hace realiza un rescalado
        while crypt1.scale()*crypt2.scale() > self.max_scales[crypt1.coeff_modulus_size()]:
            if crypt1.scale() > crypt2.scale():
                if crypt1.coeff_modulus_size() > 1 and crypt1.scale() / (2**self.coeff_modulus_array[crypt1.coeff_modulus_size()-1]) >= self.min_scale:
                    self.evaluator.rescale_to_next_inplace(crypt1)
                    crypt1.scale(2**int(math.log2(crypt1.scale())))
                    print(f"\t\tcrypt1 - rescale_to: {math.log2(crypt1.scale())} | {crypt1.size()}/{crypt1.size_capacity()} | {crypt1.coeff_modulus_size()} | {crypt1.parms_id()}")
                else:
                    raise Exception(f"Precisión mínima necesaria ({math.ceil(math.log2(crypt1.scale() / (2**self.coeff_modulus_array[crypt1.coeff_modulus_size()-1])))})")
                if crypt2.coeff_modulus_size() > 1 and crypt2.scale() / (2**self.coeff_modulus_array[crypt2.coeff_modulus_size()-1]) >= self.min_scale:
                    self.evaluator.rescale_to_next_inplace(crypt2)
                    crypt2.scale(2**int(math.log2(crypt2.scale())))
                    print(f"\t\tcrypt2 - rescale_to: {math.log2(crypt2.scale())} | {crypt2.size()}/{crypt2.size_capacity()} | {crypt2.coeff_modulus_size()} | {crypt2.parms_id()}")
                else:
                    self.evaluator.mod_switch_to_next_inplace(crypt2)
                    print(f"\t\tcrypt2 - mod_switch_to: {math.log2(crypt2.scale())} | {crypt2.size()}/{crypt2.size_capacity()} | {crypt2.coeff_modulus_size()} | {crypt2.parms_id()}")
            elif crypt1.scale() < crypt2.scale():
                if crypt1.coeff_modulus_size() > 1 and crypt1.scale() / (2**self.coeff_modulus_array[crypt1.coeff_modulus_size()-1]) >= self.min_scale:
                    self.evaluator.rescale_to_next_inplace(crypt1)
                    crypt1.scale(2**int(math.log2(crypt1.scale())))
                    print(f"\t\tcrypt1 - rescale_to: {math.log2(crypt1.scale())} | {crypt1.size()}/{crypt1.size_capacity()} | {crypt1.coeff_modulus_size()} | {crypt1.parms_id()}")
                else:
                    self.evaluator.mod_switch_to_next_inplace(crypt1)
                    print(f"\t\tcrypt1 - mod_switch_to: {math.log2(crypt1.scale())} | {crypt1.size()}/{crypt1.size_capacity()} | {crypt1.coeff_modulus_size()} | {crypt1.parms_id()}")
                if crypt2.coeff_modulus_size() > 1 and crypt2.scale() / (2**self.coeff_modulus_array[crypt2.coeff_modulus_size()-1]) >= self.min_scale:
                    self.evaluator.rescale_to_next_inplace(crypt2)
                    crypt2.scale(2**int(math.log2(crypt2.scale())))
                    print(f"\t\tcrypt2 - rescale_to: {math.log2(crypt2.scale())} | {crypt2.size()}/{crypt2.size_capacity()} | {crypt2.coeff_modulus_size()} | {crypt2.parms_id()}")
                else:
                    raise Exception(f"Precisión mínima necesaria ({math.ceil(math.log2(crypt2.scale() / (2**self.coeff_modulus_array[crypt2.coeff_modulus_size()-1])))})")
            else:
                if crypt1.coeff_modulus_size() > 1 and crypt1.scale() / (2**self.coeff_modulus_array[crypt1.coeff_modulus_size()-1]) >= self.min_scale:
                    self.evaluator.rescale_to_next_inplace(crypt1)
                    crypt1.scale(2**int(math.log2(crypt1.scale())))
                    print(f"\t\tcrypt1 - rescale_to: {math.log2(crypt1.scale())} | {crypt1.size()}/{crypt1.size_capacity()} | {crypt1.coeff_modulus_size()} | {crypt1.parms_id()}")
                else:
                    raise Exception(f"Precisión mínima necesaria ({math.ceil(math.log2(crypt1.scale() / (2**self.coeff_modulus_array[crypt1.coeff_modulus_size()-1])))})")
                if crypt2.coeff_modulus_size() > 1 and crypt2.scale() / (2**self.coeff_modulus_array[crypt2.coeff_modulus_size()-1]) >= self.min_scale:
                    self.evaluator.rescale_to_next_inplace(crypt2)
                    crypt2.scale(2**int(math.log2(crypt2.scale())))
                    print(f"\t\tcrypt2 - rescale_to: {math.log2(crypt2.scale())} | {crypt2.size()}/{crypt2.size_capacity()} | {crypt2.coeff_modulus_size()} | {crypt2.parms_id()}")
                else:
                    raise Exception(f"Precisión mínima necesaria ({math.ceil(math.log2(crypt2.scale() / (2**self.coeff_modulus_array[crypt2.coeff_modulus_size()-1])))})")
            
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
            if crypt1.scale() - 2**self.coeff_modulus_array[crypt1.coeff_modulus_size()-1] >= self.min_scale:
                self.evaluator.rescale_to_next_inplace(crypt1)
                crypt1.scale(2**int(math.log2(crypt1.scale())))
                print(f"\t\tcrypt1 - rescale_to: {math.log2(crypt1.scale())} | {crypt1.size()}/{crypt1.size_capacity()} | {crypt1.coeff_modulus_size()} | {crypt1.parms_id()}")
            elif crypt1.scale() <= self.max_scales[crypt1.coeff_modulus_size()-1]:
                self.evaluator.mod_switch_to_next_inplace(crypt1)
                print(f"\t\tcrypt1 - mod_switch_to: {math.log2(crypt1.scale())} | {crypt1.size()}/{crypt1.size_capacity()} | {crypt1.coeff_modulus_size()} | {crypt1.parms_id()}")
            else:
                raise Exception("No se puede equiparar el 'coeff_modulus_size' sin descender de la escala deseada ni superar la escala restante")
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

    def rotar_col_coment(self, crypt, k, vectDecrypt):
        k %= self.N()
        if k != 0:
            print(f"decA1({k}):\n{[1 if i % self.N() < k + self.N() - self.l else 0 for i in range(self.M() * self.N())]}")
            eA1 = self.evaluator.rotate_vector(self.multiply_vect_coment(crypt, [1 if i % self.N() < k + self.N() - self.l else 0 for i in range(self.M() * self.N())], vectDecrypt), k - self.l, self.gal_keys)
            decA1 = rotar_list(mult_list(vectDecrypt, [1 if i % self.N() < k + self.N() - self.l else 0 for i in range(self.M() * self.N())]), k - self.l)
            print(f"{vectDecrypt}\n{decA1}")
            
            eA2 = self.evaluator.rotate_vector(self.multiply_vect_coment(crypt, [1 if k <= i % self.N() < self.l else 0 for i in range(self.M() * self.N())], vectDecrypt), k, self.gal_keys)
            decA2 = rotar_list(mult_list(vectDecrypt, [1 if k <= i % self.N() < self.l else 0 for i in range(self.M() * self.N())]), k)
            print(f"decA2({k}):\n{[1 if k <= i % self.N() < self.l else 0 for i in range(self.M() * self.N())]}\n{vectDecrypt}\n{decA2}")
            return self.add_coment(eA1, eA2, decA1, decA2), sum_list(decA1, decA2)
        else:
            return crypt, vectDecrypt
# ERROR: FALTA ARREGLAR EL CASO (m == n < l) - NO FUNCIONA
    def rotar_fil_coment(self, crypt, k, vectDecrypt):
        k %= self.M()
        if k != 0:
            wB1 = self.evaluator.rotate_vector(self.multiply_vect_coment(crypt, [1 if i < (k + self.M() - self.l) * self.N() else 0 for i in range(self.M() * self.N())], vectDecrypt), (k - self.l) * self.N(), self.gal_keys)
            decB1 = rotar_list(mult_list(vectDecrypt, [1 if i < (k + self.M() - self.l) * self.N() else 0 for i in range(self.M() * self.N())]), (k - self.l) * self.N())
            #print(f"decB1({k}):\n{[1 if i < (k + self.M() - self.l) * self.N() else 0 for i in range(self.M() * self.N())]}\n{vectDecrypt}\n{decB1}")
            
            wB2 = self.evaluator.rotate_vector(self.multiply_vect_coment(crypt, [1 if k * self.N() <= i < self.l * self.N() else 0 for i in range(self.M() * self.N())], vectDecrypt), k * self.N(), self.gal_keys)
            decB2 = rotar_list(mult_list(vectDecrypt, [1 if k * self.N() <= i < self.l * self.N() else 0 for i in range(self.M() * self.N())]), k * self.N())
            #print(f"decB2({k}):\n{[1 if k * self.N() <= i < self.l * self.N() else 0 for i in range(self.M() * self.N())]}\n{vectDecrypt}\n{decB2}")
            return self.add_coment(wB1, wB2, decB1, decB2), sum_list(decB1, decB2)
        else:
            return crypt, vectDecrypt
    
    def multiplica_coment(self, a, b, vA, vB):
        c = self.multiply_coment(a, b, vA, vB)
        vC = mult_list(vA, vB)
        #print(f"vC[0]: [{vC}]")
        js = set()
        for i in range(self.t()):
            js.add(i * self.p())
            #print(f"j: {i * self.p()}")
        #print(f"js: {js}")
        for k in range(1, self.p()):
            ctA, vAk = self.rotar_col_coment(a, k, vA)
            #print(f"vAk[{k}]: {vAk}")
            ctB, vBk  = self.rotar_fil_coment(b, k, vB)
            #print(f"vBk[{k}]: {vBk}")
            ctCtemp = self.multiply_coment(ctA, ctB, vAk, vBk)
            vCtemp = mult_list(vAk, vBk)
            #print(f"vCtemp[{k}]: {vCtemp}")
            mask = []
            for i in range(self.t()):
                j = (k + i * self.p())% self.l
                #print(f"j: {j}")
                if j in js:
                    mask += [0] * self.m * self.n
                else:
                    mask += [1] * self.m * self.n
                    js.add(j)
                #print(mask)
            if any(map(lambda x: x == 0, mask)):
                ctCtemp = self.multiply_vect_coment(ctCtemp, mask, vCtemp)
                vCtemp = mult_list(vCtemp, mask)
            #print(f"js: {js}")
            c = self.add_coment(c, ctCtemp, vC, vCtemp)
            vC = sum_list(vC, vCtemp)
            #print(f"vC[{k}]: {vC}")
        #print(self.compara(c, vC))
        return c
    
    def compara(self, ct, v):
        return self.comparaDec(self.decryptor.decrypt(ct), v)
    
    def comparaDec(self, cod, v):
        cod = self.encoder.decode(cod)[:len(v)]
        t = abs(cod-v)
        return f"\t±{numpy.mean(t)}"
    
