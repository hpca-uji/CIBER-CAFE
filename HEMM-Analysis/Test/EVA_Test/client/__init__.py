import os
import psutil
import time
import math
import numpy

#from eva import Input, Output, EvaProgram, save, load, evaluate
#from eva.ckks import CKKSCompiler
from eva.seal import generate_keys
#from eva.metric import valuation_mse

process = psutil.Process(os.getpid())

class Client:
    def __init__(self, test = False):
        self.secret_ctx = None
        self.public_ctx = None
        self.signature = None
        self.test = test
    
    #def __init__(self, params, signature):
    #    self.public_ctx, self.secret_ctx = generate_keys(params)
    #    self.signature = signature

    def keyGen(self, params):
        if self.test:
            tests = {}

            m=process.memory_info().rss
            t=time.time()
            self.public_ctx, self.secret_ctx = generate_keys(params)
            tests["t_generate_keys"]=time.time()-t
            tests["m_generate_keys"]=process.memory_info().rss-m

            return self.public_ctx, tests
        else:
            self.public_ctx, self.secret_ctx = generate_keys(params)
            return self.public_ctx

    def encrypt(self, matriz1, matriz2, signature):
        if self.test:
            tests = {}
            self.signature = signature
            d = int(math.sqrt(signature.vec_size))

            m=process.memory_info().rss
            t=time.time()
            matriz1_in = [0]*signature.vec_size
            m1_col=len(matriz1[0])
            for i in range(len(matriz1)):
                for j in range(m1_col):
                    matriz1_in[i*d+j] = matriz1[i][j]
            
            matriz2_in = [0]*signature.vec_size
            m2_col=len(matriz2[0])
            for i in range(len(matriz2)):
                for j in range(m2_col):
                    matriz2_in[i*d+j] = matriz2[i][j]

            inputs = {
                'a': matriz1_in,
                'b': matriz2_in
            }
            tests["t_prep_inputs"]=time.time()-t
            tests["m_prep_inputs"]=process.memory_info().rss-m

            m=process.memory_info().rss
            t=time.time()
            encInputs = self.public_ctx.encrypt(inputs, signature)
            tests["t_encrypt"]=time.time()-t
            tests["m_encrypt"]=process.memory_info().rss-m

            return encInputs, tests
        else:
            self.signature = signature
            d = int(math.sqrt(signature.vec_size))
            matriz1_in = [0]*signature.vec_size
            m1_col=len(matriz1[0])
            for i in range(len(matriz1)):
                for j in range(m1_col):
                    matriz1_in[i*d+j] = matriz1[i][j]
            
            matriz2_in = [0]*signature.vec_size
            m2_col=len(matriz2[0])
            for i in range(len(matriz2)):
                for j in range(m2_col):
                    matriz2_in[i*d+j] = matriz2[i][j]

            inputs = {
                'a': matriz1_in,
                'b': matriz2_in
            }

            encInputs = self.public_ctx.encrypt(inputs, signature)

            return encInputs

    def decrypt(self, encOutputs, signature=None, tamOut=None):
        if self.test:
            tests = {}
            if signature:
                self.signature = signature
            m=process.memory_info().rss
            t=time.time()
            outputsDecode = self.secret_ctx.decrypt(encOutputs, self.signature)
            tests["t_decrypt"]=time.time()-t
            tests["m_decrypt"]=process.memory_info().rss-m

            m=process.memory_info().rss
            t=time.time()
            d = int(math.sqrt(self.signature.vec_size))
            if (tamOut == None):
                tamOut = (d, d)
            output = numpy.array([outputsDecode['c'][i*d:i*d+tamOut[1]] for i in range(tamOut[0])])
            tests["t_prep_outputs"]=time.time()-t
            tests["m_prep_outputs"]=process.memory_info().rss-m

            return output, tests
        else:
            if signature:
                self.signature = signature
            outputsDecode = self.secret_ctx.decrypt(encOutputs, self.signature)

            d = int(math.sqrt(self.signature.vec_size))
            if (tamOut == None):
                tamOut = (d, d)
            output = numpy.array([outputsDecode['c'][i*d:i*d+tamOut[1]] for i in range(tamOut[0])])

            return output
