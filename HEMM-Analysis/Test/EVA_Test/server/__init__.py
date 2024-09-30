import os
import psutil
import time
import math
import numpy

from eva import EvaProgram, Input, Output, save, load#, evaluate
from eva.ckks import CKKSCompiler
#from eva.seal import generate_keys
#from eva.metric import valuation_mse

process = psutil.Process(os.getpid())

MULTI_CORE = False
################################################
# https://link.springer.com/article/10.1007/s11227-022-04850-4

def replica_diag(ct, i, a_tam, hc_tam, v_size):
    for j in range(a_tam[0]):
        if j == 0:
            ct_m = ct * ([0]*hc_tam[1]*j+[1]*hc_tam[1]+[0]*(v_size-hc_tam[1]*(j+1))) << (j+i)%a_tam[1]
        else: 
            ct_m = ct_m + ct * ([0]*hc_tam[1]*j+[1]*hc_tam[1]+[0]*(v_size-hc_tam[1]*(j+1))) << (j+i)%a_tam[1]
    e = 1
    ret = ct_m
    for j in reversed(range(len(bin(hc_tam[1]))-3)):
        ret = ret + (ret << -e)
        e=e*2
        if bin(hc_tam[1])[-(j+1)]== "1":
            ret = ct_m + (ret << -1)
            e=e+1
    return ret

def muti_matriz(a, m, l, b, l2, n, v_size):
    res = -1
    if(l==l2):
        if m <= l:
            #Caso A
            t = max(l, n)
            for i in range(l):
                if i==0:
                    res = b * replica_diag(a*[1 if j<l*t and (j%t) == ((i+math.floor(j/t))%l) else 0 for j in range(v_size)], i, (m, l), (l, t), v_size)
                else:
                    if v_size == l * n:
                        res = res + (b << i*t) * replica_diag(a*[1 if j<l*t and (j%t) == ((i+math.floor(j/t))%l) else 0 for j in range(v_size)], i, (m, l), (l, t), v_size)
                    else:
                        res = res + ((b*([0]*t*i+[1]*(v_size-t*i)) << t*i) + (b*([1]*t*i+[0]*(v_size-t*i)) << -t*(l-i))) * replica_diag(a*[1 if j<l*t and (j%t) == ((i+math.floor(j/t))%l) else 0 for j in range(v_size)], i, (m, l), (l, t), v_size)
        else:
            #Caso B
            k = math.ceil(m/l)*l 
            t = max(l, n)

            
            res = k
    return res

#################################################
# https://eprint.iacr.org/2018/1041.pdf
# https://github.com/mark-schultz/EVA-matrix-multiplication

def vec_from_pred(n, pred):
    # Returns a vector v in {0,1}^n s.t.
    # v[i] = pred(i)
    return [1 if pred(ell) else 0 for ell in range(n)]

def heMatMult(a, b, d):
    n = d**2

    if MULTI_CORE:
        # Step 1-1
        ctA = [0]*d
        ctA[0] = [0]*n
        for k in range(-d-1, d):
            if k >= 0:
                uk = vec_from_pred(n, lambda ell: 0 <= ell - d * k < (d-k))
            else:
                uk = vec_from_pred(n, lambda ell: -k <= ell - (d+k) * d < d)
            if any(map(lambda e: e != 0, uk)):
                ctA[0] += (a << k) * uk

        # Step 1-2
        ctB = [0]*d
        ctB[0] = [0]*n
        for k in range(d):
            uk = vec_from_pred(n, lambda ell: ell % d == k)
            if any(map(lambda e: e != 0, uk)):
                ctB[0] += (b << (d * k)) * uk
        
        # Step 2
        for k in range(1, d):
            vk = vec_from_pred(n, lambda ell: 0 <= ell % d < d - k)
            vk_minus_d = vec_from_pred(n, lambda ell: (d-k) <= ell % d < d)
            vk_not_empty = any(map(lambda e: e != 0, vk))
            vk_minus_d_not_empty = any(map(lambda e: e != 0, vk_minus_d))
            if vk_not_empty and vk_minus_d_not_empty:
                ctA[k] = (ctA[0] << k) * vk + (ctA[0] << (k-d)) * vk_minus_d
            elif vk_not_empty:
                ctA[k] = (ctA[0] << k) * vk
            elif vk_minus_d_not_empty:
                ctA[k] = (ctA[0] << (k-d)) * vk_minus_d
            else:
                ctA[k] = [0]*n
            ctB[k] = (ctB[0] << (d * k))

        # Step 3
        ctAB = ctA[0] * ctB[0]
        for k in range(1, d):
            ctAB += ctA[k] * ctB[k]
    else:
        # Step 1-1
        ctA0 = [0]*n
        for k in range(-d-1, d):
            if k >= 0:
                uk = vec_from_pred(n, lambda ell: 0 <= ell - d * k < (d-k))
            else:
                uk = vec_from_pred(n, lambda ell: -k <= ell - (d+k) * d < d)
            if any(map(lambda e: e != 0, uk)):
                ctA0 += (a << k) * uk

        # Step 1-2
        ctB0 = [0]*n
        for k in range(d):
            uk = vec_from_pred(n, lambda ell: ell % d == k)
            if any(map(lambda e: e != 0, uk)):
                ctB0 += (b << (d * k)) * uk
        
        # Step 2 & 3
        ctAB = ctA0 * ctB0
        for k in range(1, d):
            vk = vec_from_pred(n, lambda ell: 0 <= ell % d < d - k)
            vk_minus_d = vec_from_pred(n, lambda ell: (d-k) <= ell % d < d)
            vk_not_empty = any(map(lambda e: e != 0, vk))
            vk_minus_d_not_empty = any(map(lambda e: e != 0, vk_minus_d))
            if vk_not_empty and vk_minus_d_not_empty:
                ctAB += ((ctA0 << k) * vk + (ctA0 << (k-d)) * vk_minus_d) * (ctB0 << (d * k))
            elif vk_not_empty:
                ctAB += ((ctA0 << k) * vk) * (ctB0 << (d * k))
            elif vk_minus_d_not_empty:
                ctAB += ((ctA0 << (k-d)) * vk_minus_d) * (ctB0 << (d * k))
    
    return ctAB

#################################################
class Server:
    def __init__(self, test = False):
        self.program = None
        self.test = test
    
    #def __init__(self, params, signature):
    #    self.public_ctx, self.secret_ctx = generate_keys(params)
    #    self.signature = signature

    def compileProgram(self, d, ranges, scales, security_level=128, d_max=8, load_if_can=True, save_prog=False):
        if self.test:
            if load_if_can and os.path.isfile(f'program.{d}x{d}.{d_max}.eva') and os.path.isfile(f'program.{d}x{d}.{d_max}.evaparams') and os.path.isfile(f'program.{d}x{d}.{d_max}.evasignature'):
                self.program = load(f'program.{d}x{d}.{d_max}.eva')
                params = load(f'program.{d}x{d}.{d_max}.evaparams')
                signature = load(f'program.{d}x{d}.{d_max}.evasignature')
            else:
                config = {
                    #balance_reductions - Balance trees of mul, add or sub operations. bool (default=true)
                    #rescaler           - Rescaling policy. One of: lazy_waterline (default), eager_waterline, always, minimum
                    #lazy_relinearize   - Relinearize as late as possible. bool (default=true)
                    #security_level     - How many bits of security parameters should be selected for. int (default=128)
                    #quantum_safe       - Select quantum safe parameters. bool (default=false)
                    #warn_vec_size      - Warn about possibly inefficient vector size selection. bool (default=true)
                    'warn_vec_size':'false',
                    'security_level': str(security_level)
                }
                
                tests = {}
                tests["config"] = str(config)
                tests["program"] = f'Matrix {d}x{d} Multiplication'
                m=process.memory_info().rss
                t=time.time()
                program = EvaProgram(f'Matrix {d}x{d} Multiplication', vec_size=d**2)
                with program:
                    a = Input('a')
                    b = Input('b')
                    #c = Input('c')
                    #c = c + heMatMult(a, b, d)
                    #Output('c', c)
                    Output('c', heMatMult(a, b, d))
                    #Output('c', muti_matriz(a, d, d, b, d, d, program.vec_size))

                program.set_output_ranges(ranges)
                program.set_input_scales(scales)
                tests["t_program"]=time.time()-t
                tests["m_program"]=process.memory_info().rss-m

                m=process.memory_info().rss
                t=time.time()
                compiler = CKKSCompiler(config=config)
                tests["t_compiler"]=time.time()-t
                tests["m_compiler"]=process.memory_info().rss-m
                m=process.memory_info().rss
                t=time.time()
                self.program, params, signature = compiler.compile(program)
                tests["t_compile"]=time.time()-t
                tests["m_compile"]=process.memory_info().rss-m

                #print(program.to_DOT())
                if save_prog:
                    save(program, f'program.{d}x{d}.{d_max}.eva')
                    save(params, f'program.{d}x{d}.{d_max}.evaparams')
                    save(signature, f'program.{d}x{d}.{d_max}.evasignature')

            return params, signature, tests
        else:
            if load_if_can and os.path.isfile(f'program.{d}x{d}.{d_max}.eva') and os.path.isfile(f'program.{d}x{d}.{d_max}.evaparams') and os.path.isfile(f'program.{d}x{d}.{d_max}.evasignature'):
                self.program = load(f'program.{d}x{d}.{d_max}.eva')
                params = load(f'program.{d}x{d}.{d_max}.evaparams')
                signature = load(f'program.{d}x{d}.{d_max}.evasignature')
            else:
                config = {
                    #balance_reductions - Balance trees of mul, add or sub operations. bool (default=true)
                    #rescaler           - Rescaling policy. One of: lazy_waterline (default), eager_waterline, always, minimum
                    #lazy_relinearize   - Relinearize as late as possible. bool (default=true)
                    #security_level     - How many bits of security parameters should be selected for. int (default=128)
                    #quantum_safe       - Select quantum safe parameters. bool (default=false)
                    #warn_vec_size      - Warn about possibly inefficient vector size selection. bool (default=true)
                    'warn_vec_size':'false',
                    'security_level': str(security_level)
                }
                
                program = EvaProgram(f'Matrix {d}x{d} Multiplication', vec_size=d**2)
                with program:
                    a = Input('a')
                    b = Input('b')
                    #c = Input('c')
                    #c = c + heMatMult(a, b, d)
                    #Output('c', c)
                    Output('c', heMatMult(a, b, d))
                    #Output('c', muti_matriz(a, d, d, b, d, d, program.vec_size))

                program.set_output_ranges(ranges)
                program.set_input_scales(scales)
                compiler = CKKSCompiler(config=config)
                self.program, params, signature = compiler.compile(program)

                #print(program.to_DOT())
                if save_prog:
                    save(program, f'program.{d}x{d}.{d_max}.eva')
                    save(params, f'program.{d}x{d}.{d_max}.evaparams')
                    save(signature, f'program.{d}x{d}.{d_max}.evasignature')

            return params, signature

    def execute(self, public_ctx, encInputs):
        if self.test:
            tests = {}

            m=process.memory_info().rss
            t=time.time()
            encOutputs = public_ctx.execute(self.program, encInputs)
            tests["t_execute"]=time.time()-t
            tests["m_execute"]=process.memory_info().rss-m

            return encOutputs, tests
        else:
            return public_ctx.execute(self.program, encInputs)
