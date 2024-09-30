import os
import psutil
import time
import math
import numpy as np

from .server import Server
from .client import Client


def eva_test(matriz1, matriz2, numeros, params={'enc':128, 'pre':40}):
    if not params.get("silence", False):
        print("\nSTARTING EVA...")
    
    ret={}
    ret["ini"] = time.time()

    ret["lib"]="EVA"
    ret["tam_matr1"] = f"{len(matriz1)}x{len(matriz1[0])}"
    ret["tam_matr2"] = f"{len(matriz2)}x{len(matriz2[0])}"
    #ret["poly_modulus_degree"]=params["pmd"]
    ret["encriptacion"] = params["enc"]
    ret["precision"] = params["pre"]
    ret["scheme"] = 'CKKS'
    ret["encript_ctx"] = params["enc"]
    try:
        lado_mat_cuad=2**math.ceil(math.log2(max([len(matriz1),len(matriz1[0]),len(matriz2),len(matriz2[0])])))
        process = psutil.Process(os.getpid())
        m=process.memory_info().rss
        t=time.time()
        ser = Server()
        parametros, signature = ser.compileProgram(lado_mat_cuad, params["pre"], params["pre"], params["enc"], load_if_can=False)
        cli = Client()
        public_ctx = cli.keyGen(parametros)
        encInputs = cli.encrypt(matriz1, matriz2, signature)
        encOutputs = ser.execute(public_ctx, encInputs)
        outputs = cli.decrypt(encOutputs, tamOut=(len(matriz1),len(matriz2[0])))#, signature)
        ret["t_total"]=time.time()-t
        ret["m_total"]=process.memory_info().rss-m
        matMultEsp = np.matmul(matriz1, matriz2)
        ret["e_program"] = np.mean(abs((outputs-matMultEsp))/matMultEsp)

        ################################################################################################
        
        ser = Server(test=True)

        parametros, signature, tests = ser.compileProgram(lado_mat_cuad, params["pre"], params["pre"], params["enc"], load_if_can=False)
        ret.update(tests)
        ret["poly_modulus_degree"] = parametros.poly_modulus_degree
        ret["poly_modulus_degree_ctx"] = parametros.poly_modulus_degree
        ret["rotations"] = parametros.rotations
        ret["coeff_modulus"] = sum(parametros.prime_bits)
        ret["coeff_modulus_arr"] = parametros.prime_bits
        
        cli = Client(test=True)
        public_ctx, tests = cli.keyGen(parametros)
        ret.update(tests)
        encInputs, tests = cli.encrypt(matriz1, matriz2, signature)
        ret.update(tests)
        encOutputs, tests = ser.execute(public_ctx, encInputs)
        ret.update(tests)
        outputs, tests = cli.decrypt(encOutputs, tamOut=(len(matriz1),len(matriz2[0])))#, signature)
        ret.update(tests)
        ret["e_program"] = (ret["e_program"]+np.mean(abs((outputs-matMultEsp))/matMultEsp))/2
    except Exception as error:
        ret["errors"] = ret.get("errors", []) + [f"{time.time()}: EVA–{type(error).__name__}–{error}"]
        print(f"Error en EVA: {type(error).__name__}–{error}")

    ret["fin"] = time.time()
    
    if not params.get("silence", False):
        for k in ret.keys():
            print(f"{k}: {ret[k]}")

    return ret
