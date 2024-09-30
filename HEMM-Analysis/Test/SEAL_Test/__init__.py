import os
import psutil
import time
import math
import numpy

from .seal import *
from .heMatMult import HEMatMult
from .HEGMM import HEGMM

MAX_BITS_COEFF_MODULUS = 60

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

def compara(ct, v, encoder, decryptor):
    ct = encoder.decode(decryptor.decrypt(ct))
    t = abs(ct-v)
    d = math.floor(math.sqrt(len(v)))
    return f"\n\t±{numpy.mean(t)}"#\n\t>{max(t)}\n\t<{min(t)}"#\nct:\n{numpy.array(ct).reshape(d, d)}\ntestigo:\n{numpy.array(v).reshape(d, d)}"



#################################################
#   A larger coeff_modulus implies a larger noise budget, hence more encrypted
#   computation capabilities. However, an upper bound for the total bit-length
#   of the coeff_modulus is determined by the poly_modulus_degree, as follows:
#       +----------------------------------------------------+
#       | poly_modulus_degree | max coeff_modulus bit-length |
#       +---------------------+------------------------------+
#       | 1024                | 27                           |
#       | 2048                | 54                           |
#       | 4096                | 109                          |
#       | 8192                | 218                          |
#       | 16384               | 438                          |
#       | 32768               | 881                          |
#       +---------------------+------------------------------+
encriptacion_types = {
    128: sec_level_type.tc128,
    192: sec_level_type.tc192,
    256: sec_level_type.tc256,
    'None': sec_level_type.none
}

def seal_test(matriz1, matriz2, numeros, params={'pmd':8192, 'enc':128, 'pre':40}):
    if not params.get("silence", False):
        print("\nSTARTING SEAL...")
    ret={}
    ret["ini"] = time.time()

    ret["lib"]="SEAL"
    ret["tam_matr1"] = f"{len(matriz1)}x{len(matriz1[0])}"
    ret["tam_matr2"] = f"{len(matriz2)}x{len(matriz2[0])}"
    ret["poly_modulus_degree"]=params["pmd"]
    ret["encriptacion"] = params["enc"]
    ret["precision"] = params["pre"]
    
    try:
        filas_vect1=len(matriz1)
        columnas_vect1=len(matriz1[0])
        vect1 = list(matriz1.reshape(-1))
        vect1 *= math.floor(params["pmd"]/2/(filas_vect1*columnas_vect1))
        
        filas_vect2=len(matriz2)
        columnas_vect2=len(matriz2[0])
        vect2 = list(matriz2.reshape(-1))
        vect2 *= math.floor(params["pmd"]/2/(filas_vect2*columnas_vect2))
        
        scale=pow(2.0, params["pre"])
        process = psutil.Process(os.getpid())

        # Coeff_Modulus
        m=process.memory_info().rss
        t=time.time()
        if "cmp" in params:
            coeff_modulus = list(map(lambda x: Modulus(x), params["cmp"]))
        elif "cma" in params:
            coeff_modulus = CoeffModulus.Create(params["pmd"], params["cma"])
        elif "cm" in params:
            n_coeff_modulus=math.ceil(params["cm"]/MAX_BITS_COEFF_MODULUS)
            v_coeff_modulus=[math.floor(params["cm"]/n_coeff_modulus)] * n_coeff_modulus
            for i in range(n_coeff_modulus-math.floor(params["cm"]%n_coeff_modulus), len(v_coeff_modulus)):
                v_coeff_modulus[i] += 1
            coeff_modulus = CoeffModulus.Create(params["pmd"], v_coeff_modulus)
        else:
            if 1024 <= params["pmd"] <= 32768:
                coeff_modulus = CoeffModulus.BFVDefault(params["pmd"], encriptacion_types.get(params["enc"], sec_level_type.tc128))
                #coeff_modulus = CoeffModulus.Create(params["pmd"], [60, 60, 60, 60, 60, 60, 60, 60, 60])
            elif params["pmd"] > 32768:
                b_coeff_modulus=params["pmd"]/32768*sum(map(lambda x: x.bit_count(), CoeffModulus.BFVDefault(32768, encriptacion_types.get(params["enc"], sec_level_type.tc128))))
                n_coeff_modulus=math.ceil(b_coeff_modulus/MAX_BITS_COEFF_MODULUS)
                v_coeff_modulus=[math.floor(b_coeff_modulus/n_coeff_modulus)] * n_coeff_modulus
                for i in range(n_coeff_modulus-math.floor(b_coeff_modulus%n_coeff_modulus), len(v_coeff_modulus)):
                    v_coeff_modulus[i] += 1
                coeff_modulus = CoeffModulus.Create(params["pmd"], v_coeff_modulus)
            elif params["pmd"] < 1024:
                b_coeff_modulus=params["pmd"]/1024*sum(map(lambda x: x.bit_count(), CoeffModulus.BFVDefault(1024, encriptacion_types.get(params["enc"], sec_level_type.tc128))))
                n_coeff_modulus=math.ceil(b_coeff_modulus/MAX_BITS_COEFF_MODULUS)
                v_coeff_modulus=[math.floor(b_coeff_modulus/n_coeff_modulus)] * n_coeff_modulus
                for i in range(n_coeff_modulus-math.floor(b_coeff_modulus%n_coeff_modulus), len(v_coeff_modulus)):
                    v_coeff_modulus[i] += 1
                coeff_modulus = CoeffModulus.Create(params["pmd"], v_coeff_modulus)

        ret["t_prep_coeff_modulus"]=time.time()-t
        ret["m_prep_coeff_modulus"]=process.memory_info().rss-m

        # Parametros
        m=process.memory_info().rss
        t=time.time()
        parms = EncryptionParameters(scheme_type.ckks)
        parms.set_poly_modulus_degree(params["pmd"])
        parms.set_coeff_modulus(coeff_modulus)
        ret["t_prep_params"]=time.time()-t
        ret["m_prep_params"]=process.memory_info().rss-m
        
        # Preparativos
        m=process.memory_info().rss
        t=time.time()
        context = SEALContext(parms, True, encriptacion_types.get(params["enc"], sec_level_type.tc128) if params["pmd"] <= 32768 else sec_level_type.none)
        ret["t_prep_context"]=time.time()-t
        ret["m_prep_context"]=process.memory_info().rss-m
        m=process.memory_info().rss
        t=time.time()
        keygen = KeyGenerator(context)
        ret["t_prep_keygen"]=time.time()-t
        ret["m_prep_keygen"]=process.memory_info().rss-m
        m=process.memory_info().rss
        t=time.time()
        secret_key = keygen.secret_key()
        ret["t_prep_secret_key"]=time.time()-t
        ret["m_prep_secret_key"]=process.memory_info().rss-m
        m=process.memory_info().rss
        t=time.time()
        public_key = keygen.create_public_key()
        ret["t_prep_public_key"]=time.time()-t
        ret["m_prep_public_key"]=process.memory_info().rss-m
        try:
            m=process.memory_info().rss
            t=time.time()
            relin_keys = keygen.create_relin_keys()
            ret["t_prep_relin_keys"]=time.time()-t
            ret["m_prep_relin_keys"]=process.memory_info().rss-m
        except Exception as error:
            relin_keys = None
            ret["errors"] = ret.get("errors", []) + [f"{time.time()}: Realin Key–{type(error).__name__}–{error}"]
            print(f"SEAL-No se pudo crear la clave de Realineamiento: {type(error).__name__}–{error}")
        try:
            m=process.memory_info().rss
            t=time.time()
            gal_keys = keygen.create_galois_keys()
            ret["t_prep_gal_keys"]=time.time()-t
            ret["m_prep_gal_keys"]=process.memory_info().rss-m
        except Exception as error:
            gal_keys = None
            ret["errors"] = ret.get("errors", []) + [f"{time.time()}: Galois Key–{type(error).__name__}–{error}"]
            print(f"SEAL-No se pudo crear la clave de Galois: {type(error).__name__}–{error}")
        m=process.memory_info().rss
        t=time.time()
        encryptor = Encryptor(context, public_key)
        ret["t_prep_encryptor"]=time.time()-t
        ret["m_prep_encryptor"]=process.memory_info().rss-m
        m=process.memory_info().rss
        t=time.time()
        evaluator = Evaluator(context)
        ret["t_prep_evaluator"]=time.time()-t
        ret["m_prep_evaluator"]=process.memory_info().rss-m
        m=process.memory_info().rss
        t=time.time()
        decryptor = Decryptor(context, secret_key)#
        ret["t_prep_decryptor"]=time.time()-t
        ret["m_prep_decryptor"]=process.memory_info().rss-m
        m=process.memory_info().rss
        t=time.time()
        encoder = CKKSEncoder(context)
        ret["t_prep_encoder"]=time.time()-t
        ret["m_prep_encoder"]=process.memory_info().rss-m
        ret["t_prep_total"]=sum(ret[x] for x in ret.keys() if x.startswith("t_prep_"))
        ret["m_prep_total"]=sum(ret[x] for x in ret.keys() if x.startswith("m_prep_"))
        
        context_data = context.key_context_data()
        ret["scheme"] = 'CKKS' if context_data.parms().scheme() == scheme_type.ckks else 'BFV'
        ret["poly_modulus_degree_ctx"] = context_data.parms().poly_modulus_degree()
        ret["encript_ctx"] = [k for k in encriptacion_types if encriptacion_types[k]==context_data.qualifiers().sec_level][0]
        ret["coeff_modulus"] = sum(map(lambda x: x.bit_count(), context_data.parms().coeff_modulus()))
        ret["coeff_modulus_arr"] = list(map(lambda x: x.bit_count(), context_data.parms().coeff_modulus()))
        ret["coeff_modulus_primes"] = list(map(lambda x: x.value(), context_data.parms().coeff_modulus()))
        
        #Matriz 1
        m=process.memory_info().rss
        t=time.time()
        matr1_encode = encoder.encode(vect1, scale)
        ret["t_matr1_encode"]=time.time()-t
        ret["m_matr1_encode"]=process.memory_info().rss-m
        m=process.memory_info().rss
        t=time.time()
        matr1_encrypt=encryptor.encrypt(matr1_encode)
        ret["t_matr1_encrypt"]=time.time()-t
        ret["m_matr1_encrypt"]=process.memory_info().rss-m
        ret["t_matr1_total"]=ret["t_matr1_encode"]+ret["t_matr1_encrypt"]
        ret["m_matr1_total"]=ret["m_matr1_encode"]+ret["m_matr1_encrypt"]
        ret["b_matr1_encrypt"]=math.ceil(math.log2(matr1_encrypt.scale()))

        #Matriz 2
        m=process.memory_info().rss
        t=time.time()
        matr2_encode = encoder.encode(vect2, scale)
        ret["t_matr2_encode"]=time.time()-t
        ret["m_matr2_encode"]=process.memory_info().rss-m
        m=process.memory_info().rss
        t=time.time()
        matr2_encrypt=encryptor.encrypt(matr2_encode)
        ret["t_matr2_encrypt"]=time.time()-t
        ret["m_matr2_encrypt"]=process.memory_info().rss-m
        ret["t_matr2_total"]=ret["t_matr2_encode"]+ret["t_matr2_encrypt"]
        ret["m_matr2_total"]=ret["m_matr2_encode"]+ret["m_matr2_encrypt"]
        ret["b_matr2_encrypt"]=math.ceil(math.log2(matr2_encrypt.scale()))

#        print(f"matr1_encrypt: {math.ceil(math.log2(matr1_encrypt.scale()))}b - {matr1_encrypt.parms_id()}")
#        matr1_t = matr1_encrypt
#        for i in range(1, len(ret['coeff_modulus_arr'])+1):
#            evaluator.relinearize_inplace(matr1_t, relin_keys)
#            print(f"matr1_encrypt-{i}ºrelin:\t{math.ceil(math.log2(matr1_t.scale()))}b - {matr1_t.coeff_modulus_size()} - {matr1_t.size()} - {matr1_t.size_capacity()} - {matr1_t.parms_id()}")
#            evaluator.rescale_to_next_inplace(matr1_t)
#            print(f"matr1_encrypt-{i}ºresca:\t{math.ceil(math.log2(matr1_t.scale()))}b - {matr1_t.coeff_modulus_size()} - {matr1_t.size()} - {matr1_t.size_capacity()} - {matr1_t.parms_id()}")
#            evaluator.multiply_inplace(matr1_t, matr1_t)
#            print(f"matr1_encrypt-{i+1}ºsquare:\t{math.ceil(math.log2(matr1_t.scale()))}b - {matr1_t.coeff_modulus_size()} - {matr1_t.size()} - {matr1_t.size_capacity()} - {matr1_t.parms_id()}")

#        # c(Matriz 1) Square x²
#        try:
#            m=process.memory_info().rss
#            t=time.time()
#            matr1_square=evaluator.square(matr1_encrypt)
#            ret["t_matr1_square"]=time.time()-t
#            ret["m_matr1_square"]=process.memory_info().rss-m
#            ret["b_matr1_square"]=math.ceil(math.log2(matr1_square.scale()))
#            if relin_keys!=None:
#                m=process.memory_info().rss
#                t=time.time()
#                evaluator.relinearize_inplace(matr1_square, relin_keys)
#                ret["t_matr1_square_realin"]=time.time()-t
#                ret["m_matr1_square_realin"]=process.memory_info().rss-m
#                ret["b_matr1_square_realin"]=math.ceil(math.log2(matr1_square.scale()))
#            try:
#                m=process.memory_info().rss
#                t=time.time()
#                evaluator.rescale_to_next_inplace(matr1_square)
#                ret["t_matr1_square_rescale"]=time.time()-t
#                ret["m_matr1_square_rescale"]=process.memory_info().rss-m
#                ret["b_matr1_square_rescale"]=math.ceil(math.log2(matr1_square.scale()))
#            except Exception as error:
#                ret["errors"] = ret.get("errors", []) + [f"{time.time()}: c(Matriz 1)^2–Reescalado–{type(error).__name__}–{error}"]
#                print(f"SEAL-Error Reescalando: {type(error).__name__}–{error}")
#            m=process.memory_info().rss
#            t=time.time()
#            matr1_square_dencript=decryptor.decrypt(matr1_square)
#            ret["t_matr1_square_dencript"]=time.time()-t
#            ret["m_matr1_square_dencript"]=process.memory_info().rss-m
#            m=process.memory_info().rss
#            t=time.time()
#            matr1_square_decode=encoder.decode(matr1_square_dencript)
#            ret["t_matr1_square_decode"]=time.time()-t
#            ret["m_matr1_square_decode"]=process.memory_info().rss-m
#            ret["t_matr1_square_total"]=sum(ret[x] for x in ret.keys() if x.startswith("t_matr1_square"))
#            ret["m_matr1_square_total"]=sum(ret[x] for x in ret.keys() if x.startswith("m_matr1_square"))
#            matr1_square_decode=matr1_square_decode[:filas_vect1*columnas_vect1].reshape(filas_vect1,columnas_vect1)
#            ret["e_matr1_square"]=numpy.mean(abs(matriz1*matriz1-matr1_square_decode))
#        except Exception as error:
#            ret["errors"] = ret.get("errors", []) + [f"{time.time()}: c(Matriz 1)^2–{type(error).__name__}–{error}"]
#            print(f"SEAL-c(Matriz 1) Square x²: {type(error).__name__}–{error}")
#        
#        num_encode={}
#        num_encrypt={}
#        
#        matr1_n_plain_multiply={}
#        matr1_n_plain_multiply_dencript={}
#        matr1_n_plain_multiply_decode={}
#
#        matr1_n_multiply={}
#        matr1_n_multiply_dencript={}
#        matr1_n_multiply_decode={}
#
#        try:
#            for k in numeros.keys():
#                # Numeros
#                ret[f"num_{k}"] = numeros[k]
#                m=process.memory_info().rss
#                t=time.time()
#                num_encode[k] = encoder.encode(numeros[k], scale)
#                ret[f"t_num_{k}_encode"]=time.time()-t
#                ret[f"m_num_{k}_encode"]=process.memory_info().rss-m
#                m=process.memory_info().rss
#                t=time.time()
#                num_encrypt[k] = encryptor.encrypt(num_encode[k])
#                ret[f"t_num_{k}_encrypt"]=time.time()-t
#                ret[f"m_num_{k}_encrypt"]=process.memory_info().rss-m
#                ret[f"t_num_{k}_total"]=ret[f"t_num_{k}_encode"]+ret[f"t_num_{k}_encrypt"]
#                ret[f"m_num_{k}_total"]=ret[f"m_num_{k}_encode"]+ret[f"m_num_{k}_encrypt"]
#                ret[f"b_num_{k}_encrypt"]=math.ceil(math.log2(num_encrypt[k].scale()))
#
#                # c(Matriz 1) * Numeros
#                try:
#                    m=process.memory_info().rss
#                    t=time.time()
#                    matr1_n_plain_multiply[k]=evaluator.multiply_plain(matr1_encrypt,num_encode[k])
#                    ret[f"t_matr1_{k}_plain_multiply"]=time.time()-t
#                    ret[f"m_matr1_{k}_plain_multiply"]=process.memory_info().rss-m
#                    ret[f"b_matr1_{k}_plain_multiply"]=math.ceil(math.log2(matr1_n_plain_multiply[k].scale()))
#                    try:
#                        m=process.memory_info().rss
#                        t=time.time()
#                        evaluator.rescale_to_next_inplace(matr1_n_plain_multiply[k])
#                        ret[f"t_matr1_{k}_plain_multiply_rescale"]=time.time()-t
#                        ret[f"m_matr1_{k}_plain_multiply_rescale"]=process.memory_info().rss-m
#                        ret[f"b_matr1_{k}_plain_multiply_rescale"]=math.ceil(math.log2(matr1_n_plain_multiply[k].scale()))
#                    except Exception as error:
#                        ret["errors"] = ret.get("errors", []) + [f"{time.time()}: c(Matriz 1)*Numero {k}–Reescalado–{type(error).__name__}–{error}"]
#                        print(f"SEAL-Error Reescalando: {type(error).__name__}–{error}")
#                    m=process.memory_info().rss
#                    t=time.time()
#                    matr1_n_plain_multiply_dencript[k]=decryptor.decrypt(matr1_n_plain_multiply[k])
#                    ret[f"t_matr1_{k}_plain_multiply_dencript"]=time.time()-t
#                    ret[f"m_matr1_{k}_plain_multiply_dencript"]=process.memory_info().rss-m
#                    m=process.memory_info().rss
#                    t=time.time()
#                    matr1_n_plain_multiply_decode[k]=encoder.decode(matr1_n_plain_multiply_dencript[k])
#                    ret[f"t_matr1_{k}_plain_multiply_decode"]=time.time()-t
#                    ret[f"m_matr1_{k}_plain_multiply_decode"]=process.memory_info().rss-m
#                    ret[f"t_matr1_{k}_plain_multiply_total"]=sum(ret[x] for x in ret.keys() if x.startswith(f"t_matr1_{k}_plain_multiply"))
#                    ret[f"m_matr1_{k}_plain_multiply_total"]=sum(ret[x] for x in ret.keys() if x.startswith(f"m_matr1_{k}_plain_multiply"))
#                    matr1_n_plain_multiply_decode[k]=matr1_n_plain_multiply_decode[k][:filas_vect1*columnas_vect1].reshape(filas_vect1,columnas_vect1)
#                    ret[f"e_matr1_{k}_plain_multiply"]=numpy.mean(abs(numeros[k]*matriz1-matr1_n_plain_multiply_decode[k]))
#                except Exception as error:
#                    ret["errors"] = ret.get("errors", []) + [f"{time.time()}: c(Matriz 1)*Numero {k}–{type(error).__name__}–{error}"]
#                    print(f"SEAL-c(Matriz 1) * Numeros: {type(error).__name__}–{error}")
#
#                # c(Matriz 1) * c(Numeros)
#                try:
#                    m=process.memory_info().rss
#                    t=time.time()
#                    matr1_n_multiply[k]=evaluator.multiply(matr1_encrypt,num_encrypt[k])
#                    ret[f"t_matr1_{k}_multiply"]=time.time()-t
#                    ret[f"m_matr1_{k}_multiply"]=process.memory_info().rss-m
#                    ret[f"b_matr1_{k}_multiply"]=math.ceil(math.log2(matr1_n_multiply[k].scale()))
#                    if relin_keys!=None:
#                        m=process.memory_info().rss
#                        t=time.time()
#                        evaluator.relinearize_inplace(matr1_n_multiply[k], relin_keys)
#                        ret[f"t_matr1_{k}_multiply_realin"]=time.time()-t
#                        ret[f"m_matr1_{k}_multiply_realin"]=process.memory_info().rss-m
#                        ret[f"b_matr1_{k}_multiply_realin"]=math.ceil(math.log2(matr1_n_multiply[k].scale()))
#                    try:
#                        m=process.memory_info().rss
#                        t=time.time()
#                        evaluator.rescale_to_next_inplace(matr1_n_multiply[k])
#                        ret[f"t_matr1_{k}_multiply_rescale"]=time.time()-t
#                        ret[f"m_matr1_{k}_multiply_rescale"]=process.memory_info().rss-m
#                        ret[f"b_matr1_{k}_multiply_rescale"]=math.ceil(math.log2(matr1_n_multiply[k].scale()))
#                    except Exception as error:
#                        ret["errors"] = ret.get("errors", []) + [f"{time.time()}: c(Matriz 1)*c(Numero {k})–Reescalado–{type(error).__name__}–{error}"]
#                        print(f"SEAL-Error Reescalando: {type(error).__name__}–{error}")
#                    m=process.memory_info().rss
#                    t=time.time()
#                    matr1_n_multiply_dencript[k]=decryptor.decrypt(matr1_n_multiply[k])
#                    ret[f"t_matr1_{k}_multiply_dencript"]=time.time()-t
#                    ret[f"m_matr1_{k}_multiply_dencript"]=process.memory_info().rss-m
#                    m=process.memory_info().rss
#                    t=time.time()
#                    matr1_n_multiply_decode[k]=encoder.decode(matr1_n_multiply_dencript[k])
#                    ret[f"t_matr1_{k}_multiply_decode"]=time.time()-t
#                    ret[f"m_matr1_{k}_multiply_decode"]=process.memory_info().rss-m
#                    ret[f"t_matr1_{k}_multiply_total"]=sum(ret[x] for x in ret.keys() if x.startswith(f"t_matr1_{k}_multiply"))
#                    ret[f"m_matr1_{k}_multiply_total"]=sum(ret[x] for x in ret.keys() if x.startswith(f"m_matr1_{k}_multiply"))
#                    matr1_n_multiply_decode[k]=matr1_n_multiply_decode[k][:filas_vect1*columnas_vect1].reshape(filas_vect1,columnas_vect1)
#                    ret[f"e_matr1_{k}_multiply"]=numpy.mean(abs(numeros[k]*matriz1-matr1_n_multiply_decode[k]))
#                except Exception as error:
#                    ret["errors"] = ret.get("errors", []) + [f"{time.time()}: c(Matriz 1)*c(Numero {k})–{type(error).__name__}–{error}"]
#                    print(f"SEAL-c(Matriz 1) * c(Numeros): {type(error).__name__}–{error}")
#        except Exception as error:
#            ret["errors"] = ret.get("errors", []) + [f"{time.time()}: Numeros–{type(error).__name__}–{error}"]
#            print(f"SEAL-Numeros: {type(error).__name__}–{error}")
#
#        if(len(matriz1[0])==len(matriz2[0])):
#            # c(Matriz 1) * Matriz 2
#            try:
#                m=process.memory_info().rss
#                t=time.time()
#                matr1_matr2_plain_multiply=evaluator.multiply_plain(matr1_encrypt, matr2_encode)
#                ret["t_matr1_matr2_plain_multiply"]=time.time()-t
#                ret["m_matr1_matr2_plain_multiply"]=process.memory_info().rss-m
#                ret["b_matr1_matr2_plain_multiply"]=math.ceil(math.log2(matr1_matr2_plain_multiply.scale()))
#                try:
#                    m=process.memory_info().rss
#                    t=time.time()
#                    evaluator.rescale_to_next_inplace(matr1_matr2_plain_multiply)
#                    ret["t_matr1_matr2_plain_multiply_rescale"]=time.time()-t
#                    ret["m_matr1_matr2_plain_multiply_rescale"]=process.memory_info().rss-m
#                    ret["b_matr1_matr2_plain_multiply_rescale"]=math.ceil(math.log2(matr1_matr2_plain_multiply.scale()))
#                except Exception as error:
#                    ret["errors"] = ret.get("errors", []) + [f"{time.time()}: c(Matriz 1)*Matriz 2–Reescalando–{type(error).__name__}–{error}"]
#                    print(f"SEAL-Error Reescalando: {type(error).__name__}–{error}")
#                m=process.memory_info().rss
#                t=time.time()
#                matr1_matr2_plain_multiply_dencript=decryptor.decrypt(matr1_matr2_plain_multiply)
#                ret["t_matr1_matr2_plain_multiply_dencript"]=time.time()-t
#                ret["m_matr1_matr2_plain_multiply_dencript"]=process.memory_info().rss-m#
#                m=process.memory_info().rss
#                t=time.time()
#                matr1_matr2_plain_multiply_decode=encoder.decode(matr1_matr2_plain_multiply_dencript)
#                ret["t_matr1_matr2_plain_multiply_decode"]=time.time()-t
#                ret["m_matr1_matr2_plain_multiply_decode"]=process.memory_info().rss-m
#                ret["t_matr1_matr2_plain_multiply_total"]=sum(ret[x] for x in ret.keys() if x.startswith("t_matr1_matr2_plain_multiply"))
#                ret["m_matr1_matr2_plain_multiply_total"]=sum(ret[x] for x in ret.keys() if x.startswith("m_matr1_matr2_plain_multiply"))
#                f_multiply = min(filas_vect1,filas_vect2)
#                c_multiply = (int)(min(columnas_vect1,columnas_vect2))
#                matr1_matr2_plain_multiply_decode=matr1_matr2_plain_multiply_decode[:f_multiply*c_multiply].reshape(f_multiply,c_multiply)#
#                ret["e_matr1_matr2_plain_multiply"]=numpy.mean(abs(matriz1[:min(filas_vect1,filas_vect2)]*matriz2[:min(filas_vect1,filas_vect2)]-matr1_matr2_plain_multiply_decode))
#            except Exception as error:
#                ret["errors"] = ret.get("errors", []) + [f"{time.time()}: c(Matriz 1)*Matriz 2–{type(error).__name__}–{error}"]
#                print(f"SEAL-c(Matriz 1) * Matriz 2: {type(error).__name__}–{error}")
#            
#            # c(Matriz 1) * c(Matriz 2)
#            try:
#                m=process.memory_info().rss
#                t=time.time()
#                matr1_matr2_multiply=evaluator.multiply(matr1_encrypt, matr2_encrypt)
#                ret["t_matr1_matr2_multiply"]=time.time()-t
#                ret["m_matr1_matr2_multiply"]=process.memory_info().rss-m
#                ret["b_matr1_matr2_multiply"]=math.ceil(math.log2(matr1_matr2_multiply.scale()))
#                if relin_keys!=None:
#                    m=process.memory_info().rss
#                    t=time.time()
#                    evaluator.relinearize_inplace(matr1_matr2_multiply, relin_keys)
#                    ret["t_matr1_matr2_multiply_realin"]=time.time()-t
#                    ret["m_matr1_matr2_multiply_realin"]=process.memory_info().rss-m
#                    ret["b_matr1_matr2_multiply_realin"]=math.ceil(math.log2(matr1_matr2_multiply.scale()))
#                try:
#                    m=process.memory_info().rss
#                    t=time.time()
#                    evaluator.rescale_to_next_inplace(matr1_matr2_multiply)
#                    ret["t_matr1_matr2_multiply_rescale"]=time.time()-t
#                    ret["m_matr1_matr2_multiply_rescale"]=process.memory_info().rss-m
#                    ret["b_matr1_matr2_multiply_rescale"]=math.ceil(math.log2(matr1_matr2_multiply.scale()))
#                except Exception as error:
#                    ret["errors"] = ret.get("errors", []) + [f"{time.time()}: c(Matriz 1)*c(Matriz 2)–Reescalando–{type(error).__name__}–{error}"]
#                    print(f"SEAL-Error Reescalando: {type(error).__name__}–{error}")
#                m=process.memory_info().rss
#                t=time.time()
#                matr1_matr2_multiply_dencript=decryptor.decrypt(matr1_matr2_multiply)
#                ret["t_matr1_matr2_multiply_dencript"]=time.time()-t
#                ret["m_matr1_matr2_multiply_dencript"]=process.memory_info().rss-m
#                m=process.memory_info().rss
#                t=time.time()
#                matr1_matr2_multiply_decode=encoder.decode(matr1_matr2_multiply_dencript)
#                ret["t_matr1_matr2_multiply_decode"]=time.time()-t
#                ret["m_matr1_matr2_multiply_decode"]=process.memory_info().rss-m
#                ret["t_matr1_matr2_multiply_total"]=sum(ret[x] for x in ret.keys() if x.startswith(f"t_matr1_matr2_multiply"))
#                ret["m_matr1_matr2_multiply_total"]=sum(ret[x] for x in ret.keys() if x.startswith(f"m_matr1_matr2_multiply"))
#                f_multiply = min(filas_vect1,filas_vect2)
#                c_multiply = min(columnas_vect1,columnas_vect2)
#                matr1_matr2_multiply_decode=matr1_matr2_multiply_decode[:f_multiply*c_multiply].reshape(f_multiply,c_multiply)
#                ret["e_matr1_matr2_multiply"]=numpy.mean(abs(matriz1[:min(filas_vect1,filas_vect2)]*matriz2[:min(filas_vect1,filas_vect2)]-matr1_matr2_multiply_decode))
#            except Exception as error:
#                ret["errors"] = ret.get("errors", []) + [f"{time.time()}: c(Matriz 1)*c(Matriz 2)–{type(error).__name__}–{error}"]
#                print(f"SEAL-c(Matriz 1) * c(Matriz 2): {type(error).__name__}–{error}")
        
        # c(Matriz 1) matrix multiply c(Matriz 2)
        try:
            if not (filas_vect1 == columnas_vect1 == filas_vect2 == columnas_vect2):
                raise Exception(f"Para ejecutar HEMatMult es necesario que las matrices sean cuadaradas del mismo tamaño, actuales: A({ret['tam_matr1']}), B({ret['tam_matr2']})")
            if relin_keys!=None and gal_keys!=None:
                m=process.memory_info().rss
                t=time.time()
                scale=pow(2.0, params["pre"])
                filas_vect1 = len(matriz1)
                if "cmp" in params:
                    coeff_modulus = list(map(lambda x: Modulus(x), params["cmp"]))
                elif "cma" in params:
                    coeff_modulus = CoeffModulus.Create(params["pmd"], params["cma"])
                elif "cm" in params:
                    n_coeff_modulus=math.ceil(params["cm"]/MAX_BITS_COEFF_MODULUS)
                    v_coeff_modulus=[math.floor(params["cm"]/n_coeff_modulus)] * n_coeff_modulus
                    for i in range(n_coeff_modulus-math.floor(params["cm"]%n_coeff_modulus), len(v_coeff_modulus)):
                        v_coeff_modulus[i] += 1
                    coeff_modulus = CoeffModulus.Create(params["pmd"], v_coeff_modulus)
                else:
                    if 1024 <= params["pmd"] <= 32768:
                        coeff_modulus = CoeffModulus.BFVDefault(params["pmd"], encriptacion_types.get(params["enc"], sec_level_type.tc128))
                        #coeff_modulus = CoeffModulus.Create(params["pmd"], [60, 60, 60, 60, 60, 60, 60, 60, 60])
                    elif params["pmd"] > 32768:
                        b_coeff_modulus=params["pmd"]/32768*sum(map(lambda x: x.bit_count(), CoeffModulus.BFVDefault(32768, encriptacion_types.get(params["enc"], sec_level_type.tc128))))
                        n_coeff_modulus=math.ceil(b_coeff_modulus/MAX_BITS_COEFF_MODULUS)
                        v_coeff_modulus=[math.floor(b_coeff_modulus/n_coeff_modulus)] * n_coeff_modulus
                        for i in range(n_coeff_modulus-math.floor(b_coeff_modulus%n_coeff_modulus), len(v_coeff_modulus)):
                            v_coeff_modulus[i] += 1
                        coeff_modulus = CoeffModulus.Create(params["pmd"], v_coeff_modulus)
                    elif params["pmd"] < 1024:
                        b_coeff_modulus=params["pmd"]/1024*sum(map(lambda x: x.bit_count(), CoeffModulus.BFVDefault(1024, encriptacion_types.get(params["enc"], sec_level_type.tc128))))
                        n_coeff_modulus=math.ceil(b_coeff_modulus/MAX_BITS_COEFF_MODULUS)
                        v_coeff_modulus=[math.floor(b_coeff_modulus/n_coeff_modulus)] * n_coeff_modulus
                        for i in range(n_coeff_modulus-math.floor(b_coeff_modulus%n_coeff_modulus), len(v_coeff_modulus)):
                            v_coeff_modulus[i] += 1
                        coeff_modulus = CoeffModulus.Create(params["pmd"], v_coeff_modulus)
                parms = EncryptionParameters(scheme_type.ckks)
                parms.set_poly_modulus_degree(params["pmd"])
                parms.set_coeff_modulus(coeff_modulus)
                context = SEALContext(parms, True, encriptacion_types.get(params["enc"], sec_level_type.tc128) if params["pmd"] <= 32768 else sec_level_type.none)
                keygen = KeyGenerator(context)
                secret_key = keygen.secret_key()
                public_key = keygen.create_public_key()
                relin_keys = keygen.create_relin_keys()
                gal_keys = keygen.create_galois_keys()
                encryptor = Encryptor(context, public_key)
                evaluator = Evaluator(context)
                decryptor = Decryptor(context, secret_key)#
                encoder = CKKSEncoder(context)
                
                if math.log2(filas_vect1) == int(math.log2(filas_vect1)):
                    vect1 = list(matriz1.reshape(-1)) * math.floor(params["pmd"]/2/(filas_vect1**2))
                    vect2 = list(matriz2.reshape(-1)) * math.floor(params["pmd"]/2/(filas_vect1**2))
                    lado_HEMM = filas_vect1
                else:
                    lado_HEMM = 2**math.ceil(math.log2(filas_vect1))
                    def matr_HEMM(matr):
                        return (sum([list(matr[i]) + [0] * (lado_HEMM-filas_vect1) for i in range(filas_vect1)], []) + [0] * lado_HEMM * (lado_HEMM-filas_vect1)) * math.floor(params["pmd"]/2/(lado_HEMM**2))
                    vect1 = matr_HEMM(matriz1)
                    vect2 = matr_HEMM(matriz2)
                matr1_encrypt_HEMM = encryptor.encrypt(encoder.encode(vect1, scale))
                matr2_encrypt_HEMM = encryptor.encrypt(encoder.encode(vect2, scale))

                heMatMult = HEMatMult(int(context_data.parms().poly_modulus_degree()/2), scale, scale, list(map(lambda x: x.bit_count(), context_data.parms().coeff_modulus())), evaluator, encoder, relin_keys, gal_keys)
                #heMatMult = HEMatMult(int(context_data.parms().poly_modulus_degree()/2), scale, scale, list(map(lambda x: x.bit_count(), context_data.parms().coeff_modulus())), evaluator, encoder, relin_keys, gal_keys, decryptor)

                matr1_matr2_matrix_multiply = heMatMult.multiplica(matr1_encrypt_HEMM, matr2_encrypt_HEMM, lado_HEMM)
                #matr1_matr2_matrix_multiply = heMatMult.multiplica_coment(matr1_encrypt_HEMM, matr2_encrypt_HEMM, lado_HEMM, vect1, vect2)
                
                matr1_matr2_matrix_multiply_decode=encoder.decode(decryptor.decrypt(matr1_matr2_matrix_multiply))[:lado_HEMM**2].reshape(lado_HEMM,lado_HEMM)[:filas_vect1, :filas_vect1]
                ret["t_matr1_matr2_matrix_multiply_total"]=time.time()-t
                ret["m_matr1_matr2_matrix_multiply_total"]=process.memory_info().rss-m
                matMultEsp = numpy.matmul(matriz1, matriz2)
                ret["e_matr1_matr2_matrix_multiply"]=numpy.mean(abs(matr1_matr2_matrix_multiply_decode-matMultEsp)/matMultEsp)
                
                
                m=process.memory_info().rss
                t=time.time()
                if math.log2(filas_vect1) == int(math.log2(filas_vect1)):
                    vect1 = list(matriz1.reshape(-1)) * math.floor(params["pmd"]/2/(filas_vect1**2))
                    vect2 = list(matriz2.reshape(-1)) * math.floor(params["pmd"]/2/(filas_vect1**2))
                    lado_HEMM = filas_vect1
                else:
                    lado_HEMM = 2**math.ceil(math.log2(filas_vect1))
                    def matr_HEMM(matr):
                        return (sum([list(matr[i]) + [0] * (lado_HEMM-filas_vect1) for i in range(filas_vect1)], []) + [0] * lado_HEMM * (lado_HEMM-filas_vect1)) * math.floor(params["pmd"]/2/(lado_HEMM**2))
                    vect1 = matr_HEMM(matriz1)
                    vect2 = matr_HEMM(matriz2)
                ret["t_matr1_matr2_matrix_multiply_prep_inputs"]=time.time()-t
                ret["m_matr1_matr2_matrix_multiply_prep_inputs"]=process.memory_info().rss-m
                m=process.memory_info().rss
                t=time.time()
                matr1_encrypt_HEMM = encryptor.encrypt(encoder.encode(vect1, scale))
                matr2_encrypt_HEMM = encryptor.encrypt(encoder.encode(vect2, scale))
                ret["t_matr1_matr2_matrix_multiply_encrypt"]=time.time()-t
                ret["m_matr1_matr2_matrix_multiply_encrypt"]=process.memory_info().rss-m
                m=process.memory_info().rss
                t=time.time()
                heMatMult = HEMatMult(int(context_data.parms().poly_modulus_degree()/2), scale, scale, list(map(lambda x: x.bit_count(), context_data.parms().coeff_modulus())), evaluator, encoder, relin_keys, gal_keys)
                #heMatMult = HEMatMult(int(context_data.parms().poly_modulus_degree()/2), scale, 2**round(params["pre"]/2), list(map(lambda x: x.bit_count(), context_data.parms().coeff_modulus())), evaluator, encoder, relin_keys, gal_keys, decryptor)
                ret["t_matr1_matr2_matrix_multiply_prep"]=time.time()-t
                ret["m_matr1_matr2_matrix_multiply_prep"]=process.memory_info().rss-m
                m=process.memory_info().rss
                t=time.time()
                #matr1_matr2_matrix_multiply = heMatMult.multiplica(matr1_encrypt, matr2_encrypt, filas_vect1)
                matr1_matr2_matrix_multiply = heMatMult.multiplica(matr1_encrypt_HEMM, matr2_encrypt_HEMM, lado_HEMM)
                #matr1_matr2_matrix_multiply = heMatMult.multiplica_coment(matr1_encrypt, matr2_encrypt, filas_vect1, vect1, vect2)
                #matr1_matr2_matrix_multiply = heMatMult(matr1_encrypt, matr2_encrypt, filas_vect1, int(context_data.parms().poly_modulus_degree()/2), scale, 2**20, list(map(lambda x: x.bit_count(), context_data.parms().coeff_modulus())), evaluator, encoder, relin_keys, gal_keys, decryptor, vect1, vect2)
                ret["t_matr1_matr2_matrix_multiply"]=time.time()-t
                ret["m_matr1_matr2_matrix_multiply"]=process.memory_info().rss-m
                ret["b_matr1_matr2_matrix_multiply"]=math.ceil(math.log2(matr1_matr2_matrix_multiply.scale()))
                m=process.memory_info().rss
                t=time.time()
                matr1_matr2_matrix_multiply_decode=encoder.decode(decryptor.decrypt(matr1_matr2_matrix_multiply))
                ret["t_matr1_matr2_matrix_multiply_decrypt"]=time.time()-t
                ret["m_matr1_matr2_matrix_multiply_decrypt"]=process.memory_info().rss-m
                m=process.memory_info().rss
                t=time.time()
                matr1_matr2_matrix_multiply_decode=matr1_matr2_matrix_multiply_decode[:lado_HEMM**2].reshape(lado_HEMM,lado_HEMM)[:filas_vect1, :filas_vect1]
                ret["t_matr1_matr2_matrix_multiply_prep_outputs"]=time.time()-t
                ret["m_matr1_matr2_matrix_multiply_prep_outputs"]=process.memory_info().rss-m
                ret["e_matr1_matr2_matrix_multiply"]=(ret["e_matr1_matr2_matrix_multiply"]+numpy.mean(abs(matr1_matr2_matrix_multiply_decode-matMultEsp)/matMultEsp))/2
                
                
                #err =abs(numpy.matmul(matriz1, matriz2)-matr1_matr2_matrix_multiply_decode)
                #vAB = numpy.array(vAB[:f_multiply*c_multiply]).reshape(f_multiply,c_multiply)
                #err2 =abs(numpy.matmul(matriz1, matriz2)-vAB)
                
                #print("Matriz 1:\n{0}\nMatriz 2:\n{1}\nheMatMult:\n{2}\nMatrix Multiplic:\n{3}\nError:\n\t±{4}\n\t>{5}\n\t<{6}".format(matriz1, matriz2, matr1_matr2_matrix_multiply_decode, numpy.matmul(matriz1, matriz2), ret["e_matr1_matr2_matrix_multiply"], numpy.max(err), numpy.min(err)))
                #print("Matriz 1:\n{0}\nMatriz 2:\n{1}\nheMatMult:\n{2}\nMatrix Multiplic:\n{3}\nError:\n\t±{4}\n\t>{5}\n\t<{6}\nTestigo:\n{7}\nError:\n\t±{8}\n\t>{9}\n\t<{10}".format(matriz1, matriz2, matr1_matr2_matrix_multiply_decode, numpy.matmul(matriz1, matriz2), ret["e_matr1_matr2_matrix_multiply"], numpy.max(err), numpy.min(err), vAB, numpy.mean(err2), numpy.max(err2), numpy.min(err2)))
                #print(f"t: {ret['t_matr1_matr2_matrix_multiply_total']}")
                #try:
                #    print("heMatMult_list:\n{0}".format(numpy.array(heMatMult_list(sum_list([0]*(int(context_data.parms().poly_modulus_degree()/2)), vect1), sum_list([0]*(int(context_data.parms().poly_modulus_degree()/2)), vect2), filas_vect1)[:filas_vect1**2]).reshape(filas_vect1,filas_vect1)))
                #except Exception as error:
                #    print(error)
                #try:
                #    print("heMatMult_list:\n{0}".format(numpy.array(heMatMult_list(list(vect1), list(vect2), filas_vect1)[:filas_vect1**2]).reshape(filas_vect1,filas_vect1)))
                #except Exception as error:
                #    print(error)
        except Exception as error:
            ret["errors"] = ret.get("errors", []) + [f"{time.time()}: c(Matriz 1)*c(Matriz 2)–{type(error).__name__}–{error}"]
            print(f"SEAL-c(Matriz 1) matrix multiply c(Matriz 2): {type(error).__name__}–{error}")
        
        # c(Matriz 1) general matrix multiply c(Matriz 2)
        #try:
        #    if columnas_vect1 != filas_vect2:
        #        raise Exception(f"Para poder ejecutar la HEGMM es necesario que la matiz 1 tenga el mismo nº de columnas que la matriz 2 nº de filas, actuales: A({ret['tam_matr1']}), B({ret['tam_matr2']})")
        #    if relin_keys!=None and gal_keys!=None:
        #        m=process.memory_info().rss
        #        t=time.time()
        #        #hEGMM = HEGMM(filas_vect1, columnas_vect1, columnas_vect2, int(context_data.parms().poly_modulus_degree()/2), scale, 2**round(params["pre"]/2), list(map(lambda x: x.bit_count(), context_data.parms().coeff_modulus())), evaluator, encoder, relin_keys, gal_keys)
        #        hEGMM = HEGMM(filas_vect1, columnas_vect1, columnas_vect2, int(context_data.parms().poly_modulus_degree()/2), scale, 2**round(params["pre"]/2), list(map(lambda x: x.bit_count(), context_data.parms().coeff_modulus())), evaluator, encoder, relin_keys, gal_keys, decryptor)
        #        ret["t_matr1_matr2_general_matrix_multiply_prep"]=time.time()-t
        #        ret["m_matr1_matr2_general_matrix_multiply_prep"]=process.memory_info().rss-m
        #        if hEGMM.m*hEGMM.n > ret["poly_modulus_degree"]/2:
        #            raise Exception(f"Para poder ejecutar la HEGMM es necesario un tamaño de vector de ({hEGMM.m}x{hEGMM.n})={hEGMM.m*hEGMM.n}, tam_vector: {ret['poly_modulus_degree']/2}")
        #        print(f"m: {hEGMM.m}, l: {hEGMM.l}, n:{hEGMM.n}")
        #        pm=process.memory_info().rss
        #        t=time.time()
        #        prepMatr1 = hEGMM.oA(vect1[:filas_vect1*columnas_vect1])
        #        prepMatr2 = hEGMM.tB(vect2[:filas_vect2*columnas_vect2])
        #        ret["t_matr1_matr2_general_matrix_multiply_prep_matr"]=time.time()-t
        #        ret["m_matr1_matr2_general_matrix_multiply_prep_matr"]=process.memory_info().rss-m
        #        print(f"prepMatr1({len(prepMatr1)}): {prepMatr1}")
        #        print(f"prepMatr2({len(prepMatr2)}): {prepMatr2}")
        #        m=process.memory_info().rss
        #        t=time.time()
        #        prepMatr1_encrypt = encryptor.encrypt(encoder.encode(prepMatr1, scale))
        #        prepMatr2_encrypt = encryptor.encrypt(encoder.encode(prepMatr2, scale))
        #        ret["t_matr1_matr2_general_matrix_multiply_encr_prep_matr"]=time.time()-t
        #        ret["m_matr1_matr2_general_matrix_multiply_encr_prep_matr"]=process.memory_info().rss-m
        #        m=process.memory_info().rss
        #        t=time.time()
        #        #matr1_matr2_general_matrix_multiply = hEGMM.multiplica(prepMatr1_encrypt, prepMatr2_encrypt)
        #        matr1_matr2_general_matrix_multiply = hEGMM.multiplica_coment(prepMatr1_encrypt, prepMatr2_encrypt, prepMatr1, prepMatr2)
        #        ret["t_matr1_matr2_general_matrix_multiply"]=time.time()-t
        #        ret["m_matr1_matr2_general_matrix_multiply"]=process.memory_info().rss-m
        #        ret["b_matr1_matr2_general_matrix_multiply"]=math.ceil(math.log2(matr1_matr2_general_matrix_multiply.scale()))
        #        m=process.memory_info().rss
        #        t=time.time()
        #        matr1_matr2_general_matrix_multiply_dencript=decryptor.decrypt(matr1_matr2_general_matrix_multiply)
        #        ret["t_matr1_matr2_general_matrix_multiply_dencript"]=time.time()-t
        #        ret["m_matr1_matr2_general_matrix_multiply_dencript"]=process.memory_info().rss-m
        #        m=process.memory_info().rss
        #        t=time.time()
        #        matr1_matr2_general_matrix_multiply_decode=encoder.decode(matr1_matr2_general_matrix_multiply_dencript)
        #        ret["t_matr1_matr2_general_matrix_multiply_decode"]=time.time()-t
        #        ret["m_matr1_matr2_general_matrix_multiply_decode"]=process.memory_info().rss-m
        #        ret["t_matr1_matr2_general_matrix_multiply_total"]=sum(ret[x] for x in ret.keys() if x.startswith("t_matr1_matr2_general_matrix_multiply"))
        #        ret["m_matr1_matr2_general_matrix_multiply_total"]=sum(ret[x] for x in ret.keys() if x.startswith("m_matr1_matr2_general_matrix_multiply"))
        #        f_multiply, c_multiply = hEGMM.getTamMatC()
        #        matr1_matr2_general_matrix_multiply_decode=numpy.array(hEGMM.arregloC(matr1_matr2_general_matrix_multiply_decode[:hEGMM.M()*hEGMM.N()])).reshape(f_multiply,c_multiply)
        #        ret["e_matr1_matr2_general_matrix_multiply"]=numpy.mean(abs(numpy.matmul(matriz1, matriz2)-matr1_matr2_general_matrix_multiply_decode))
        #        err =abs(numpy.matmul(matriz1, matriz2)-matr1_matr2_general_matrix_multiply_decode)
        #        #vAB = numpy.array(vAB[:f_multiply*c_multiply]).reshape(f_multiply,c_multiply)
        #        #err2 =abs(numpy.matmul(matriz1, matriz2)-vAB)
        #        
        #        #print("Matriz 1:\n{0}\nMatriz 2:\n{1}\nHEGMM:\n{2}\nMatrix Multiplic:\n{3}\nError:\n\t±{4}\n\t>{5}\n\t<{6}".format(matriz1, matriz2, matr1_matr2_general_matrix_multiply_decode, numpy.matmul(matriz1, matriz2), ret["e_matr1_matr2_general_matrix_multiply"], numpy.max(err), numpy.min(err)))
        #        #print(f"t: {ret['t_matr1_matr2_general_matrix_multiply_total']}")
        #        
        #except Exception as error:
        #    ret["errors"] = ret.get("errors", []) + [f"{time.time()}: c(Matriz 1)*c(Matriz 2)–{type(error).__name__}–{error}"]
        #    print(f"SEAL-c(Matriz 1) general matrix multiply c(Matriz 2): {type(error).__name__}–{error}")

        if not params.get("silence", False):
            for k in ret.keys():
                print(f"{k}: {ret[k]}")
    
    except Exception as error:
        ret["errors"] = ret.get("errors", []) + [f"{time.time()}: {type(error).__name__}–{error}"]
        print(f"ERROR en SEAL: {type(error).__name__}–{error}")
    finally:
        ret["fin"] = time.time()

        return ret
