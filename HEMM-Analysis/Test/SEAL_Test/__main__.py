import os
import psutil
import time
import json
import math
import numpy

from scipy import optimize

from seal import *
from heMatMult import HEMatMult

MAX_BITS_COEFF_MODULUS = 60
MIN_BITS_COEFF_MODULUS = 15

encriptacion_types = {
    128: sec_level_type.tc128,
    192: sec_level_type.tc192,
    256: sec_level_type.tc256,
    'None': sec_level_type.none
}

def obtener_in_out(tam_lado, distribucion_in, prec_in):
    tam_list = tam_lado**2
    list1 = distribucion_in(tam_list)
    list2 = distribucion_in(tam_list)
    if prec_in != None:
        list1 = numpy.around(list1, -1*prec_in)
        list2 = numpy.around(list2, -1*prec_in)
    matr1 = list1.reshape(tam_lado, tam_lado)
    matr2 = list2.reshape(tam_lado, tam_lado)
    prodMatr = numpy.matmul(matr1, matr2)
    #print(f"Matriz 1:\n{matr1}\nMatriz 2:\n{matr2}\nProducto de matrices:\n{prodMatr}\n")
    return list(list1), list(list2), prodMatr, matr1, matr2

def poly_modulus_degrees(encriptacion, tam_list, pmd_estandar):
    pmd = max(16384 if encriptacion > 192 else 8192, tam_list*2)
    while True:
        if pmd > 32768:
            if pmd_estandar:
                return
            else:
                print(f"WARNING: El poly_modulus_degree que se está evaluando es de {pmd}, no se encuentra dentro de los estandarizados, este proceso puede consumir mucho tiempo y memoria.")
        yield pmd
        pmd *= 2

def creaContextKeys(pmd, coeff_modulus, enc, repeticiones, tam_lado, distribucion_in, prec_in, prec_out, pre_rec = set(), rapido = False, process = psutil.Process(os.getpid())):
    try:
        testeos = {"rapido": True} if rapido else {}
        cml = list(map(lambda x: x.bit_count(), coeff_modulus))
        m = process.memory_info().rss
        t = time.time()
        parms = EncryptionParameters(scheme_type.ckks)
        parms.set_poly_modulus_degree(pmd)
        parms.set_coeff_modulus(coeff_modulus)
        testeos["t_params"] = time.time()-t
        testeos["m_params"] = process.memory_info().rss-m
        m = process.memory_info().rss
        t = time.time()
        context = SEALContext(parms, True, encriptacion_types.get(enc, sec_level_type.tc128) if pmd <= 32768 else sec_level_type.none)
        testeos["t_context"] = time.time()-t
        testeos["m_context"] = process.memory_info().rss-m
        m = process.memory_info().rss
        t = time.time()
        evaluator = Evaluator(context)
        testeos["t_evaluator"] = time.time()-t
        testeos["m_evaluator"] = process.memory_info().rss-m
        m = process.memory_info().rss
        t = time.time()
        encoder = CKKSEncoder(context)
        testeos["t_encoder"] = time.time()-t
        testeos["m_encoder"] = process.memory_info().rss-m
        m = process.memory_info().rss
        t = time.time()
        keygen = KeyGenerator(context)
        testeos["t_keygen"] = time.time()-t
        testeos["m_keygen"] = process.memory_info().rss-m
        m = process.memory_info().rss
        t = time.time()
        encryptor = Encryptor(context, keygen.create_public_key())
        testeos["t_encryptor"] = time.time()-t
        testeos["m_encryptor"] = process.memory_info().rss-m
        m = process.memory_info().rss
        t = time.time()
        decrypytor = Decryptor(context, keygen.secret_key())
        testeos["t_decrypytor"] = time.time()-t
        testeos["m_decrypytor"] = process.memory_info().rss-m
        m = process.memory_info().rss
        t = time.time()
        relin_keys = keygen.create_relin_keys()
        testeos["t_relin_keys"] = time.time()-t
        testeos["m_relin_keys"] = process.memory_info().rss-m
        m = process.memory_info().rss
        t = time.time()
        galois_keys = keygen.create_galois_keys()
        testeos["t_galois_keys"] = time.time()-t
        testeos["m_galois_keys"] = process.memory_info().rss-m
        testeos["t_total"] = sum(testeos[x] for x in testeos.keys() if x.startswith("t_"))
        testeos["m_total"] = sum(testeos[x] for x in testeos.keys() if x.startswith("m_"))
        tam_vect = int(pmd/2)

        precisiones = {}
        pre_max = sum(cml[:-1])-2
        pre_min = 1

        pre_valid = set()               # No dan error
        pre_pass = set()                # Pasan la restricción prec_out
        pre_best = (None, float("inf")) # La precisión con mejor resultado

        def pre_min_pass():
            if len(pre_pass) > 0 and len(pre_valid) > 0:
                min_pass= min(pre_pass)
                ret = list(filter(lambda pre: pre < min_pass, pre_valid))
                if len(ret) > 0:
                    return max(ret)
            return pre_min
            
        def pre_max_pass():
            if len(pre_pass) > 0 and len(pre_valid) > 0:
                max_pass= max(pre_pass)
                ret = list(filter(lambda pre: pre > max_pass, pre_valid))
                if len(ret) > 0:
                    return min(ret)
            return pre_max

        testPrec_l = lambda pre: testPrec(pre, prec_in, prec_out, distribucion_in, tam_lado, tam_vect, cml,
                                        evaluator, encoder, encryptor, decrypytor, relin_keys, galois_keys,
                                        1 if rapido else repeticiones, process)

        if len(list(pre_rec)) > 0:
            pre_rec = list(filter(lambda pre: pre_min < pre < pre_max, pre_rec))
            pre_rec.sort()

            for pre in pre_rec:
                if (rapido and len(pre_pass) > 0):
                    break
                if pre_min < pre < pre_max and pre not in precisiones:
                    precisiones[pre], pre_max, pre_min, pre_valid, pre_pass, pre_best = optTestPrec(pre, testPrec_l, prec_out, pre_max, pre_min, pre_valid, pre_pass, pre_best)
                    
        while True:
            if rapido and len(pre_pass) > 0:
                break
            pre = int((pre_max_pass()-pre_min_pass())/2+pre_min_pass())
            if pre in precisiones:
                for i in range(5):
                    pre = int((pre_max_pass()-pre_min_pass())/2**(i+1)+pre_min_pass())
                    if pre not in precisiones:
                        break
            if pre in precisiones:
                for i in range(5):
                    pre = int((pre_max_pass()-pre_min_pass())*(2**(i+1)-1)/2**(i+1)+pre_min_pass())
                    if pre not in precisiones:
                        break
            if pre_min < pre < pre_max and pre not in precisiones:
                precisiones[pre], pre_max, pre_min, pre_valid, pre_pass, pre_best = optTestPrec(pre, testPrec_l, prec_out, pre_max, pre_min, pre_valid, pre_pass, pre_best)
            else:
                break

        for pre in filter(lambda pre: pre not in precisiones, range(pre_min_pass()+1,pre_max_pass())):
            if rapido and len(pre_pass) > 0:
                break
            if pre_min < pre < pre_max and pre not in precisiones:
                precisiones[pre], pre_max, pre_min, pre_valid, pre_pass, pre_best = optTestPrec(pre, testPrec_l, prec_out, pre_max, pre_min, pre_valid, pre_pass, pre_best)
        
        testeos["precisiones"] = precisiones
        testeos["precisiones_validas"] = list(pre_valid)
        testeos["precisiones_cumplen_requisitos"] = list(pre_pass)
        testeos["mejor_precision"] = pre_best[0]

        #print(f"pre_valid: {pre_valid}")
        #print(f"pre_pass: {pre_pass}")
        #print(f"pre_best: {pre_best}")
        #print(f"pre_min: {pre_min}")
        #print(f"pre_max: {pre_max}")
        #print(f"pre_min_pass: {pre_min_pass()}")
        #print(f"pre_max_pass: {pre_max_pass()}")
    except Exception as error:
        testeos["error"] = f"[{time.time()}]{type(error).__name__}: {error}"
    finally:
        return testeos

def optTestPrec(pre, testPrec, prec_out, pre_max, pre_min, pre_valid, pre_pass, pre_best):
    ret = testPrec(pre)
    if "error" not in ret and ret["prec_out_enc_dec"]:
        pre_valid.add(pre)
        if ret.get("diff_max", float("inf")) < 10**prec_out:
            pre_pass.add(pre)
        if ret.get("diff_max", float("inf")) < pre_best[1]:
            pre_best = (pre, ret.get("diff_max", float("inf")))
    if not ret["prec_out_enc_dec"]:
        pre_min = max(pre_min, pre)
    if "error" in ret and type(ret["error"]) == list and \
    any(map(lambda e: ret["error"][0].endswith(e), \
            ["ValueError: scale out of bounds", "ValueError: end of modulus switching chain reached", "Exception: No se puede continuar la multiplicación sin descender de la escala mínima."])):
        pre_max = min(pre_max, pre)
    return ret, pre_max, pre_min, pre_valid, pre_pass, pre_best

def testPrec(pre, prec_in, prec_out, distribucion_in, tam_lado, tam_vect, cml, evaluator, encoder, encryptor, decrypytor, relin_keys, galois_keys, repeticiones, process=psutil.Process(os.getpid())):
    test = {"ini_time": time.time()}
    datoTest = [10.0**(prec_in if prec_in != None else prec_out-1)]*tam_vect
    scale = 2**pre
    test["max_dif_enc_dec"] = max(abs(encoder.decode(decrypytor.decrypt(encryptor.encrypt(encoder.encode(datoTest, scale)))) - datoTest))
    test["prec_out_enc_dec"] = bool(test["max_dif_enc_dec"] < 10.0**prec_out)
    if test["prec_out_enc_dec"]:
        #for j in map(lambda x: x/10, reversed(range(5, 11))):
        #    test["prec_min"] = round(pre*j)
        print(f"\t\tprecision: {pre}")#, min: {test['prec_min']}, j: {j}")
        heMatMult = HEMatMult(tam_vect, scale, scale, cml, evaluator, encoder, relin_keys, galois_keys)
        test["tests"] = []
        for i in range(repeticiones):
            tes = {"ini_time": time.time()}
            list1, list2, prodMatr, _, _ = obtener_in_out(tam_lado, distribucion_in, prec_in)
            list1 *= int(tam_vect/len(list1))
            list2 *= int(tam_vect/len(list2))
            m = process.memory_info().rss
            t = time.time()
            list1_enc = encryptor.encrypt(encoder.encode(list1, scale))
            tes["t_enc_list1"] = time.time()-t
            tes["m_enc_list1"] = process.memory_info().rss-m
            m = process.memory_info().rss
            t = time.time()
            list2_enc = encryptor.encrypt(encoder.encode(list2, scale))
            tes["t_enc_list2"] = time.time()-t
            tes["m_enc_list2"] = process.memory_info().rss-m
            try:
                m = process.memory_info().rss
                t = time.time()
                prodMatr_enc = heMatMult.multiplica(list1_enc, list2_enc, tam_lado)
                tes["t_multiplica"] = time.time()-t
                tes["m_multiplica"] = process.memory_info().rss-m
                m = process.memory_info().rss
                t = time.time()
                prodMatr_he = encoder.decode(decrypytor.decrypt(prodMatr_enc))[:tam_lado**2]
                tes["t_decry"] = time.time()-t
                tes["m_decry"] = process.memory_info().rss-m
                tes["t_total"] = sum([tes[k] for k in tes.keys() if k.startswith("t_")])
                tes["m_total"] = sum([tes[k] for k in tes.keys() if k.startswith("m_")])
                prodMatr_he = numpy.array(prodMatr_he).reshape(tam_lado, tam_lado)
                diff_prodMatr = abs(prodMatr - prodMatr_he)
                tes["diff_mean"] = numpy.mean(diff_prodMatr)
                tes["diff_min"] = numpy.min(diff_prodMatr)
                tes["diff_max"] = numpy.max(diff_prodMatr)
                tes["prec_out"] = bool(tes["diff_max"] < 10**prec_out)
                tes["fin_time"] = time.time()
                test["tests"].append(tes)
            except Exception as error:
                tes["error"] = f"[{time.time()}]{type(error).__name__}: {error}"
                tes["fin_time"] = time.time()
                test["tests"].append(tes)
                break
            #if "error" not in tes:
            #    break
        try:
            if any(map(lambda tes: "error" in tes, test["tests"])):
                test["error"] = [tes["error"] for tes in test["tests"] if "error" in tes]
            if any(map(lambda tes: "t_total" in tes, test["tests"])):
                test["t_mean"] = numpy.mean([tes["t_total"] for tes in test["tests"] if "t_total" in tes])
            if any(map(lambda tes: "m_total" in tes, test["tests"])):
                test["m_mean"] = numpy.mean([tes["m_total"] for tes in test["tests"] if "m_total" in tes])
            if any(map(lambda tes: "diff_mean" in tes, test["tests"])):
                test["diff_mean"] = numpy.mean([tes["diff_mean"] for tes in test["tests"] if "diff_mean" in tes])
            if any(map(lambda tes: "diff_min" in tes, test["tests"])):
                test["diff_min"] = min([tes["diff_min"] for tes in test["tests"] if "diff_min" in tes])
            if any(map(lambda tes: "diff_max" in tes, test["tests"])):
                test["diff_max"] = max([tes["diff_max"] for tes in test["tests"] if "diff_max" in tes])
                print(f"\t\t\terror max: ±{test['diff_max']}")
            if any(map(lambda tes: "prec_out" in tes, test["tests"])):
                test["prec_out"] = bool(all([tes["prec_out"] for tes in test["tests"] if "prec_out" in tes]))
        except Exception as error:
            test["error"] = f"[{time.time()}]{type(error).__name__}: {error}"
    return test
      
        

#################################################
# Función para obtener los parametros de entrada optimos para
# la multiplicación de matrizes con el algoritmo HE Matrix Multiplication 
# segun los parametros de entrada y salida deseados
def optimizador_heMatMult(tam_lado = 2,        # Tamaño del lado de la matriz cuadrada a optimizar (esta optimización se podrá usar para matrizes de tamaño inferior manteniendo la precisión pero teniendo un sobrecoste de tiempo)
                          encriptacion = 128,   # 
                          distribucion_in = (lambda size: numpy.random.default_rng().random(size)), # Función para obtener una lista de 'size' datos en una distribución lo mas paracida a la de los datos para los que se optimiza
                          prec_in = None,       # Precision de los datos de entrada, expresada en 10^x. None: para mantener la precision maxima
                          prec_out = -10,        # Precision mínima de los datos de salida, expresada en 10^x.
                          repeticiones = 10,    # Repeticiones que se ejecutarán para comprobar la pecisión
                          pmd_estandar = True,  # Solo usar los poly_modbus degree estandares
                          save = True
                          ):
    if math.log2(tam_lado) != int(math.log2(tam_lado)):
        raise Exception(f"El tam_lado ({tam_lado}) para la optimización deve ser una potencia de 2")
    tam_list = tam_lado**2
    #list1, list2, prodMatr, matr1, matr2, _, _ = obtener_in_out(tam_lado, distribucion_in, prec_in, prec_out)
    testeos = {encriptacion: {}}
    minBitCount = MIN_BITS_COEFF_MODULUS*2
    for pmd in poly_modulus_degrees(encriptacion, tam_list, pmd_estandar):
        testeos[encriptacion][pmd] = {}
        print(f"poly_modulus_degree: {pmd}")
        
        #rranges = (slice(MIN_BITS_COEFF_MODULUS-1, MAX_BITS_COEFF_MODULUS+1, 1), slice(MIN_BITS_COEFF_MODULUS-1, MAX_BITS_COEFF_MODULUS+1, 1), slice(MIN_BITS_COEFF_MODULUS-1, MAX_BITS_COEFF_MODULUS+1, 1), slice(MIN_BITS_COEFF_MODULUS-1, MAX_BITS_COEFF_MODULUS+1, 1), slice(MIN_BITS_COEFF_MODULUS-1, MAX_BITS_COEFF_MODULUS+1, 1))
        #params =(pmd, encriptacion, repeticiones, tam_lado, distribucion_in, prec_in, prec_out)
        #t = optimize.brute(optCoeffModulus, rranges, args=params, full_output=True)
        #print(t)
        maxBitCount = CoeffModulus.MaxBitCount(pmd, encriptacion_types.get(encriptacion, sec_level_type.tc128))
        cml_valid = set()               # No dan error
        cml_pass = set()                # Pasan la restricción prec_out
        cml_best = (None, float("inf")) # La coeff_modulus con mejor resultado
        test = {}
        cml_l = math.ceil(maxBitCount/MAX_BITS_COEFF_MODULUS)
        cml = [math.floor(maxBitCount/cml_l)]*(cml_l-maxBitCount%cml_l) + [math.ceil(maxBitCount/cml_l)]*(maxBitCount%cml_l)
        pre_rec = set()
        continuar = True
        rapido = True
        diff_max_ant = None
        s_cml = ",".join(map(lambda x: str(x), cml))
        while True: #bucle de coeff_modulus
            try:
                if s_cml not in test or ("rapido" in test[s_cml] and not rapido):
                    coeff_modulus = CoeffModulus.Create(pmd, cml)
                    #cml = list(map(lambda x: x.bit_count(), coeff_modulus))
                    print(f"\tcoeff_modulus: {sum(cml)} => {s_cml}{' - Rapido' if rapido else ''}")
                    test[s_cml] = creaContextKeys(pmd, coeff_modulus, encriptacion, repeticiones, tam_lado, distribucion_in, prec_in, prec_out, pre_rec, rapido)
                    if len(test[s_cml]["precisiones_validas"]) > 0:
                        cml_valid.add(s_cml)
                        if len(test[s_cml]["precisiones_cumplen_requisitos"]) > 0:
                            cml_pass.add(s_cml)
                            pre_rec = test[s_cml]["precisiones_cumplen_requisitos"]
                            t_max = test[s_cml]["t_total"] + max(map(lambda prec: test[s_cml]["precisiones"][prec]["t_mean"], pre_rec))
                            #t_max = test[s_cml]["precisiones"][test[s_cml]["mejor_precision"]]["diff_max"]
                            if cml_best[1] > t_max:
                                cml_best = (s_cml, t_max)
                            print(f"\tt_max: {t_max}\n")
            except Exception as error:
                print(f"No se encontro resultado: {type(error).__name__} – {error}")
                test[s_cml] = {"error": f"{type(error).__name__} – {error}"}
                if len(cml_pass) == 0:
                    break
                elif len([cm for cm in map(lambda i: ",".join(map(lambda j: str(j), [math.floor(sum(cml)/i)]*(i-sum(cml)%i) + [math.ceil(sum(cml)/i)]*(sum(cml)%i))), range(math.ceil(sum(cml)/MAX_BITS_COEFF_MODULUS), math.floor(sum(cml)/MIN_BITS_COEFF_MODULUS)+1)) if cm in cml_pass]) == 0:
                    minBitCount = max(minBitCount,sum(cml))
                    continuar = True
                else:
                    maxBitCount = min(maxBitCount,sum(cml))
                    continuar = True
            finally:
                if len(cml_pass) > 0:
                    #if "error" not in test[s_cml] and (len(test[s_cml]["precisiones_cumplen_requisitos"]) == 0 or test[s_cml]["precisiones"][test[s_cml]["mejor_precision"]]["diff_max"] > 10**prec_out):
                    #    minBitCount = min(minBitCount,sum(cml))
                    #    print(f"minBitCount: {minBitCount}")
                    #    continuar = True
                    if continuar:
                        i = 0
                        while i < 5:
                            bitCount = round((maxBitCount-minBitCount)/2**(i+1))+minBitCount
                            diff_max_ant = None
                            for j in range(math.ceil(bitCount/MAX_BITS_COEFF_MODULUS), math.floor(bitCount/MIN_BITS_COEFF_MODULUS)+1):
                                cml = [math.floor(bitCount/j)]*(j-bitCount%j) + [math.ceil(bitCount/j)]*(bitCount%j)
                                s_cml = ",".join(map(lambda x: str(x), cml))
                                #cml = [math.floor(bitCount/len(cml))]*(len(cml)-bitCount%len(cml)) + [math.ceil(bitCount/len(cml))]*(bitCount%len(cml))
                                if s_cml not in test:
                                    break
                                elif "error" not in test[s_cml] and "mejor_precision" in test[s_cml] and test[s_cml]["mejor_precision"] != None:
                                    diff_max = test[s_cml]["precisiones"][test[s_cml]["mejor_precision"]]["diff_max"]
                                    if diff_max_ant == None:
                                        diff_max_ant = diff_max
                                    elif diff_max_ant >= diff_max:
                                        diff_max_ant = diff_max
                                    else:
                                        if len([cm for cm in map(lambda i: ",".join(map(lambda j: str(j), [math.floor(sum(cml)/i)]*(i-sum(cml)%i) + [math.ceil(sum(cml)/i)]*(sum(cml)%i))), range(math.ceil(sum(cml)/MAX_BITS_COEFF_MODULUS), math.floor(sum(cml)/MIN_BITS_COEFF_MODULUS)+1)) if cm in cml_pass]) == 0:
                                            if minBitCount < sum(cml):
                                                minBitCount = sum(cml)
                                                i = -1
                                        else:
                                            if maxBitCount > sum(cml):
                                                maxBitCount = sum(cml)
                                                i = -1
                                        break
                            if ",".join(map(lambda x: str(x), cml)) not in test:
                                continuar = True
                                rapido = False
                                break
                            i += 1
                        if ",".join(map(lambda x: str(x), cml)) in test and "rapido" not in test[",".join(map(lambda x: str(x), cml))]:
                            continuar = False
                    else:
                        break
                else:
                    diff_max = None
                    if "precisiones_validas" in test[s_cml] and len(test[s_cml]["precisiones_validas"]) > 0:
                        diff_max = min([test[s_cml]["precisiones"][prec]["diff_max"] for prec in test[s_cml]["precisiones_validas"]])
                    if diff_max_ant != None and (diff_max == None or diff_max_ant < diff_max):
                        minBitCount = sum(cml)
                        break
                    elif diff_max != None:
                        diff_max_ant = diff_max
                    continuar = True
                    pre_rec = set()
                    cml_l = len(cml)+1
                    minBitCount = MIN_BITS_COEFF_MODULUS*cml_l
                    cml = [math.floor(maxBitCount/cml_l)]*(cml_l-maxBitCount%cml_l) + [math.ceil(maxBitCount/cml_l)]*(maxBitCount%cml_l)
                    s_cml = ",".join(map(lambda x: str(x), cml))
                    #for i in reversed(range(1, len(cml))):
                    #    for j in reversed(range(i)):
                    #        if (cml[i]-1) > cml[j]:
                    #            cml[j] += 1
                    #            cml[i] -= 1
                    #            break
                    #        elif cml[i] == cml[j]:
                    #            break
                    #    if s_cml != ",".join(map(lambda x: str(x), cml)):
                    #        break
                    #if s_cml == ",".join(map(lambda x: str(x), cml)):
                    #    if MIN_BITS_COEFF_MODULUS * (len(cml)+1) > maxBitCount:
                    #        break
                    #    cml = [MIN_BITS_COEFF_MODULUS] * (len(cml)+1)
                    #    for i in reversed(range(len(cml))):
                    #        if sum(cml)+MAX_BITS_COEFF_MODULUS-MIN_BITS_COEFF_MODULUS <= maxBitCount:
                    #            cml[i] = MAX_BITS_COEFF_MODULUS
                    #        elif sum(cml) < maxBitCount:
                    #            cml[i] = maxBitCount - (sum(cml) - MIN_BITS_COEFF_MODULUS)
                    #        else:
                    #            break
                    #
                    #if s_cml == ",".join(map(lambda x: str(x), cml)):
                    #    break
        
        testeos[encriptacion][pmd]["coeff_modulus"] = test
        testeos[encriptacion][pmd]["coeff_modulus_validos"] = list(cml_valid)
        testeos[encriptacion][pmd]["coeff_modulus_cumplen_requisitos"] = list(cml_pass)
        testeos[encriptacion][pmd]["mejor_coeff_modulus"] = cml_best[0]
        print(f"\t\tcml_valid: {cml_valid}")
        print(f"\t\tcml_pass: {cml_pass}")
        print(f"\t\tcml_best: {cml_best}")
        if len(testeos[encriptacion][pmd]["coeff_modulus_cumplen_requisitos"]) > 0:
            break

    if save:
        with open(f"test_HEMatMult_{tam_lado}x{tam_lado}_{encriptacion}b_{repeticiones}rep{'_est' if pmd_estandar else ''}.json", "w") as outfile: 
            json.dump(testeos, outfile)
    
    


if __name__ == "__main__":
    try:
        optimizador_heMatMult()
    except Exception as error:
            print(f"No se encontro resultado: {type(error).__name__} – {error}")
