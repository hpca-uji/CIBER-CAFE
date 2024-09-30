import importlib
import math
import numpy
import sys
import csv
from functools import reduce
from columnar import columnar

ENCRIPTACION = 128
PRECISION = 40
MAT1_F = 64
MAT1_C = 64
MAT2_F = 64
MAT2_C = 64

def tests(args):
    params={}
    for a in ["-pmd", "--poly_modulus_degree"]:
        if a in args:
            params["pmd"] = int(args[args.index(a)+1])
            pmd = math.log2(params["pmd"])
            if int(pmd) != pmd:
                raise Exception(f"Se ha seleccionado un \"poly_modulus_degree\" de {params['pmd']}, este valor no es una potencia de dos.")
            break
    
    params["enc"] = ENCRIPTACION
    for a in ["-enc", "--encriptacion"]:
        if a in args:
            params["enc"] = int(args[args.index(a)+1])
            if params["enc"] not in [128, 192, 256]:
                raise Exception(f"Se ha seleccionado una encriptación de {params['enc']}bits y solo se permite encriptación de 128bits, 192bits y 256bits.")
            break
    
    params["pre"] = PRECISION
    for a in ["-pre", "--precision"]:
        if a in args:
            params["pre"] = int(args[args.index(a)+1])
            if params["pre"] <= 0:
                raise Exception(f"Se ha seleccionado una precisión de {params['pre']}bits y esta debe ser positiva.")
            break

    numeros = {
        #"13": 13.0,
        #"π+10": math.pi+10,
        #"53": 53.0,
        #"π+50": math.pi+50,
        #"103": 103.0,
        #"π+100": math.pi+100
    }

    params["fm1"] = MAT1_F
    for a in ["-fm1", "--fil_mat1"]:
        if a in args:
            params["fm1"] = int(args[args.index(a)+1])
            break
    params["cm1"] = MAT1_C
    for a in ["-cm1", "--col_mat1"]:
        if a in args:
            params["cm1"] = int(args[args.index(a)+1])
            break
    params["fm2"] = MAT2_F
    for a in ["-fm2", "--fil_mat2"]:
        if a in args:
            params["fm2"] = int(args[args.index(a)+1])
            break
    params["cm2"] = MAT2_C
    for a in ["-cm2", "--col_mat2"]:
        if a in args:
            params["cm2"] = int(args[args.index(a)+1])
            break
    if min([params["fm1"],params["cm1"],params["fm2"],params["cm2"]]) <= 0:
        raise Exception("La dimensión de las matrices no puede ser zero o inferior.")
            
    if not (any(map(lambda x: x in args, ["-lm1", "--load_mat1"])) and any(map(lambda x: x in args, ["-lm2", "--load_mat2"]))):
        rg = numpy.random.default_rng()
    
    if any(map(lambda x: x in args, ["-lm1", "--load_mat1"])):
        for a in ["-lm1", "--load_mat1"]:
            if a in args:
                with open(args[args.index(a)+1]) as f:
                    matriz1 = numpy.array(list(map(lambda x: list(map(lambda y: numpy.float64(y), x.split(";")[:params["cm1"]])), f.readlines()))[:params["fm1"]])
                if len(matriz1) != params["fm1"] or len(matriz1[0]) != params["cm1"]:
                    raise Exception(f"El fichero '{args[args.index(a)+1]}' no tiene suficientes numeros para cargar una matriz 1 de {params['fm1']}x{params['cm1']}, ha podido cargar una matriz de {len(matriz1)}x{len(matriz1[0])}")
                break
    else:
        matriz1 = rg.random((params["fm1"], params["cm1"]))
        #matriz1 = numpy.array(range(1, params["fm1"] * params["cm1"]+1)).reshape(params["fm1"], params["cm1"])

    if any(map(lambda x: x in args, ["-lm2", "--load_mat2"])):
        for a in ["-lm2", "--load_mat2"]:
            if a in args:
                with open(args[args.index(a)+1]) as f:
                    matriz2 = numpy.array(list(map(lambda x: list(map(lambda y: numpy.float64(y), x.split(";")[:params["cm2"]])), f.readlines()))[:params["fm2"]])
                if len(matriz2) != params["fm2"] or len(matriz2[0]) != params["cm2"]:
                    raise Exception(f"El fichero '{args[args.index(a)+1]}' no tiene suficientes numeros para cargar una matriz 2 de {params['fm2']}x{params['cm2']}, ha podido cargar una matriz de {len(matriz2)}x{len(matriz2[0])}")
                break
    else:
        matriz2 = rg.random((params["fm2"], params["cm2"]))
        #matriz2 = numpy.array(range(1, params["fm2"] * params["cm2"]+1)).reshape(params["fm2"], params["cm2"])
    
    for a in ["-sm1", "--sav_mat1"]:
        if a in args:
            with open(args[args.index(a)+1], "w") as f:
                f.write("\n".join(map(lambda x: ";".join(map(lambda y: str(y), x)), matriz1)))
            break

    for a in ["-sm2", "--sav_mat2"]:
        if a in args:
            with open(args[args.index(a)+1], "w") as f:
                f.write("\n".join(map(lambda x: ";".join(map(lambda y: str(y), x)), matriz2)))
            break

    for a in ["-cm", "--coeff_modulus"]:
        if a in args:
            params['cm'] = int(args[args.index(a)+1])
            if params["cm"] <= 0:
                raise Exception(f"El coeff modulus selecionado es {params['cm']} y no puede ser zero o inferior")
            break

    for a in ["-cma", "--coeff_modulus_array"]:
        if a in args:
            params['cma'] = list(map(lambda x: int(x), args[args.index(a)+1].split(",")))
            if any(map(lambda x: x <= 0, params['cma'])):
                raise Exception(f"En la lista de coeff modulus no puede haber ningun zero o inferior. coeff_modulus_array: {params['cma']}")
            if "cm" in params.keys() and params["cm"] != sum(params['cma']):
                raise Exception(f"Se ha definido un coeff modulus de {params['cm']} y és inconsistente con la lista de coeff modulus definida. Suma({params['cma']}) = {sum(params['cma'])} != {params['cm']}")
            if "cm" not in params:
                params["cm"] = sum(params['cma'])
            break
    
    for a in ["-cmp", "--coeff_modulus_primes"]:
        if a in args:
            params['cmp'] = list(map(lambda x: int(x), args[args.index(a)+1].split(",")))
            cma = list(map(lambda x: math.ceil(math.log2(x)), params['cmp']))
            if any(map(lambda x: x <= 0, params['cmp'])):
                raise Exception(f"En la lista de coeff_modulus_primes no puede haber ningun zero o inferior. coeff_modulus_primes: {params['cma']}")
            if "cm" in params.keys() and sum(cma) != params['cm']:
                raise Exception(f"Los coeff_modulus_primes definidos y el coeff_modulus no son consistentes. coeff_modulus: {params['cm']} coeff_modulus_primes: {params['cmp']}=>{cma}=>{sum(cma)} != {params['cm']}")
            if "cma" in params.keys() and (len(params['cmp'] ) != len(params['cma']) or any([x != y for x, y in zip(cma, params['cma'])])):
                raise Exception(f"Los coeff_modulus_primes definidos y el coeff_modulus_array no son consistentes. coeff_modulus_array: {params['cma']} coeff_modulus_primes: {params['cmp']}")
            if "cm" not in params:
                params["cm"] = sum(cma)
            if "cma" not in params:
                params["cma"] = cma
            break
    
    if "-s" not in args:
        print(f"CONFIGURATION:")
        print(f"Tamaño de la matriz 1: {params['fm1']}x{params['cm1']}")
        print(f"Tamaño de la matriz 2: {params['fm2']}x{params['cm2']}")
        print(f"encriptacion: {params['enc'] : >6}b")
        print(f"precision: {params['pre'] : >9}b\t=> escala: 2^{params['pre']} = {pow(2, params['pre'])}")
        if any(map(lambda x: x in params, ["cm", "cma", "cmp"])):
            print(f"coeff_modulus: {params.get('cm', None) : >6}b => {params.get('cma', None)} => {params.get('cmp', None)}")
        else:
            print("coeff_modulus:\tCalculado por las libreria")

        print("\nTESTING VALUES:")
        for k in numeros.keys():
            print(f"{k}:\t{numeros[k]}")
        print(f"Matriz 1 ({len(matriz1)}x{len(matriz1[0])}):\n{matriz1}")
        print(f"Matriz 2 ({len(matriz2)}x{len(matriz2[0])}):\n{matriz2}")
    else:
        params["silence"] = True
    
    if any(map(lambda x: x in ["-a", "--all", "-P", "--plain", "-S", "--seal", "-E", "--eva"],args)):
        if "pmd" not in params.keys():
            params["pmd"] = max(16384 if params["enc"] == 256 else 8192, 2**math.ceil(math.log2(max([params["fm1"]*params["cm1"], params["fm2"]*params["cm2"],params["fm1"]*params["cm2"]])*2)))

        res = []
        
        if any(map(lambda x: x in ["-P", "--plain", "-a", "--all"], args)):
            try:
                res += [importlib.import_module("Plain_computation_Test").muestra_test(matriz1, matriz2, numeros, params)]
            except Exception as error:
                res += [{'lib':'Plain_computation','tam_matr1':f'{params["fm1"]}x{params["cm1"]}','tam_matr2':f'{params["fm2"]}x{params["cm2"]}','parametros':params,'errores':[f'{type(error).__name__}–{error}']}]
                print(f"El test sin encriptar falló: {type(error).__name__}–{error}")
        if any(map(lambda x: x in ["-E", "--eva", "-a", "--all"], args)):
            try:
                res += [importlib.import_module("EVA_Test").eva_test(matriz1, matriz2, numeros, params)]
                if any(map(lambda x: x in ["-S", "--seal", "-a", "--all"],args)):
                    #log2_pmd_max = int(math.log2(res[-1]["poly_modulus_degree_ctx"]))
                    try:
                        params_t = params.copy()
                        params_t.update({'pmd':res[-1]["poly_modulus_degree_ctx"], 'cma':res[-1]["coeff_modulus_arr"]})
                        res += [importlib.import_module("SEAL_Test").seal_test(matriz1, matriz2, numeros, params_t)]
                    except Exception as error:
                        res += [{'lib':'SEAL','tam_matr1':f'{params["fm1"]}x{params["cm1"]}','tam_matr2':f'{params["fm2"]}x{params["cm2"]}','parametros':params_t,'errores':[f'{type(error).__name__}–{error}']}]
                        print(f"El test de SEAL falló: {type(error).__name__}–{error}")
                    
                    for pmd in map(lambda x: 2**x, range(int(math.log2(params["pmd"])), 16)):
                        try:
                            params_t = params.copy()
                            params_t.update({'pmd':pmd})
                            res += [importlib.import_module("SEAL_Test").seal_test(matriz1, matriz2, numeros, params_t)]
                        except Exception as error:
                            res += [{'lib':'SEAL','tam_matr1':f'{params["fm1"]}x{params["cm1"]}','tam_matr2':f'{params["fm2"]}x{params["cm2"]}','parametros':params_t,'errores':[f'{type(error).__name__}–{error}']}]
                            print(f"El test de SEAL falló: {type(error).__name__}–{error}")
            except Exception as error:
                res += [{'lib':'EVA','tam_matr1':f'{params["fm1"]}x{params["cm1"]}','tam_matr2':f'{params["fm2"]}x{params["cm2"]}','parametros':params,'errores':[f'{type(error).__name__}–{error}']}]
                print(f"El test de EVA falló: {type(error).__name__}–{error}")

        if any(map(lambda x: x in ["-S", "--seal", "-a", "--all"], args)) and not any(map(lambda x: x in ["-E", "--eva"], args)):
            try:
                res += [importlib.import_module("SEAL_Test").seal_test(matriz1, matriz2, numeros, params)]
            except Exception as error:
                res += [{'lib':'SEAL','tam_matr1':f'{params["fm1"]}x{params["cm1"]}','tam_matr2':f'{params["fm2"]}x{params["cm2"]}','parametros':params,'errores':[f'{type(error).__name__}–{error}']}]
                print(f"El test de SEAL falló: {type(error).__name__}–{error}")

        if "-csv" in args:
            n_cab = list(map(lambda x: list(x.keys()), res))
            n_cab.sort(key=len, reverse=True)
            try:
                f = open(args[args.index("-csv")+1], "x")
                writer = csv.DictWriter(f, fieldnames=reduce(lambda x, y: x + [e for e in y if e not in x], n_cab), delimiter=";")
                writer.writeheader()
            except Exception as error:
                f = open(args[args.index("-csv")+1], "r+")
                cab = f.readline().replace("\n", "").split(";")
                n_cab = reduce(lambda x, y: x + [e for e in y if e not in x], ([cab]+n_cab))
                if(len(cab) == len(n_cab)):
                    f.read()
                    writer = csv.DictWriter(f, fieldnames=cab, delimiter=";")
                else:
                    f.seek(0)
                    res = [row for row in csv.DictReader(f, delimiter=";")] + res
                    f.seek(0)
                    f.truncate()
                    n_cab = list(map(lambda x: list(x.keys()), res))
                    n_cab.sort(key=len, reverse=True)
                    writer = csv.DictWriter(f, fieldnames=reduce(lambda x, y: x + [e for e in y if e not in x], n_cab), delimiter=";")
                    writer.writeheader()
            finally:
                writer.writerows(res)
                f.close()
        
if __name__ == "__main__":
    if len(sys.argv[1::])==0 or \
            any(map(lambda x: x in sys.argv[1::], ["-h", "--help"])) or \
            all(map(lambda x: x not in sys.argv[1::], ["-a", "--all", "-P", "--plain", "-S", "--seal", "-E", "--eva", "-HL", "--helib","-HN", "--heaan","-O", "--openfhe", "-sm1", "--sav_mat1", "-sm2", "--sav_mat2"])):
        help = [
            ["-pmd {int}, --poly_modulus_degree {int}", "indica el grado del modulo de polinomio usado, este numero debe ser un entero positivo potencia de 2, por defecto se calcula en base a el tamaño maximo necesario para almacenar las matrizes."],
            ["-enc {int}, --encriptacion {int}", f"indica el nivel de seguridad de la encriptación en bits, solo son validos los valores \"128\", \"192\" y \"256\", por defecto {ENCRIPTACION}."],
            ["-pre {int}, --precision {int}", f"numero de bits de precisión para los valores a calcular, por defecto {PRECISION}."],
            ["-cm {int}, --coeff_modulus {int}", ", por defecto calcula los mas adecuados cada libreria."],
            ["-cma {vector<int>}, --coeff_modulus_array {vector<int>}", "indica el vector de numeros primos que se usarán para el coeff_modulus, los numeros deben ser primos positivos separados por ',', por defecto calcula los mas adecuados cada libreria."],
            ["-cmp {vector<int>}, --coeff_modulus_primes {vector<int>}", "indica el vector de numeros primos que se usarán para el coeff_modulus, los numeros deben ser primos positivos separados por ',', por defecto calcula los mas adecuados cada libreria."],
            ["-fm1 {int}, --fil_mat1 {int}", f"indica el numero de filas de la matriz 1 con la que ejecutar la pruebas, por defecto {MAT1_F}."],
            ["-cm1 {int}, --col_mat1 {int}", f"indica el numero de columnas de la matriz 1 con la que ejecutar la pruebas, por defecto {MAT1_C}."],
            ["-sm1 {str}, --sav_mat1 {str}", "indica que guarde la matriz 1 generada en el fichero indicado."],
            ["-lm1 {str}, --load_mat1 {str}", "indica que carge la matriz 1 desde el fichero indicado."],
            ["-fm2 {int}, --fil_mat2 {int}", f"indica el numero de filas de la matriz 2 con la que ejecutar la pruebas, por defecto {MAT2_F}."],
            ["-cm2 {int}, --col_mat2 {int}", f"indica el numero de columnas de la matriz 2 con la que ejecutar la pruebas, por defecto {MAT2_C}."],
            ["-sm2 {str}, --sav_mat2 {str}", "indica que guarde la matriz 2 generada en el fichero indicado."],
            ["-lm2 {str}, --load_mat2 {str}", "indica que carge la matriz 2 desde el fichero indicado."],
            ["-a, --all", "ejecuta todas las pruebas."],
            ["-P, --plain", "ejecuta la prueba si ninguna encriptación"],
            ["-S, --seal", "ejecuta la prueba de SEAL."],
            ["-E, --eva", "ejecuta la prueba de EVA - SEAL."],
            ["-s", "modo silencioso, no imprime nada por pantalla."],
            ["-csv {str}", "almacena los resultados en una nueva linea del fichero indicado con formato CSV"],
            ["-h, --help", "muestra esta ayuda."]
        ]
        print(f"Testeo de diferentes librerias de cifrado homomorfico:{columnar(help, no_borders=True)}")
    else:
        try:
            tests(sys.argv[1::])
        except Exception as error:
            print(f"La inicialización de los tests falló: {type(error).__name__}–{error}")
