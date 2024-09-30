import os
import psutil
import time
import numpy

def muestra_test(matriz1, matriz2, numeros, params={'pmd':8192, 'enc':128, 'pre':40}):
    if not params.get("silence", False):
        print("\nSTARTING Plain_computation...")
    
    ret={}
    ret["ini"] = time.time()

    ret["lib"]="Plain computation"
    ret["tam_matr1"] = f"{len(matriz1)}x{len(matriz1[0])}"
    ret["tam_matr2"] = f"{len(matriz2)}x{len(matriz2[0])}"
    ret["poly_modulus_degree"] = params["pmd"]
    ret["encriptacion"] = params["enc"]
    ret["precision"] = params["pre"]

    try:
        process = psutil.Process(os.getpid())
        m=process.memory_info().rss
        t=time.time()
        matr1_square=matriz1*matriz1
        ret["t_matr1_square_total"]=time.time()-t
        ret["m_matr1_square_total"]=process.memory_info().rss-m

        matr1_n_multiply={}

        for k in numeros.keys():
            m=process.memory_info().rss
            t=time.time()
            matr1_n_multiply[k]=numeros[k]*matriz1
            ret[f"t_matr1_{k}_multiply_total"]=time.time()-t
            ret[f"m_matr1_{k}_multiply_total"]=process.memory_info().rss-m

        n_filas=min(len(matriz1),len(matriz2))
        n_columnas=min(len(matriz1[0]),len(matriz2[0]))
        temp_matriz1=numpy.array(list(map(lambda x: x[:n_columnas], matriz1[:n_filas])))
        temp_matriz2=numpy.array(list(map(lambda x: x[:n_columnas], matriz2[:n_filas])))
        m=process.memory_info().rss
        t=time.time()
        matr1_matr2_multiply=temp_matriz1*temp_matriz2
        ret["t_matr1_matr2_multiply_total"]=time.time()-t
        ret["m_matr1_matr2_multiply_total"]=process.memory_info().rss-m

        if len(matriz1[0]) == len(matriz2):
            m=process.memory_info().rss
            t=time.time()
            matr1_matr2_matrix_multiply=numpy.matmul(matriz1, matriz2)
            ret["t_matr1_matr2_matrix_multiply_total"]=time.time()-t
            ret["m_matr1_matr2_matrix_multiply_total"]=process.memory_info().rss-m

        
        if not params.get("silence", False):
            print(f"\nOUTPUTS:")
            print(f"Elementos de la matriz 1 al cuadrado (x²):")
            print(f"    Total x² matriz 1:   {ret.get('t_matr1_square_total', 0) : <22.20f}s - {ret.get('m_matr1_square_total', 0)/1048576 : >9.4f}MiB")
            #print(f"Resultado:\n{matr1_square}"")

            print("\nMatriz 1 por escalares:")
            for k in numeros.keys():
                print(f"    {k}:\t{numeros[k]}")
                print(f"        Total:       {ret.get(f't_matr1_{k}_multiply_total') : <22.20f}s - {ret.get(f'm_matr1_{k}_multiply_total', 0)/1048576 : >9.4f}MiB")
                #print(f"Resultado:\n{matr1_n_multiply[k]}"")
            
            if(len(matriz1[0])==len(matriz2[0])):
                print("\nMatriz 1 por matriz 2, elemento a elemento:")
                print(f"    Total matr1*matr2:   {ret.get('t_matr1_matr2_multiply_total', 0) : <22.20f}s - {ret.get('m_matr1_matr2_multiply_total', 0)/1048576 : >9.4f}MiB")
                #print(f"Resultado:\n{matr1_matr2_multiply}")
            else:
                print("No se han podido multiplicar las matrices elemento a elemento ya que el tamaño de sus columnas no coincide.")
            
            if len(matriz1[0]) == len(matriz2):
                print("\nMatriz 1 por matriz 2, multiplicación de matrices:")
                print(f"    Total matr1*matr2:   {ret.get('t_matr1_matr2_matrix_multiply_total', 0) : <22.20f}s - {ret.get('m_matr1_matr2_matrix_multiply_total', 0)/1048576 : >9.4f}MiB")
            else:
                print(f"No se han podido realizar la multiplicación de matrices porque el numero de columnas de matr1 ({len(matriz1[0])}) y de filas de matr2 ({len(matriz2)}) no coincide.")
    
    except Exception as error:
        ret["errors"] = ret.get("errors", []) + [f"{time.time()}: {type(error).__name__}–{error}"]
        print(f"ERROR en Plain_computation: {type(error).__name__}–{error}")
    finally:
        ret["fin"] = time.time()
        return ret
