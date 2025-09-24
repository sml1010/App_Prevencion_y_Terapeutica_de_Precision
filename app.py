import pandas as pd

def lectura_csv(path):
    # Lectura del archivo csv con los datos de cada paciente.
    df = pd.read_csv("Genotype Matrix.csv", sep=';')
    # Los datos se pasan a un diccionario
    datos_dict = {df["Sample/Assay"][i]: {df.columns[j]: df[df.columns[j]][i] for j in range(1, len(df.columns))} for i in range(len(df))}
    # Tabla con las variantes y sus mutaciones asociadas
    tabla = {"CYP2D6": {"3": ["_"], "4": ["T"], "6": [ "_"], "7": ["G"],
                    "8": ["A"], "9": ["_"], "10*4": ["A"], "10": ["C"],
                    "12": ["T"], "14": ["T"], "15": ["_"], "17": ["A"],
                    "19": ["_"], "29": ["A"], "41": ["T"], "56B": ["A"],
                    "59": ["T"]},
         "UGT1A1": {"80": ["T"]},
         "DPYD": {"2A": [ "G","T" ], "13": ["C", "T"], "HapB3": ["T"], "D949V": ["A"]}}

    dict_pacientes = {}


    for i in datos_dict:
        dict_pacientes[i] = {}
        for j in datos_dict[i]:
            # Separacion de nombre de gen y variantes
            if "*" in j:
                separador = "*"
            elif "_" in j:
                separador = "_"
            snp = j.split(separador, 1)


            if snp[0] not in dict_pacientes[i]:
                dict_pacientes[i][snp[0]] = []

            # Manejo de Indeterminados
            if datos_dict[i][j] == "UND":
                dict_pacientes[i][snp[0]].append(("*1", "*1"))
                continue
            # Determinar genotipo 
            variante = []
            if datos_dict[i][j].split("/")[0] in tabla[snp[0]][snp[1]]:
                variante.append(f"{separador}{snp[1]}")
            else:
                variante.append(f"*1")

            if datos_dict[i][j].split("/")[1] in tabla[snp[0]][snp[1]]:
                variante.append(f"{separador}{snp[1]}")
            else:
                variante.append(f"*1")

            dict_pacientes[i][snp[0]].append(tuple(variante))

    # Eliminar duplicados
    for i in dict_pacientes:
        for j in dict_pacientes[i]:
            dict_pacientes[i][j] = list(set(dict_pacientes[i][j]))

    return dict_pacientes

def determinar_genotipo_definitivo(datos_pacientes):
    """
    Determina el genotipo definitivo para cada gen de cada paciente
    """
    resultados = {}
    
    for paciente, genes in datos_pacientes.items():
        resultados[paciente] = {}
        
        for gen, alelos in genes.items():
            # Recoger todos los alelos únicos, excluyendo *1 cuando hay otros
            alelos_unicos = list()
            
            for alelo_materno, alelo_paterno in alelos:
                # Si ambos son *1, es el caso base (sin mutaciones)
                if alelo_materno == '*1' and alelo_paterno == '*1':
                    continue
                
                # Añadir alelos no silvestres
                
                if alelo_materno != '*1':
                    alelos_unicos.append(alelo_materno)
                if alelo_paterno != '*1':
                    alelos_unicos.append(alelo_paterno)

                # Manejo de Variante *10*4
                if gen == "CYP2D6":
                    if "*10*4" in alelos_unicos and "*10" in alelos_unicos:
                        alelos_unicos.remove("*10*4")
                        if "*4" in alelos_unicos:
                            alelos_unicos.remove("*10")

                # Manejo de variante *80 como *28
                if gen == "UGT1A1" and alelos_unicos:
                    for i in range(len(alelos_unicos)):
                        alelos_unicos[i] = "*28"
               
            
            # Si no hay mutaciones, el genotipo es *1/*1
            if not alelos_unicos:
                resultados[paciente][gen] = ('*1', '*1')
            # Si hay una mutación, es heterocigoto *1/mutación
            elif len(alelos_unicos) == 1:
                mutacion = list(alelos_unicos)[0]
                resultados[paciente][gen] = ('*1', mutacion)
            # Si hay dos mutaciones, es heterocigoto mutación1/mutación2
            elif len(alelos_unicos) == 2:
                mutaciones = sorted(list(alelos_unicos))
                resultados[paciente][gen] = (mutaciones[0], mutaciones[1])
            
    
    return resultados

# Función para formatear el resultado como string
def formatear_genotipos(resultados):
    formateados = {}
    for paciente, genes in resultados.items():
        formateados[paciente] = {}
        for gen, alelos in genes.items():
            formateados[paciente][gen] = f"{alelos[0]}/{alelos[1]}"
    return formateados


# Procesar
dict_pacientes = lectura_csv("Genotype Matrix.csv")
resultados = determinar_genotipo_definitivo(dict_pacientes)
resultados_formateados = formatear_genotipos(resultados)

print(resultados_formateados)
