

import numpy as np
import random


def Calculo_incr_Tml(T1e, T1s, T2e,T2s):  # T1e,T1s,T2e,T2s: Temperaturas de entrada(e) y salida(s) de los fluidos caliente(1)
    # y frío(2). Sólo para flujo en contracorriente.

    incr_1 = T1s - T2e
    incr_2 = T1e - T2s
    if incr_1 == incr_2:
        incr_Tml = incr_1
    else:
        incr_Tml = (incr_1 - incr_2) / np.log(incr_1 / incr_2)

    return incr_Tml


def Cp_AguaLiquida(T):

    T_Kelvin=(T+273.15)/1000
    A = -203.6060
    B = 1523.290
    C = -3196.413
    D = 2474.455
    E = 3.855326
    Cp = A + B * T_Kelvin + C * T_Kelvin ** 2. + D * T_Kelvin ** 3. + E / (T_Kelvin ** 2.)  # En J/mol·K
    Cp_Masa = Cp / 18  # En kJ/kg·K
    return Cp_Masa


def generador_parametros():  # generación de parámetros del problema

    m1= random.randint(10000., 35000.) ; m2= random.randint(10000., 35000.) # Caudal másico kg/h


    while True:  # Generador T fluido caliente
        T1e = random.uniform(0, 100.)
        T1s = random.uniform(0, 100.)
        if T1e > T1s:
            break

    while True:  # Generador T fluido frío
        T2e = random.uniform(0, 100.)
        T2s = random.uniform(0, 100.)
        if T2e < T2s and T1e>T2s and T1s>T2e:
            break

    U = random.randint(1500., 2900.)  # Coeff. Global kJ/kg·m2·K



    ValoresInicialesTodos = np.array([1000., 1000., 50., 25., 25., 40.]) # Valores iniciales para caudales y Ts

    ValoresIniciales = np.array([100])  # Valor inicial para el Área de intercambio

    Area = float()
    incognitas = [Area]
    incognitasNombres = ['Area']

    variables = [m1, m2, T1e, T1s, T2e, T2s]
    variablesNombres = ['m1', 'm2', 'T1e', 'T1s', 'T2e', 'T2s']
    variablesConocidas = variables.copy()
    variablesConocidasNombres = variablesNombres.copy()

    variableDesconocida = random.randint(0, len(variables))

    print(variableDesconocida)

    incognitasNombres.append(variablesConocidasNombres[variableDesconocida])
    incognitas.append(variables[variableDesconocida])

    del variablesConocidas[variableDesconocida]
    del variablesConocidasNombres[variableDesconocida]


    variables.append(U)
    variablesConocidasNombres.append('U')

    ValoresIniciales = np.insert(ValoresIniciales, 1, ValoresInicialesTodos[variableDesconocida])

    valoresVariables = []
    valoresVariables.append(variables)
    valoresVariables.append(variableDesconocida)

    return variablesConocidasNombres,variables,incognitasNombres,incognitas,ValoresIniciales, variableDesconocida,valoresVariables


variablesConocidasNombres,variables,incognitasNombres,incognitas,ValoresIniciales, variableDesconocida,valoresVariables= generador_parametros()

print('Variable conocidas: ', variablesConocidasNombres, '\n', 'Incognitas: ', incognitasNombres,
      '\n Valores iniciales', ValoresIniciales, '\n Valores variables: ', variables,'\n Valores incognitas: ', incognitas,'\n Variable desconocida: ', variableDesconocida)

print(valoresVariables)


print (valoresVariables[1])

def Sistema_Ecuaciones(incognitas, valoresVariables):


    m1, m2, T1e, T1s, T2e, T2s, U, index=valoresVariables

    variables= [m1, m2, T1e, T1s, T2e, T2s, U]

    Area, variables[index],Q = incognitas

    T1Media=np.average(np.array([T1e, T1s]))
    T2Media = np.average(np.array([T2e, T2s]))

    Cp1=Cp_AguaLiquida(T1Media)
    Cp2 = Cp_AguaLiquida(T2Media)


    incr_Tfrio = T2s - T2e
    incr_Tcaliente = T1e - T1s
    incr_Tml = Calculo_incr_Tml(T1e, T1s, T2e, T2s)

    values=Q-m1 * Cp1 * incr_Tcaliente  # Ec. 1
    values.append = Q - m2 * Cp2 * incr_Tfrio  # Ec. 2
    values.append = Q - U * Area * incr_Tml  # Ec. 1

    return values

