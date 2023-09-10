
import numpy as np
from scipy import optimize
import random
import pandas as pd

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

    m1 = random.randint(10., 200.);
    m2 = random.randint(10., 200.)  # Caudal másico kg/h

    while True:  # Generador T fluido caliente
        T1e = random.uniform(0, 100.)
        T1s = random.uniform(0, 100.)
        if T1e > T1s:
            break

    while True:  # Generador T fluido frío
        T2e = random.uniform(0, 100.)
        T2s = random.uniform(0, 100.)
        if T2e < T2s and T1e > T2s and T1s > T2e:
            break

    U = random.uniform(0.9, 1.5)  # Coeff. Global kW/m2·K

    ValoresInicialesTodos = np.array([10., 10., 50., 25., 25., 40.]) # Valores iniciales para caudales y Ts

    ValoresIniciales = np.array([10])  # Valor inicial para el Área de intercambio

    Area = float()
    incognitas = [Area]
    incognitasNombres = ['Area (m2)']

    variables = [m1, m2, T1e, T1s, T2e, T2s]
    variablesNombres = ['m1 (kg/s)', 'm2 (kg/s)', 'T1e (ºC)', 'T1s (ºC)', 'T2e (ºC)', 'T2s (ºC)']
    variablesConocidas = variables.copy()
    variablesConocidasNombres = variablesNombres.copy()

    variableDesconocida = random.randint(0, len(variables))

    incognitasNombres.append(variablesConocidasNombres[variableDesconocida])
    incognitas.append(variables[variableDesconocida])
    incognitasNombres.append('Q (kW)')

    del variablesConocidas[variableDesconocida]
    del variablesConocidasNombres[variableDesconocida]


    variables.append(U)
    variablesConocidasNombres.append('U (kW/m2·ºC)')

    ValoresIniciales = np.insert(ValoresIniciales, 1, ValoresInicialesTodos[variableDesconocida])
    ValoresIniciales = np.insert(ValoresIniciales, 2, 10)

    m1, m2, T1e, T1s, T2e, T2s,U=variables
    valoresVariables = [m1, m2, T1e, T1s, T2e, T2s,U,variableDesconocida]

    del variables[variableDesconocida]

    valoresVariables=tuple(valoresVariables)

    return variablesConocidasNombres,variables,incognitasNombres,incognitas,ValoresIniciales, variableDesconocida,valoresVariables




def Sistema_Ecuaciones(incognitas, m1, m2, T1e, T1s, T2e, T2s, U, index):

    variables= [m1, m2, T1e, T1s, T2e, T2s, U]
    Area, variables[index],Q = incognitas
    m1, m2, T1e, T1s, T2e, T2s,U=variables

    T1Media=np.average(np.array([T1e, T1s]))
    T2Media = np.average(np.array([T2e, T2s]))

    Cp1=Cp_AguaLiquida(T1Media)
    Cp2 = Cp_AguaLiquida(T2Media)


    incr_Tfrio = T2s - T2e
    incr_Tcaliente = T1e - T1s
    incr_Tml = Calculo_incr_Tml(T1e, T1s, T2e, T2s)

    values=[Q-m1 * Cp1 * incr_Tcaliente]  # Ec. 1
    values.append(Q - m2 * Cp2 * incr_Tfrio)  # Ec. 2
    values.append(Q - U * Area * incr_Tml)  # Ec. 1

    return values




variablesConocidasNombres,variables,incognitasNombres,incognitas,ValoresIniciales, variableDesconocida,valoresVariables= generador_parametros() 
m1, m2, T1e, T1s, T2e, T2s, U, index=valoresVariables

print(variables)

ValoresMostrados=[]

for i in range(len(variables)):
    if variableDesconocida>1:
        if i>1 and i<(len(variables)-1):
            ValoresMostrados.append(np.round(variables[i],1))
        elif i==(len(variables)-1):
                 ValoresMostrados.append(np.round(variables[i],3))
        else:
            ValoresMostrados.append(variables[i])
    else:
        if i>0 and i<(len(variables)-1):
            ValoresMostrados.append(np.round(variables[i],1))
        elif i==(len(variables)-1):
                 ValoresMostrados.append(np.round(variables[i],3))
        else:
            ValoresMostrados.append(variables[i])

data = dict(zip(variablesConocidasNombres, ValoresMostrados))

values = pd.DataFrame(data,index=['Datos'], columns=variablesConocidasNombres)


print(values,'\n')


titulos=['Incognita 1', 'Incognita 2','Incognita 3' ]


data1 = dict(zip(titulos,incognitasNombres))

values1 = pd.DataFrame(data1,index=[' '], columns=titulos)


print(values1,'\n')



Resultado = optimize.fsolve(Sistema_Ecuaciones, ValoresIniciales,args=valoresVariables, xtol=1e-06, maxfev=500)

Area,valorVariableDesconocida,Calor=Resultado

ValoresMostrar=[m1, m2, T1e, T1s, T2e, T2s]
ValoresMostrar[variableDesconocida]=valorVariableDesconocida
m1, m2, T1e, T1s, T2e, T2s=ValoresMostrar


incr_Tml=Calculo_incr_Tml(T1e, T1s, T2e,T2s)
T1Media=(T1e+T1s)/2
T2Media = (T2e+T2s)/2
Cp1=Cp_AguaLiquida(T1Media)
Cp2 = Cp_AguaLiquida(T2Media)
parametrosValores=[np.round(Cp1,3),np.round(Cp2,3),np.round(incr_Tml,2)]

    
valores_resultado=[int(np.round(Area))]

if variableDesconocida<=1:
    valores_resultado.append(int(np.round(valorVariableDesconocida)))
else:
    valores_resultado.append(float(np.round(valorVariableDesconocida,1)))
                                   
valores_resultado.append(int(np.round(Calor)))                              


parametrosNombres=['Cp1 (kJ/kg·ºC)','Cp2 (kJ/kg·ºC)','Incr. Temp. Media log. (ºC)']

data_parametros = dict(zip(parametrosNombres,parametrosValores))
values_parametros = pd.DataFrame(data_parametros ,index=['Parámetros'], columns=parametrosNombres)
print(values_parametros,'\n' )




data_resultado = dict(zip(incognitasNombres,valores_resultado))
values_resultado = pd.DataFrame(data_resultado ,index=['Resultados'], columns=incognitasNombres)
print(values_resultado,'\n')



