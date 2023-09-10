
import numpy as np
from scipy import optimize
import pandas as pd
from scipy.integrate import ode
import matplotlib.pyplot as mpl
import matplotlib.patches as patches
import random
import math

Cp_disolucion = 1  # Cp en kcal/kg·ºC


def generador_parametros():  # generación de parámetros del problema de un evaporador de doble efecto

    CaudalAlimento = round(random.uniform(20000., 35000.))  # Caudal másico alimento
    W_SolutoAlimento = round(random.uniform(0.01, 0.05), 3)  # Fracción másica soluto en la alimentación
    W_SolutoProducto = round(random.uniform(0.35, 0.55), 3)  # CFracción másica soluto en la salida del segundo elemento
    Temperatura_Alimento = round(random.uniform(85., 95.), 1)  # Temperatura del alimento
    T_SaturacionVapor = round(random.uniform(110., 125.), 1)  # Temperatura de saturación del vapor de agua
    T_SaturacionCondensador = round(random.uniform(45., 75.), 1)  # Temperatura de saturación condensador
    U = round(random.uniform(2500, 2900.))  # Coeff. Global Efecto
    Teb = round(random.uniform(2, 10.), 1)  # Elevación del punto de ebullición

    return CaudalAlimento, W_SolutoAlimento, W_SolutoProducto, Temperatura_Alimento, T_SaturacionVapor, \
           T_SaturacionCondensador, U, Teb


def Regnault(T):
    value = 606.5 - 0.695 * T  # T en ºC y value en kcal/kg
    return value


def Balances(incognitas, A, Wa_s, Wl_s, Ta, T1, T3, Cp, U, teb):
    L, V, W, T2, Area = incognitas
    CalorT1 = Regnault(T1)
    CalorT2 = Regnault(T2)

    # Balance de materia
    values = [A * Wa_s - L * Wl_s]  # Ec. 1
    values.append(A - (L + V))  # Ec. 2

    # B. Entalpico
    values.append(W * CalorT1 + A * Cp * (Ta - T2) - CalorT2* V)  # Ec. 3
    # Ec. Trans. Calor
    values.append(W * CalorT1 - U * Area * (T1 - T2))  # Ec. 4

    # Ec. Inc. Temp. Ebullición
    values.append(T3 - (T2 - teb))  # Ec. 5

    return values


CaudalAlimento,W_SolutoAlimento,W_SolutoProducto,Temperatura_Alimento,T_SaturacionVapor,\
           T_SaturacionCondensador,U, Teb=generador_parametros()


print ("A(kg/h): ",CaudalAlimento, "; Soluto_A(%): ",'%.1f' % (W_SolutoAlimento*100),"; Temp. A(ºC): ",Temperatura_Alimento,\
       "; Soluto_L(%): ", '%.1f' % (W_SolutoProducto*100),"; T1(ºC): ",T_SaturacionVapor,"; T Saturación V (ºC): ", \
       T_SaturacionCondensador,"\n Teb(ºC): ",Teb,"; U(kcal/h·m2·ºC): ",U)




Resultado = optimize.fsolve(Balances, [1000,  1000, 1000, 50,  50],\
                                args=(CaudalAlimento, W_SolutoAlimento ,W_SolutoProducto,Temperatura_Alimento, \
                                T_SaturacionVapor,T_SaturacionCondensador,Cp_disolucion,U, Teb), xtol=1e-06, maxfev=500)


L, V, W, T2, Area = Resultado





print ("\n L: ", round(L), " V: ", round(V), "\n W: ", round(W), "T2: ", round(T2,1), "\n Area: ",round(Area), \
       "\n 'Calor latente(T1) (kcal/kg): ",int(np.round(Regnault(T_SaturacionVapor))),\
       " 'Calor latente(T2) (kcal/kg): ", int(np.round(Regnault(T2))))
