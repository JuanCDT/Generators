
import numpy as np
from scipy import optimize
from scipy.integrate import ode
import matplotlib.pyplot as mpl
import matplotlib.patches as patches
import random
import math


Cp_vapor=0.46 #Cp en kcal/kg·ºC
Cp_disolucion=1 #Cp en kcal/kg·ºC


def generador_parametros():  # generación de parámetros del problema de un evaporador de doble efecto

    CaudalAlimento = round(random.uniform(20000., 35000.))  # Caudal másico alimento
    W_SolutoAlimento = round(random.uniform(0.01, 0.05), 3)  # Fracción másica soluto en la alimentación
    W_SolutoProducto = round(random.uniform(0.35, 0.55), 3)  # CFracción másica soluto en la salida del segundo elemento
    Temperatura_Alimento = round(random.uniform(85., 95.), 1)  # Temperatura del alimento
    T_SaturacionVapor = round(random.uniform(110., 125.), 1)  # Temperatura de saturación del vapor de agua
    T_SaturacionCondensadorFinal = round(random.uniform(45., 75.), 1)  # Temperatura de saturación condensador final
    U1 = round(random.uniform(2500, 2900.))  # Coeff. Global Efecto 1
    U2 = round(random.uniform(1500., 2200.))  # Coeff. Global Efecto 2
    Teb1 = round(random.uniform(0., 10.), 1)  # Elevación del punto de ebullición en el primer efecto
    Teb2 = round(random.uniform(0., 15.), 1)  # Elevación del punto de ebullición en el segundo efecto

    return CaudalAlimento, W_SolutoAlimento, W_SolutoProducto, Temperatura_Alimento, T_SaturacionVapor, \
           T_SaturacionCondensadorFinal, U1, U2, Teb1, Teb2

# Calor latente
def Regnault(T):
    value=606.5-0.695*T #T en ºC y value en kcal/kg
    return value


#ecuaciones de balance de materia

def Balances(incognitas, A, Wa_s ,Wl2_s,Ta,T1,T5,Cp,U1,U2,teb1,teb2):
    L1, L2, V1, V2, Wl1_s, W, T2, T3, T4,  Area =incognitas
    CalorT1=Regnault(T1)
    CalorT2=Regnault(T2)
    CalorT4=Regnault(T4)

    # Tutil=T1-T5-(teb1+teb2)
    #Total
    values=[A*Wa_s-L2*Wl2_s] #Ec. 1
    values.append(A-(L2+V1+V2))  #Ec. 2

    # Evaporador 1
    values.append(A-(L1+V1))  #Ec. 3
    values.append(A*Wa_s-(L1*Wl1_s))  #Ec. 4

    # Entalpico Evaporador 1
    values.append(W*CalorT1+A*Cp*(Ta-T2)-CalorT2*V1)  #Ec. 5
    # Entalpico Evaporador 2
    values.append(V1*CalorT2+L1*Cp*(T3-T4)-CalorT4*V2)  #Ec. 6


    values.append(W*CalorT1 - U1 *Area*(T1-T2))  #Ec. 7
    values.append(V1*CalorT2 - U2 * Area * (T3 - T4))  #Ec. 8
    values.append(T3-(T2-teb1))  #Ec. 9
    values.append(T5-(T4-teb2) ) #Ec. 10


    return values


CaudalAlimento,W_SolutoAlimento,W_SolutoProducto,Temperatura_Alimento,T_SaturacionVapor,\
           T_SaturacionCondensadorFinal,U1,U2, Teb1,Teb2=generador_parametros()

print ("A(kg/h): ",CaudalAlimento, "; Soluto_A(%): ",'%.1f' % (W_SolutoAlimento*100),"; Temp. A(ºC): ",Temperatura_Alimento,\
       "; Soluto_L2(%): ", '%.1f' % (W_SolutoProducto*100),"; T1(ºC): ",T_SaturacionVapor,"; T5(ºC): ",\
       T_SaturacionCondensadorFinal,"\n Teb1(ºC): ",Teb1,"; Teb2(ºC): ",Teb2,"; U1(kcal/h·m2·ºC): ",\
       U1,"; U2(kcal/h·m2·ºC): ",U2)




Resultado = optimize.fsolve(Balances, [1000, 1000, 1000, 1000, 0.9, 1000, 50, 50, 50,  50],\
                                args=(CaudalAlimento, W_SolutoAlimento ,W_SolutoProducto,Temperatura_Alimento, \
                                T_SaturacionVapor,T_SaturacionCondensadorFinal,Cp_disolucion,U1,U2, \
                                      Teb1,Teb2), xtol=1e-06, maxfev=500)


L1, L2, V1, V2, Wl1_s, W, T2, T3, T4,  Area=Resultado





print ("\n L1:", round(L1), " L2:",round(L2)," V1:", round(V1)," V2:", round(V2), " Wl1_s:", round(Wl1_s,3), \
    "\n W:", round(W), " T2:", round(T2,1), " T3:", round(T3,1), " T4:", round(T4,1)," Area:",round(Area))
