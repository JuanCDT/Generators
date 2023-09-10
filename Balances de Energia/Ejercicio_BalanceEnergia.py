
import numpy as np
from scipy import optimize
import random
import pandas as pd
import matplotlib.pyplot as mpl
from scipy.interpolate import interp1d



def P_Saturacion(T):
    #Temperatura en ºC y P en bar
    T_Kelvin=T + 273.15
    if T_Kelvin>=255.9 and T_Kelvin<=373.0:
        A=4.6543
        B=1435.264
        C=-64.848
        P = 10 ** (A - (B / (T_Kelvin + C)))

    elif  T_Kelvin>=379.0 and  T_Kelvin<=573.0:
        A= 3.55959
        B=643.748
        C=-198.043
        P = 10 ** (A - (B / (T_Kelvin + C)))

    elif T_Kelvin>373.0 and  T_Kelvin<379.0:
        x1=np.linspace(363, 373, num=21)
        x2=np.linspace(379, 389, num=21)
        x = np.append(x1, x2)
        y=P_saturacion(x)
        f = interp1d(x, y, kind='quadratic')
        P=f(T_Kelvin)
    else:
        raise ValueError('Temperatura fuera del intervalo válido')


    return P


def T_Saturacion(P):
    #P en kPa
    T_sat=(42.6776-(3892.7/(np.log(P/1000)-9.48654)))-273.15
    return T_sat


def Entalpia_AguaLiquida(T):
    #Correlación válida de 0.01 a 200ºC
    H_agua_Saturada=-0.033635409+4.207557011*T-6.200339E-4*T**2+4.459374E-6*T**3 #Entalpia kJ/kg
    return H_agua_Saturada

def Entalpia_Vapor(T):
    # Correlación válida de 0.01 a 200ºC
    H_vapor_Saturado=2501.6898445+1.806916015*T+5.087717E-4*T**2-1.1221E-5*T**3
    return H_vapor_Saturado



def Calor_Latente_Vaporizacion(T):
    # Correlación válida de 0.01 a 200ºC
    Lamda=2501.897149-2.407064037*T+1.192217E-3*T**2-1.5863E-5*T**3 #kJ/kg
    return Lamda

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


def Cp_AguaVapor(T):
    T_Kelvin=(T+273.15)/1000
    A = -203.6060
    B = 1523.290
    C = -3196.413
    D = 2474.455
    E = 3.855326
    Cp = A + B * T_Kelvin + C * T_Kelvin ** 2. + D * T_Kelvin ** 3. + E / (T_Kelvin ** 2.)  # En J/mol·K
    Cp_Masa = Cp / 18  # En kJ/kg·K
    return Cp_Masa


def Cp_vapor(T,P):
    # T en ºC y P en Pa
    # Correlación válida de 0  a 200 ºC y 200 kPa

    if T<50:
        AE = 1877.2
        BE = -0.49545
        CE = 8.1818E-3
        AF = 22.537
        BF= 0.49321
        CF=0.048927
    else:
        AE = 1856.1
        BE = 0.28056
        CE = 6.9444E-4
        AF = 22.537
        BF= 0.49321
        CF=0.048927
    Ptr=611.657
    Cp_Masa =(AE+BE*T +CE*T**2+(1/(AF+BF*T +CF*T**2))*(P-Ptr))/1000
     # En kJ/kg·K
    return Cp_Masa

# mV es el caudal de vapor (kg/h)
# mL es el caudal del medio (propiedades agua líquida) (kg/h)
# Tv es la temperatura del vapor saturado (ºC)
# Pv es la presión del vapor saturado (kPa)
# Qp es el caudal de pérdidas de calor (kW)
# T1 es la temperatura de entrada del medio (ºC)
# T2 es la temperatura de salida del medio (ºC)
# lamdaV es el calor latente de vaporización del agua (kJ/kg)
# CpV es la capacidad calorífica del vapor  (kJ/kg ºC)
# CpL es la capacidad calorífica del agua líquida (kJ/kg ºC)



def generador_parametros():  # generación de parámetros del problema

    mV=random.randint(100, 300.)
    Tv =round(random.uniform(100, 150),1)
    Pv=round(P_Saturacion(Tv),3)
    lamdaV=round(Calor_Latente_Vaporizacion(Tv))
    h2=round(Entalpia_Vapor(Tv))
    h1= round(Entalpia_AguaLiquida(Tv))
    Qcedido=round((mV*lamdaV))
    Qp=round((round(random.uniform(0.05, 0.2),3)*Qcedido))
    Qganado=round(Qcedido-Qp)
    T4 = round(random.uniform(40, 100), 1)
    incrT=round(random.uniform(10, 40), 1)
    T3=round(T4-incrT,1)
    Tmedia=incrT/2
    CpL=round(Cp_AguaLiquida(Tmedia),2)
    h3=round(Entalpia_AguaLiquida(T3))
    h4=round(Entalpia_AguaLiquida(T4))
    mL=round(Qganado/(CpL*incrT))
    return mV,Tv,Pv,lamdaV,h2,h1,Qcedido,Qp,Qganado,T4,incrT,T3,CpL,h3,h4,mL


#
#
# def generador_conversion():
#     mv = random.uniform(100, 1.0)
#     return  conversion
#
#
#
# def Sistema_Ecuaciones(incognitas,variables_conocidas,coeficientes,Pesos_moleculares,conversion_A):
#     #Variables conocidas y desconocidas dependen del caso pero se resuelve siempre para las mismas conocidas
#     # y generadas.
#     #Variables del sistema: C1,C2,C3,C4,C5
#     # W1_a, W1_b, W1_c, W1_d
#     # W2_a, W2_b, W2_c, W2_d
#     # W3_a, W3_b, W3_c, W3_d
#     # W4_a, W4_b, W4_c, W4_d
#     # W5_a, W5_b, W5_c, W5_d
#
#     # W3_a=W4_a=W5_a Por lo que se toma en todas las ecuaciones w3_a
#     # W3_b=W4_b=W5_b Por lo que se toma en todas las ecuaciones w3_b
#     # W3_c=W4_c=W5_c Por lo que se toma en todas las ecuaciones w3_c
#     # W3_d=W4_d=W5_d Por lo que se toma en todas las ecuaciones w3_d
#
#     # 14 variables de entrada, 5 de ellas definidas, y 9 ecuaciones y 9 incognitas
#     # Siempre se general, además de la estequiometría y la conversión, el caudal y 3 composiciones de C2 y C5
#     # Una vez resuelto, según el caso, se muestran al estudiante unos valores de partida u otros
#
#     X=conversion_A
#
#     a=coeficientes['a'];b=coeficientes['b'];c=coeficientes['c']
#     PM_a=Pesos_moleculares['a'];PM_b=Pesos_moleculares['b'];PM_c=Pesos_moleculares['c']
#
#     C1,C3,C4,W1_a, W1_b, W1_c,W3_a, W3_b, W3_c =incognitas
#     C2, C5, W2_a, W2_b, W2_c=variables_conocidas
#
#     #Ecuaciones en el mezclador
#     values=[C1*W1_a+C5*W3_a-C2*W2_a] #Ecuación 1
#     values.append(C1*W1_b+C5*W3_b-C2*W2_b)#Ecuación 2
#     values.append(C1*W1_c+C5*W3_c-C2*W2_c) #Ecuación 3
#     values.append(C1+C5-C2)#Ecuación 4 - balance Global
#
#     #Ecuaciones en el divisor
#     values.append(C3-C4-C5) #Ecuación 5 Balance global
#
#
#     #Ecuaciones en el reactor
#     values.append(C3*W3_a/PM_a+C2*W2_a*X/PM_a-C2*W2_a/PM_a) #Ecuación 6
#     values.append(C3*W3_b/PM_b+((b/a)*C2*W2_a*X)/PM_a-C2*W2_b/PM_b) #Ecuación 7
#     values.append(C3*W3_c/PM_c-((c/a)*C2*W2_a*X)/PM_a-C2*W2_c/PM_c) #Ecuación 8
#     values.append(C3-C2) #Ecuación 9 Balance global
#
#     return values
#
# def generador_problema():
#
#     comprobador=False
#
#     while comprobador==False:
#         C2, C5, W2_a, W2_b, W2_c,W2_d =generador_parametros()
#
#         coeficientes,Pesos_moleculares=generador_estequiometria()
#
#         Limitante=Reactivo_limitante(C2, W2_a,W2_b,coeficientes,Pesos_moleculares)
#
#         #print ("El reactivo limitante es: ",Limitante)
#
#         conversion=generador_conversion(Limitante, W2_a,W2_b,coeficientes,Pesos_moleculares)
#
#         variables_conocidas = [C2, C5, W2_a, W2_b, W2_c]
#         #
#         Resultado = optimize.fsolve(Sistema_Ecuaciones, [100, 100, 100, 0.3, 0.3, 0.1,0.3, 0.3, 0.1],\
#                                   args=(variables_conocidas,coeficientes,Pesos_moleculares,conversion),\
#                                  xtol=1e-06, maxfev=500)
#         #
#         #
#         C1,C3,C4,W1_a, W1_b, W1_c,W3_a, W3_b, W3_c  = Resultado
#
#         W1_d = 1-(W1_a + W1_b + W1_c)
#         W3_d = 1-(W3_a + W3_b + W3_c)
#         R=C5/C3
#         corrientes=[C1,C2, C3,C4,C5]
#         fracciones=[ W1_a, W1_b, W1_c,W1_d , W2_a, W2_b, W2_c,W2_d ,W3_a, W3_b, W3_c,  W3_d,R ]
#         corrientes=[int(elem) for elem in corrientes]
#         fracciones=[round(float(elem), 3) for elem in fracciones]
#         comprobador=all(i >= 0 for i in fracciones)
#
#     return fracciones,corrientes, coeficientes,Pesos_moleculares,round(conversion,3)
#
#

mV,Tv,Pv,lamdaV,h2,h1,Qcedido,Qp,Qganado,T4,incrT,T3,CpL,h3,h4,mL=generador_parametros()

print ("\n mv: ", mV, "\n Tv: ", Tv, "\n Pv: ", Pv, "\n lamdaV: ", lamdaV, "\n h2: ",h2, \
      "\n h1: ",h1,"\n Qcedido: ", Qcedido, "\n h3: ",h3, "\n h4: ",h4,  \
       "\n Qp: ", Qp, "\n Qganado: ", Qganado, "\n T4: ", T4, "\n T3: ", T3, "\n CpL: ", CpL,"\n mL: ", mL)


variablesConocidasNombres = ['mV (kg/h)', 'T1 (ºC)', ' h agua líquida, h1 (kJ/kg) ', ' \
h vapor saturado, h2 (kJ/kg) ', 'T2 (ºC)', 'T3 (ºC)', 'Pérdidas de calor (kJ/h)']
ValoresMostrados = [mV, Tv, h1, h2, T3, T4, Qp]

data = dict(zip(variablesConocidasNombres, ValoresMostrados))


values = pd.DataFrame(data,index=['Coef. Estequiométricos'], columns=variablesConocidasNombres)

print (values)
