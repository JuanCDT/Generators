import seaborn
import numpy as np
import matplotlib.pyplot as mpl
from scipy.optimize import curve_fit
import random
import math
import pandas as pd


def generador_intervaloGamma():
    gamma0 = round(random.uniform(1., 10.))
    gamma_final = round(random.uniform(100., 1000.))

    return gamma0,gamma_final


def generador_concentraciones(n_conc):
    concentraciones = []
    for i in range(n_conc):
        concentraciones.append(round(random.uniform(4, 20), 1))  # concentración en g/L

    Concentraciones_array = np.asarray(concentraciones)

    return np.sort(Concentraciones_array)



def ModeloPotencia(x, m,n):
    return m*x**n

def generador_concentracion():
    A_m= round(random.uniform(0.001, 0.01), 4)
    A_n= round(random.uniform(0.1, 0.95), 3)
    B_m= round(random.uniform(0.25,3.0), 2)
    B_n= round(random.uniform(0.05, 0.2), 2)
    parametros=[A_m,A_n,B_m,B_n]

    return parametros


def Calculo_R2( y, yteorica):
    residuals = y - yteorica  # Residuo: Experimentales menos predichos
    ss_res = np.sum(residuals ** 2)  # SSresiduo: Suma del cuadrado de los residuos
    ss_tot = np.sum((y - np.mean(y)) ** 2)  # Suma de cuadrados total
    r_squared = 1 - (ss_res / ss_tot)  # Cálculo de R2

    return r_squared


def parametros_Concentracion(C,A_m,A_n,B_m,B_n):  # Obtención de parámetros del modelo potencial: m y n, para una concentración
    global m, n
    m = round(A_m*C**B_m,1)
    n = round(A_n*C**B_n,3)

    return m, n


def generador_valores(n_ptos,C,parametros,gamma0,gamma_final):
    global m, n, tau, gamma
    A_m, A_n, B_m, B_n=parametros
    m, n = parametros_Concentracion(C,A_m,A_n,B_m,B_n)
    gamma = np.linspace(gamma0, gamma_final, n_ptos)
    tau = (m+random.uniform(-0.05,0.05)*m)*gamma**(n+random.uniform(-0.05,0.05)*n)

    return tau, gamma,m,n

#Generación de concentraciones
n_conc=int(input("Nº de concentraciones (máx. 5): "))

concentraciones=generador_concentraciones(n_conc)

gamma0,gamma_final=generador_intervaloGamma()

Valores_reograma=dict()


conc=dict()
conc['Tensión (N/m2)']=tau1.tolist()
conc['Velocidad de deformación (1/s)']=gamma1.tolist()


values = pd.DataFrame(data, columns=['Tensión (N/m2)', 'Velocidad de deformación (1/s)'])

#generación de valores para cada concentración

n_ptos = int(input("Nº de puntos para concentración: "))

tau, gamma, m, n = generador_valores(n_ptos)

eta = m * gamma ** (n - 1)

tau1 = np.copy(tau)
tau1 = np.round(tau1, 1)

gamma1 = np.round(gamma, 1)

data = dict()
data['Tensión (N/m2)'] = tau1.tolist()
data['Velocidad de deformación (1/s)'] = gamma1.tolist()

values = pd.DataFrame(data, columns=['Tensión (N/m2)', 'Velocidad de deformación (1/s)'])
# values.style.set_properties(**{'text-align': 'center'})
