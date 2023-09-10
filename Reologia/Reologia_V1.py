import seaborn
import numpy as np
import matplotlib.pyplot as mpl
from scipy.optimize import curve_fit
import random
import math
import pandas as pd
import IPython.core.display as di
pd.set_option('display.notebook_repr_html', True)
import notebook
from IPython.display import clear_output, display, HTML, Image,Math, Latex
from IPython.external import mathjax
FigureSize=(15,5)




def ModeloPotencia(x, m,n):
    return m*x**n



def polyfit(x, y, degree):
    results = {}

    coeffs = np.polyfit(x, y, degree)
     # Polynomial Coefficients
    results['polynomial'] = coeffs.tolist()

    correlation = np.corrcoef(x, y)[0,1]


    results['correlation'] = correlation
    r_squared=correlation**2
    results['determination'] = correlation**2

    return r_squared,coeffs


def generador_parametros():  # generarción de parámetros del modelo potencial: m y n
    global m, n
    m = round(random.uniform(10., 1000.), 1)
    n = round(random.uniform(0.2, 1.4), 3)
    #n = round(random.uniform(1.1, 1.4), 3)
    #n = 1

    return m, n



def generador_valores(n_ptos):
    global m, n, tau, gamma
    m, n = generador_parametros()

    gamma0 = round(random.uniform(1., 3.))
    gamma_final = round(random.uniform(4., 20.))

    gamma = np.linspace(gamma0, gamma_final, n_ptos)
    tau = (m+random.uniform(-0.02,0.02)*m)*gamma**(n+random.uniform(-0.02,0.02)*n)
    #tau=tau+random.uniform(-0.05,0.05)*tau

    return tau, gamma,m,n



n_ptos = int(input("Nº de puntos: "))

tau, gamma, m, n = generador_valores(n_ptos)

eta = m * gamma ** (n - 1)

tau1 = np.copy(tau)
tau1 = np.round(tau1, 1)

gamma1 = np.round(gamma, 1)

data = dict()
data['Tensión (N·m2)'] = tau1.tolist()
data['Velocidad de deformación (1/s)'] = gamma1.tolist()

values = pd.DataFrame(data, columns=['Tensión (N·m2)', 'Velocidad de deformación (1/s)'])
# values.style.set_properties(**{'text-align': 'center'})

print (values)


fig1 = mpl.figure(figsize=FigureSize);

ax1 = fig1.add_subplot(121);
mpl.plot(gamma, tau, 'ro', label='Datos Experimentales')
mpl.xlabel('Velocidad de deformación, (1/s)')
mpl.ylabel('Tensión (N·m2)')
mpl.legend(loc='best')

ax2 = fig1.add_subplot(122);

mpl.plot(eta, tau, 'ro', label='Datos Experimentales')
mpl.xlabel('Velocidad de deformación, (1/s)')
mpl.ylabel('Viscosidad (Pa·s)')
mpl.legend(loc='best')
mpl.show()

fig2=mpl.figure()
mpl.plot(gamma,tau, 'bo',label = 'Datos Experimentales')
mpl.xlabel('Velocidad de deformación, (1/s)')
mpl.ylabel('Tensión (N·m2)')
mpl.legend(loc = 'best')
mpl.yscale('log')
mpl.xscale('log')
mpl.show()

popt, pcov = curve_fit(ModeloPotencia, gamma, tau)
m_experimental = popt[0]
n_experimental = popt[1]
Predichos = ModeloPotencia(gamma, m_experimental, n_experimental)  # Valores predichos por el modelo


perr = np.sqrt(np.diag(pcov))  # Error de cada uno de los parámetros del modelo: m y n
residuals = tau - Predichos  # Residuo: Experimentales menos predichos
ss_res = np.sum(residuals ** 2)  # SSresiduo: Suma del cuadrado de los residuos
ss_tot = np.sum((tau - np.mean(tau)) ** 2)  # Suma de cuadrados total
r_squared = 1 - (ss_res / ss_tot)  # Cálculo de R2

gamma_100 = np.linspace(gamma[0], gamma[-1], 100)
tau_Predichos=ModeloPotencia(gamma_100, m_experimental,n_experimental) # Valores predichos  para 100 valores
Viscosidad_predicha=m_experimental*gamma_100**(n_experimental-1)

# Representación Experimentales y predichos

fig3 = mpl.figure(figsize=FigureSize);

ax1 = fig3.add_subplot(121);
mpl.plot(gamma, tau, 'ro', label='Datos Experimentales')
mpl.plot(gamma_100 , tau_Predichos, 'k-', linewidth=1, label='$\mathrm{R^2}$ = ' + str(round(r_squared, 3)))

mpl.xlabel('Velocidad de deformación, (1/s)')
mpl.ylabel('Tensión (N·m2)')
mpl.legend(loc='best')

ax2 = fig3.add_subplot(122);
mpl.plot(gamma, eta, 'ro', label='Datos Experimentales')
mpl.plot(gamma_100 , Viscosidad_predicha, 'k-', linewidth=1, label='Viscosidad predicha')
mpl.xlabel('Velocidad de deformación, (1/s)')
mpl.ylabel('Viscosidad (Pa·s)')
mpl.legend(loc='best')

mpl.show()


print('m calculado: ' + str(round(m_experimental, 3)) +'+-'+ str(round(perr[0], 3)) \
     + '; n calculado: ' + str(round(n_experimental, 3)) +'+-'+ str(round(perr[1], 3)))

print (perr)
print('m: ', m, '\n n:', n)


