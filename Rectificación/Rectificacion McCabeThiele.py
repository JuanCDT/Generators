import numpy as np
from scipy import optimize
import random
#import pandas as pd
#import seaborn
import matplotlib.pyplot as mpl

#Enunciado del problema
'''Ha de proyectarse una columna de rectificación para separar 15000 kg/día de una
 mezcla de 40% de benceno y 60% de tolueno en un destilado conteniendo 97% de
 benceno y un residuo con 98% de tolueno. Todos los porcentajes son en masa.
 Se usará una razón de reflujo de 3,5 moles por cada mol de destilado.'''

A=15000/(24) #Caudal másico alimentación (kg/h)
wd = 0.97 # Fracción  másica del componente ligero en  el destilado
wr = 0.02 # Fracción  másica del componente ligero en la corriente residuo
wa = 0.4  # Fracción  másica del componente ligero en la alimentación




# Fracción de vapor en la alimentación
q = 1

#Pesos moleculares de los compuestos

PM_c1,PM_c2=78,92


#Función a partir de la cual determinar y para la mezcla conocida su volatilidad relativa

def equilib(x,y,alfa):
    f=y-alfa*x/(1+(alfa-1)*x)
    return f

def InterseccionRalimentoCurvaEquilibrio(x,q,xa,alfa):
    value=(alfa*x/(1+(alfa-1)*x))- RectaAlimentacion(x,q,xa)
    return value
    
    
def RectaAlimentacion(x,q,xa):
    y=-q/(1-q)*x+xa/(1-q)
    return y

def RazonReflujoMinima(q,xa,alfa,xd_a):
    # Cálculo de la razón de reflujo mínima para alcanzar los requisitos
    if (q<1 and q>0):
        x_AlEq=optimize.fsolve(InterseccionRalimentoCurvaEquilibrio,\
                               [0.5],args=(q,xa,alfa)) 
        x_AlEq=round(float(x_AlEq),3)
        y_AlEq=RectaAlimentacion(x_AlEq,q,xa)
        y_AlEq=round(float(y_AlEq),3)
    elif q==0:
        y_AlEq=xa
        y_AlEq=round(float(y_AlEq),3)
        x_AlEq=optimize.fsolve(equilib,[0.5],args=(y_AlEq,alfa))
        x_AlEq=round(float(x_AlEq),3)
    elif q==1:
        x_AlEq=xa
        x_AlEq=round(float(x_AlEq),3)
        y_AlEq=alfa*x_AlEq/(1+(alfa-1)*x_AlEq)
        y_AlEq=round(float(y_AlEq),3)
        
    AR=(xd_a-y_AlEq)/(xd_a-x_AlEq)
    RazonMinima=AR/(1-AR)
    
    return RazonMinima
    

def RectaOperativaEnriquecimiento(x,Ln,D,xd):
    y=Ln/(Ln+D)*x+(D*xd)/(Ln+D)
    return y


def RectaOperativaAgotamiento(x,Ln,D,xd,xa,A,q):
    y=((Ln+q*A)/(Ln+D-(1-q)*A))*x+(D*xd-A*xa)/(Ln+D-(1-q)*A)
    return y
    

def InterseccionRO(x,Ln,D,xd,xa,A,q):
    value=RectaOperativaEnriquecimiento(x,Ln,D,xd)-\
           RectaOperativaAgotamiento(x,Ln,D,xd,xa,A,q)
    return value

def Balances_Materia(incognitas,A,xa,xd,xr):
    D,R=incognitas
    
    values=[D+R-A]
    values.append(D*xd+R*xr-A*xa)
    
    return values
    

def generador_valores():  # generación de PM de los compuestos

    PM_c1=random.randint(50.,100.)
    PM_c2=random.randint(50.,100.)
    A=random.randint(10000,20000)
    q=round(random.random(),2)
    R=random.uniform(0,25)
    wa=round(random.random(),3)
    wd=round(random.uniform(0.75,0.95),3)
    wr=round(random.uniform(0.02,0.15),3)
    
    return PM_c1, PM_c2,A,q,R,wa, wd,wr

#PM_c1, PM_c2,A,q,R,wa= generador_PM()
    
def composicionesMolaresyCaudal(caudalMasico,wa,PM_c1,PM_c2):
    
    caudalMolar=round(((caudalMasico*wa)/PM_c1+(caudalMasico*(1-wa))/PM_c2),3)
    xa=round(((caudalMasico*wa)/PM_c1/caudalMolar),3)
    xb=round((1-xa),3)
    
    return caudalMolar,xa,xb
    
def composicionesMolares(wa,PM_c1,PM_c2):

    xa=round(((wa/PM_c1)/(wa/PM_c1+(1-wa)/PM_c2)),3)
    xb=round((1-xa),3)
    
    return xa,xb


caudalMolarA,xa,xb=composicionesMolaresyCaudal(A,wa,PM_c1,PM_c2)               


xd_a,xd_b=composicionesMolares(wd,PM_c1,PM_c2)
xr_a,xr_b=composicionesMolares(wr,PM_c1,PM_c2)



# Determinación de los puntos de la curva de equilibro
# Alfa es la volatilidad relativa
    
ye=[]
xe=[]

alfa = 2.45

RazonMinima=RazonReflujoMinima(q,xa,alfa,xd_a)
# Razón de reflujo de la columna
RazonReflujo = 1.5


for i in range(21):
    y = round(0.05 * i,3)
    ye.append(y)
    xe_i=optimize.fsolve(equilib,[0.5],args=(y,alfa)) 
    xe.append(round(float(xe_i[0]),3))
    
    
#
D,R=optimize.fsolve(Balances_Materia,[5,5],args=(caudalMolarA,xa,xd_a,xr_a))
#
D=round(float(D),3)
R=round(float(R),3)
#
fig1=mpl.figure()

mpl.plot(xe,ye, 'g-',label = 'Curva equilibrio')
mpl.plot(ye,ye, 'k-')
mpl.xlabel('Fracción molar en el líquido (xe)')
mpl.ylabel('Fracción molar en el vapor (ye)')

mpl.grid(b=True, which='both', color='0.65',linestyle='-')
mpl.xlim(0,1)
mpl.ylim(0,1)




Ln=round(RazonReflujo*D,3)
# Cálculo del punto de intersección entre las rectas de operación 

x_interseccion=optimize.fsolve(InterseccionRO,[0.5],args=(Ln,D,xd_a,xa,A,q))

x_interseccion=round(float(x_interseccion),3)
y_interseccion=round(RectaOperativaEnriquecimiento(x_interseccion,Ln,D,xd_a),3)

xAgotamiento=[xr_a,x_interseccion]
yAgotamiento=[xr_a,y_interseccion]
mpl.plot(xAgotamiento,yAgotamiento, 'r-',label = 'ROA')

xEnriquecimiento=[x_interseccion,xd_a]
yEnriquecimiento=[y_interseccion,xd_a]
mpl.plot(xEnriquecimiento,yEnriquecimiento, 'b-',label = 'ROE')
mpl.plot([xa,x_interseccion],[xa,y_interseccion],'k-',label = 'Recta q')




numeroEtapasEnriquecimiento=0
numeroEtapasAgotamiento=0
xDiagonal,yDiagonal=xd_a,xd_a #Valor inicial en la sección de enriquecimiento



while True:
    xequilibrio=optimize.fsolve(equilib,[0.5],args=(yDiagonal,alfa))  
    xequilibrio=float(xequilibrio)
    yequilibrio=yDiagonal
    mpl.plot([xDiagonal,xequilibrio],[yDiagonal,yequilibrio], 'c-')
    numeroEtapasEnriquecimiento+=1
    if xequilibrio>x_interseccion:
        xDiagonal=xequilibrio
        yDiagonal=RectaOperativaEnriquecimiento(xDiagonal,Ln,D,xd_a)
        mpl.plot([xequilibrio,xDiagonal],[yequilibrio,yDiagonal], 'c-')
    else:
        xDiagonal=xequilibrio
        yDiagonal=RectaOperativaAgotamiento(xDiagonal,Ln,D,xd_a,xa,caudalMolarA,q)
        mpl.plot([xequilibrio,xDiagonal],[yequilibrio,yDiagonal], 'c-')
        break

while True:
    xequilibrio=optimize.fsolve(equilib,[0.5],args=(yDiagonal,alfa))  
    xequilibrio=float(xequilibrio)
    yequilibrio=yDiagonal
    mpl.plot([xDiagonal,xequilibrio],[yDiagonal,yequilibrio], 'c-')
    numeroEtapasAgotamiento+=1
    if xequilibrio>=xr_a:
        xDiagonal=xequilibrio
        yDiagonal=RectaOperativaAgotamiento(xDiagonal,Ln,D,xd_a,xa,caudalMolarA,q)
        mpl.plot([xequilibrio,xDiagonal],[yequilibrio,yDiagonal], 'c-')
    else:
        break
    
    
mpl.legend(loc = 'best')
mpl.show()

numeroPisosTeoricos=numeroEtapasEnriquecimiento+numeroEtapasAgotamiento
PisoAlimentacion=numeroEtapasEnriquecimiento+1

