# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 10:02:20 2020

@author: Agnes
"""


import numpy as np
import matplotlib.pyplot as plt
import math
from mpl_toolkits import mplot3d

#Runge-Kutta
def rungekutta (f,x,function_parameters, h):
    """f: Funktion
        x: letzter Zeitschritt
        function_parameters: Parameter der rechten Handseite
        h: Schrittweite"""
    k1 = f(x, *function_parameters)
    k2 = f(x +(h/2)*k1, *function_parameters)
    k3 = f(x +(h/2)*k2, *function_parameters)
    k4 = f(x+ h*k3, *function_parameters)
    xnew = x + (h/6) *(k1 + 2*k2 +2*k3 +k4)
    return (xnew)

#Saltzmann
#Parameter
p= 0.1
q= 1.5
r= 1.5
# Startwerte
x =- math.sqrt(r-p) +0.0001     #Startwert x (proportional Eismasse)
y = math.sqrt(r-p)          #Startwert y (proportional CO2)
z= -math.sqrt(r-p)          #Satrtwert z (proportional Tiefsee-Temperatur)

#Durchführung des Runge-Kuttaverfahren
t=0 #Startpunkt in Zeit
tmax=100 #Endpunkt in Zeit
h = 0.001 #Schrittweite
xlist=[x]
ylist= [y]
zlist = [z]
tlist=[t]

def rhs1(x,y):      #Rechtehandseite der x-Änderung
   return(-x-y)

def rhs2(y,p,z,r):          #Rechtehandseite der y-Änderung
    return(p*z+r*y-(z**2)*y)

def rhs3(z,q,x):
   return(q*(x-z))     #Rechtehandseite der z-Änderung

while t<= tmax:

    t+=h
    tlist.append(t)
    #neues x,y,z ermitteln
    xnew=rungekutta(rhs1,x,[y],h)
    xlist.append(xnew)
    ynew=rungekutta(rhs2,y,[p,z,r],h)
    ylist.append(ynew)
    znew= rungekutta (rhs3,z,[q,x],h)
    zlist.append(znew)
   
    #Variablen updaten
    x= xnew*1
    y= ynew*1
    z= znew*1
    
#Diagramm
#Zeitreihe
plt.plot(tlist,xlist)
plt.xlabel("$t$") 
plt.ylabel("$x$")
plt.show()
plt.plot(tlist,ylist)
plt.xlabel("$t$") 
plt.ylabel("$y$")
plt.show()
plt.plot(tlist,zlist)
plt.xlabel("$t$") 
plt.ylabel("$z$")
plt.show()

#Phasenraumtrajektorie
plt.axes(projection='3d') 
plt.plot(xlist, ylist, zlist) 
plt.show()
   
    
    
