# -*- coding: utf-8 -*-
"""
Created on Sun May 24 14:02:37 2020

@author: Agnes
"""
import numpy as np
import matplotlib.pyplot as plt

#Runge-Kutta-Verfahren
def rungekutta (f,x,function_parameters,h):
    """ f: Funktion
        x: letzter Zeitschritt
        function_parameters: Liste der Parameter
        h: Schrittweite"""
    k1 = f(x,*function_parameters)
    k2= f(x+(h/2)*k1, *function_parameters)
    k3 = f(x+(h/2)*k2, *function_parameters)
    k4= f(x+ h*k3, *function_parameters)
    xnew = x + (h/6)*(k1+ 2*k2 +2*k3+ k4)
    return xnew

#Räuberbeutemodell
#Parameter
e1 = 2.0
e2 = 0.8
y1= 0.02
y2= 0.0002
h = 0.025 # Schrittweite
p1= 100 #Startwert der Beutepopulation
p2= 50 #Startwert der Räuberpopulation
tmax = 100


#Ausführung des RK-Verfahren
t=0
p1list = [p1]
p2list = [p2]
tlist = [t]

while t<=tmax:
    def rhs1(p1,e1,y1,p2):  #Definiton Rechtehandseite Beutepopulation
        return p1*(e1-y1*p2)
    def rhs2(p2,e2,y2,p1new):  #Defintion Rechtehanseite Räuberpopulation
        return (-p2)*(e2-y2*p1new)
    
    p1new = rungekutta(rhs1,p1,[e1,y1,p2],h)
    p1list.append(p1new)
    
    p2new = rungekutta(rhs2,p2,[e2,y2,p1new],h)
    p2list.append(p2new)
    
    
    p1 = p1new +0
    p2 = p2new +0

    t+=h #Zeitschritt weitergehen
    tlist.append(t)
    
#Diagramm
# p(t)
plt.plot(tlist,p1list)
plt.plot(tlist,p2list)
plt.xlabel ("$t$")
plt.ylabel("$p$")
plt.show()   
#Phasenraumtrajektorie
plt.plot(p1list,p2list)
plt.xlabel("$p1$")
plt.ylabel("$p2$")
plt.show()