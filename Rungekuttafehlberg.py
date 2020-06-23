# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 18:05:48 2020

@author: Agnes
"""


import numpy as np
import matplotlib.pyplot as plt
import math
from mpl_toolkits import mplot3d

#Butcher-Tableau
alpha = np.array([0.0, 1.0/4.0, 3.0/8.0, 12.0/13.0, 1.0, 1.0/2.0])
beta = np.array([[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1.0/4.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [3.0/32.0, 9.0/32.0, 0.0, 0.0, 0.0, 0.0],
        [1932.0/2197.0, (-7200.0)/2197.0, 7296.0/2197.0, 0.0, 0.0, 0.0],
        [439.0/216.0, -8.0, 3680.0/513.0, (-845.0)/4104.0, 0.0, 0.0],
        [(-8.0)/27.0, 2.0, (-3544.0)/2565.0, 1859.0/4104.0, (-11.0)/40.0, 0.0]])
c = np.array([25.0/216.0, 0.0, 1408.0/2565.0, 2197.0/4104.0, (-1.0)/5.0, 0.0]) # coefficients for 4th order method
c_star = np.array([16.0/135.0, 0.0, 6656.0/12825.0, 28561.0/56430.0, (-9.0)/50.0, 2.0/55.0]) # coefficients for 5th order method
cerr=c-c_star # build the difference of both c and c_star for the error estimation


# Runge-Kutta-Fehlberg
def rkf45 (f,x,function_parameters, h):
    """f: Funktion
        x: letzter Zeitschrit
        function_parameters: Parameter der rechten Handseite
        h: Schrittweite"""
    k1 = f(x, *function_parameters)
    k2 = f(x + h*k1*beta[1][0], *function_parameters)
    k3 = f(x +h*(k1*beta[2][0]+k2*beta[2][1]), *function_parameters)
    k4 = f(x+ h*(k1*beta[3][0]+k2*beta[3][1]+k3*beta[3][2]), *function_parameters)
    k5 = f(x+ h*(k1*beta[4][0]+k2*beta[4][1]+k3*beta[4][2]+k4*beta[4][3]),  *function_parameters)
    k6 = f(x+ h*(k1*beta[5][0]+k2*beta[5][1]+k3*beta[5][2]+k4*beta[5][3]+k5*beta[5][4]), *function_parameters)
    xnew4 = x + h*(c[0]*k1+c[1]*k2+c[2]*k3+c[3]*k4+c[4]*k5+c[5]*k6)
    xnew5 = x + h*(c_star[0]*k1+ c_star[1]*k2+ c_star[2]*k3+ c_star[3]*k4+ c_star[4]*k5+ c_star[5]*k6)
    epsilon = abs(cerr[0]+k1 + cerr[1]*k2 +cerr[2]*k3 + cerr[3]*k4 + cerr[4]*k5 +cerr[5]*k6)  #Berechnung des Fehlers
    if epsilon>= epsilon_tol:       #Schätzung der neuen Schrittweite
        hnew = h*(epsilon_tol/epsilon)**0.2
    else:
        hnew = h* (epsilon_tol/epsilon)**0.25
    return (xnew4, xnew5,hnew,epsilon)

#Parameter
p= 0.1
q=1.5
r=1.5
# Startwerte
x =- math.sqrt(r-p) +0.0001     #Startwert x (proportional Eismasse)
y = math.sqrt(r-p)          #Startwert y (proportional CO2)
z= -math.sqrt(r-p)          #Satrtwert z (proportional Tiefsee-Temperatur)

epsilon_tol= 0.01 # tolerierter/optimaler Fehler
tx=0 #Startpunkt in Zeit
ty =0
tz=0
tmax=100 #Endpunkt in Zeit
h = 0.26 #Schrittweite
x4list=[x]
x5list = [x]
y4list= [y]
y5list= [y]
z4list = [z]
z5list = [z]
txlist=[tx]
tylist =[ty]
tzlist=[tz]
epsilonxlist = [0]
epsilonylist = [0]
epsilonzlist =[0]
#Durchführung
def rhs1(x,y):      #Rechtehandseite der x-Änderung
   return(-x-y)

def rhs2(y,p,z,r):          #Rechtehandseite der y-Änderung
    return(p*z+r*y-(z**2)*y)

def rhs3(z,q,x):
   return(q*(x-z))     #Rechtehandseite der z-Änderung

while tx <=tmax and ty<= tmax and tz<= tmax:

  
    #neues x,y,z ermitteln
    xnew=rkf45(rhs1,x,[y],h)   #welchen wert in reihe x4 oder x5???
    x4list.append(xnew[0])
    x5list.append(xnew[1])
    tx += xnew[2]
    txlist.append(tx)
    epsilonxlist.append(xnew[3])
    
    ynew=rkf45(rhs2,y,[p,z,r],h)
    y4list.append(ynew[0])
    y5list.append(ynew[1])
    ty += ynew[2]
    tylist.append(ty)
    epsilonylist.append(ynew[3])
    
    znew= rkf45(rhs3,z,[q,x],h)
    z4list.append(znew[0])
    z5list.append(znew[1])
    tz += znew[2]
    tzlist.append(tz)
    epsilonzlist.append(znew[3])
   
    #Variablen updaten
    x= xnew[0]*1
    y= ynew[0]*1
    z= znew[0]*1
 
    
    
#Diagramme
#Zeitreihe
plt.plot(txlist,x4list)
plt.plot(txlist,x5list)
plt.xlabel("$t$") 
plt.ylabel("$x$")
plt.show()
plt.plot(tylist,y4list)
plt.plot(tylist,y5list)
plt.xlabel("$t$") 
plt.ylabel("$y$")
plt.show()
plt.plot(tzlist,z4list)
plt.plot(tzlist,z5list)
plt.xlabel("$t$") 
plt.ylabel("$z$")
plt.show()

#Phasenraumtrajektorie
fig =plt.figure()
ax = fig.gca(projection='3d')
ax.plot(x4list, y4list, z4list)
ax.plot(x5list,y5list,z5list)
ax.legend()
ax.set_xlabel("$X$")
ax.set_ylabel("$y$")
ax.set_zlabel(r"$z$", fontsize=10, rotation=90)
plt.show()

#Fehler gegen Zeit
plt.plot(txlist, epsilonxlist)
plt.plot(tylist, epsilonylist)
plt.plot(tzlist,epsilonzlist)
plt.xlabel("$t$")
plt.ylabel("$Fehler$")
plt.show()