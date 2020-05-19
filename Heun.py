# -*- coding: utf-8 -*-


import matplotlib.pyplot as plt
import numpy as np
import math 

# Parameter definieren
w=1 #Parameter der rechten Handseite
xpunkt= 1 #Startwert der ersten Ableitung der zuintegriernden Funktion
x = 0 #Startwert der Funktion
rhs1 = (-w**2) #Rechtehandseite des ersten Integrationsschritts

h= 0.001 #Schrittweite
t = 0 # Startpunkt der Zeit

#Euler
tlist = [t]
xlist = [x]
xpunktlist = [xpunkt]
while t <= 20* math.pi:
    t+=h
    tlist.append(t)
    xpunkt += rhs1*x
    xpunktlist.append(xpunkt)
    x+= xpunkt
    xlist.append(x)
    
#Diagramm
plt.plot(tlist,xlist)    
plt.xlabel("$t$")
plt.ylabel("$x$")
plt.show()


#Heun
#neue Parameter benennnen
y=0
ypunkt=1
t2=0
i=1

t2list =[t2]
ypunktlist = [ypunkt]
ylist = [y]
while t2 <= 20* math.pi:
    t2+=h
    t2list.append(t2)
    ypunkt+=(h/2)*(rhs1*y+rhs1*xpunktlist[i]) #Rückgriff auf Euler von oben
    ypunktlist.append(ypunkt)
    y += (h/2)*(ypunkt + xpunktlist[i])
    ylist.append(y)
    i+=1
 
#Diagramm
  
plt.plot(t2list,ylist)    
plt.xlabel("$t$")
plt.ylabel("$y$")
plt.show()  