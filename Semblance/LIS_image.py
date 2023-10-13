# Autor: Jorge Uriel Perez Romero
# Fecha de inicio del programa: 13 de junio de 2022
# Fecha de final:
# Descripcion de lo que hace el prgrama:
# Este programa calcula los diferentes espectros, los cuales son:
# Vos y Marius
# Maurin
# W.R. Webber and J.A. Lockwood
# Burguer
# Della Torre
# Garcia Mugnoz
# Moskalenko
# Exportando bibliotecas
import numpy as np
import math
import matplotlib.pyplot as plt
from sympy import sympify
# *Para este caso el vector de rigidez es el mismo que el de la energia cinetica*
# -----------------------------------------------------------------------------------------------
# LIS de Vos y Potgieter (2015)
def LISVP(Ek):
    c = 299792.458  # [m/s]
    v = 299592.458  # [Velocidad de la particula]
    beta = v/c
    JVp = []
    for i in range(len(Ek)):
        jvp=2.7*((Ek[i]**1.12)/beta**2)*(((Ek[i]+0.67)/1.67)**-3.93)
        JVp.append(jvp)
    Vp = np.array(JVp)
    return 1000*Vp #GeV
# -----------------------------------------------------------------------------------------------
# LIS de Maurin Protones
def LISJMP(Ek):
    JMH = []
    ciH = [3.4617, -4.1310, -4.6403, -1.4058, -4.7537, 8.5077, 32.637, -28.383, -58.203, 48.129, 33.946, -29.586, 0.61683]
    for i in range(len(Ek)):
        jmH=0
        if Ek[i] < 800:
            for j in range(len(ciH)):
                jmH += ciH[j]*((math.log(Ek[i], 10))/(math.log(800, 10)))**j
        else:
            jmH = -3.7995-2.7040*(math.log(Ek[i]/800, 10))
        JMH.append(jmH)
    JMP=np.array(JMH)
    return (Ek**2)*(10**JMP)
# -----------------------------------------------------------------------------------------------
# LIS de Maurin Helio
def LISJMHe(Ek):
    JMH = np.zeros(len(Ek))
    ciH = (2.2784, -4.5746, -4.865, 0.39567, -1.1578, 4.9893, 16.511, -20.521, -28.367, 31.850, 15, -17.083, 0.60486)
    for i in range(len(Ek)):
        jmH=0
        if Ek[i] <= 800:
            for j in range(len(ciH)):
                JMH[i] += ciH[j]*((np.log10(Ek[i]))/(np.log10(800)))**j
        else:
            JMH[i] = -4.9261-2.7140*(np.log10(Ek[i]/800))
    return 2*(Ek**2)*(10**JMH)
# -----------------------------------------------------------------------------------------------
# LIS de Burguer Protones
def LISB(Ek):
    JBp = []
    for i in range(len(Ek)):
        Jbp=0.
        Jbp=(1.9E+4*Ek[i]**-2.78)/(1+0.4866*Ek[i]**-2.51)
        JBp.append(Jbp)
    JB = np.array(JBp)
    return JB
# -----------------------------------------------------------------------------------------------
# LIS de Burguer Helio
def LISBHe(Ek):
    JBp = []
    for i in range(len(Ek)):
        Jbp=0.
        Jbp=(3.8E+4*Ek[i]**-2.78)/(1+0.9732*Ek[i]**-2.51)
        JBp.append(Jbp)
    JB = np.array(JBp)
    return JB
# -----------------------------------------------------------------------------------------------
# LIS Della Torre para Protones
def LISDT(P):
    JDtp = []
    aip = [94.1 -831, 0, 16700, -10200, 0]
    aip.sort(reverse=True)
    b = 10800
    c = 8590
    d1 = -4230000
    d2 = 3190
    c1 = 274000
    c2 = 17.4
    f1 = -39400
    f2 = 0.464
    g = 0
    for i in range(len(P)):
        jdtp=0.
        if  P[i] <= 1.:
            for j in range(len(aip)):
                jdtp+=(aip[j]*P[i]**j)/(P[i]**-0.7)
        else: 
            jdtp=(b+(c/(P[i]))+(d1/(d2+P[i]))+(c1/(c2+P[i]))+(f1/(f2+P[i]))+(g*P[i]))/P[i]**2.7
        JDtp.append(jdtp)
    JDT=np.array(JDtp)
    return JDT

# -----------------------------------------------------------------------------------------------
# LIS Della Torre para Helio
def LISDTHe(P):
    JDtp = []
    aip = [1.14, 0, -118, 578, 0, -87]
    aip.sort(reverse=True)
    b = 3120
    c = -5530
    d1 = 3370
    d2 = 1.29
    c1 = 134000
    c2 = 88.5
    f1 = -1170000
    f2 = 861
    g = 0.03
    for i in range(len(P)):
        jdtp = 0.
        if  P[i] <= 1.5:
            for j in range(len(aip)):
                jdtp+=(aip[j]*P[i]**j)/(P[i]**-0.7)
        else:
            jdtp=(b+(c/(P[i]))+(d1/(d2+P[i]))+(c1/(c2+P[i]))+(f1/(f2+P[i]))+(g*P[i]))/P[i]**2.7
        JDtp.append(jdtp)
    JDT=np.array(JDtp)
    return JDT
# -----------------------------------------------------------------------------------------------
# LIS de Garcia-Mugnoz para protones
def LISJGM(Ek):
    JGmp = []
    for i in range(len(Ek)):
        jgmp = 9.9E+8*(Ek[i]+780*math.exp(-2.5E-4*Ek[i]))**-2.65
        JGmp.append(jgmp)
    JGM=np.array(JGmp)
    return JGM
# -----------------------------------------------------------------------------------------------
# LIS Garcia-Mugnoz para He
def LISJGMHe(Ek):
    JGmHe = []
    for i in range(len(Ek)):
        jgmHe = 1.8E+8*(Ek[i]+660*math.exp(-1.4E-4*Ek[i]))**-2.77
        JGmHe.append(jgmHe)
    JGM = np.array(JGmHe)
    return 1000*JGM
# -----------------------------------------------------------------------------------------------
# LIS Moskalenko
def LISJMK(Ek):
    JMk = []
    for i in range(len(Ek)):
        jmk=0.
        if Ek[i] <= 1.:
            jmk = np.exp(4.64 -0.08*np.emath.log(Ek[i])**2 -2.91*(Ek[i])**(1/2))
        else:
            jmk = np.exp(3.22 -2.86* np.emath.log(Ek[i]) -1.5/Ek[i])
        JMk.append(jmk)
    JMK = np.array(JMk)
    return 1000*JMK
# -----------------------------------------------------------------------------
# LIS de Webber & Lookwood de Protones
def LISWL(Ek):
    JWlp = []
    for i in range(len(Ek)):
        jwlp = (2.1*Ek[i]**-2.8)/(1+5.85*Ek[i]**-1.22+1.18*Ek[i]**-2.54)
        JWlp.append(jwlp)
    JWL=np.array(JWlp)
    return 1000*JWL

# -----------------------------------------------------------------------------
# LIS de Webber & Lookwood de Helio
def LISWLHe(Ek):
    JWlh = []
    for i in range(len(Ek)):
        jwlh = (1.075*Ek[i]**-2.8)/(1+3.9*Ek[i]**-1.09+0.9*Ek[i]**-2.54)
        JWlh.append(jwlh)
    JWL=np.array(JWlh)
    return 1000*JWL

# Haciendo graficas      
# Declarando el vector de la energía cinetica
Ek1 = np.logspace(-1, 2, 25) 
Ek2 = np.logspace(0, 3, 25) 
Ek3 = np.logspace(0, 1, 25) 
Ek4 = np.logspace(-2, 5, 25) 
Ek5 = np.logspace(1, 4, 25) 
Ek6 = np.logspace(-1, 3, 25) 
Ek7 = np.logspace(-2, 2, 25) 
Vos = plt.loglog(Ek1, LISVP(Ek1), linewidth = 1.5, label = 'Vos and Potgieter LIS in 2015')     # Vos
MP = plt.loglog(Ek2, LISJMP(Ek2), linewidth = 1.5, label = 'Ghelfi, Barao, Derome and Maurin Proton LIS in 2017')    # Maurin Protones
MH = plt.loglog(Ek2, LISJMHe(Ek2), linewidth = 1.5, label = 'Ghelfi, Barao, Derome and Maurin Helium LIS in 2017')   # Maurin Helio
BP = plt.loglog(Ek3, LISB(Ek3), linewidth = 1.5, label = 'Burguer and Potgieter Proton LIS in 2000')       # Burguer Protnes
BH = plt.loglog(Ek3, LISBHe(Ek3), linewidth = 1.5, label = 'Burguer and Potgieter Helium LIS in 2000')    # Burguer Helio
DP = plt.loglog(Ek4, LISDT(Ek4), linewidth = 1.5, label = 'Boschini, Della Torre, Gervasi and more Proton LIS in 2018')     ## Della Torre Protones  
DH = plt.loglog(Ek4, LISDTHe(Ek4), linewidth = 1.5, label = 'Boschini, Della Torre, Gervasi and more Helium LIS in 2018')   ## Della Torre Helio   
GMP = plt.loglog(Ek5, LISJGM(Ek5), linewidth = 1.5, label = 'Garcia-Munoz, Mason and Simpson Proton LIS in 1975')    # Garcia Muñoz Protones
GMH = plt.loglog(Ek5, LISJGMHe(Ek5), linewidth = 1.5, label = 'Garcia-Munoz, Mason and Simpson Helium LIS in 1975')  # Garcia Muñoz Helio
M = plt.loglog(Ek6, LISJMK(Ek6), linewidth = 1.5, label = 'Moskalenko, Strong, Ormes and Potgieter LIS in 2002')    # Moskalenko 
WLP = plt.loglog(Ek7, LISWL(Ek7), linewidth = 1.5, label = 'Lagner, Potgieter and Webber Proton LIS in 2003')     # Webber & Lookwood Protones
WLH = plt.loglog(Ek7, LISWLHe(Ek7), linewidth = 1.5, label = 'Lagner, Potgieter and Webber Helium LIS in 2003')   # Webber & Lookwood He
plt.title('Local Interstellar Spectra', fontsize = 15)
plt.legend(loc = 'center left', bbox_to_anchor = (1, 0.5))
plt.ylim(1E-3, 4E+6)
plt.xlabel('Kinetic energy [Gev/n]')
plt.xlim(1E-2, 1E+4)
plt.ylabel('Flux particle [1/sr  1/m^2 1/s GeV/n]')
plt.grid()
plt.savefig("LIS", bbox_inches = 'tight')
plt.show()
