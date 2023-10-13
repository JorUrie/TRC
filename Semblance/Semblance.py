import os
import numpy as np
import pandas as pd
from scipy import signal
import matplotlib.pyplot as plt
##
os.system('cls')
def Semblance(y1, y2, nscales):
    # Demo analysis
    if nscales == 0:
        t = np.linspace(0, 512, 513)
        y1 = np.sin(t*0.05) + np.sin(t*0.15)
        y2 = np.zeros(len(t))
        for i in range(len(t)):
            if i <= 170:
                y2[i] = -np.sin(t[i]*0.05) + np.sin(t[i]*0.15)
            if i > 170:
                y2[i] = np.sin(t[i]*0.05) + np.cos(t[i]*0.15)
            if i >= 341:
                y2[i] = np.cos(t[i]*0.05) - np.sin(t[i]*0.15)
        nscales = 100

    # semblance start
    np.isnan(y1)
    np.isnan(y2)
    m1 = np.mean(y1)
    m2 = np.mean(y2)
    y1 = y1-m1
    y2 = y2-m2
    nscales = round(abs(nscales))
    scales = np.arange(1, nscales)
    c1 = signal.cwt(y1, signal.morlet2, scales)
    c2 = signal.cwt(y2, signal.morlet2, scales)
    ctc = c1*np.conj(c2)
    spt = np.arctan2(ctc.imag, ctc.real)
    s = np.cos(spt)
    # Plot images
    plt.subplot(3, 1, 1)
    plt.plot(y1)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Signal I')
    plt.subplot(3, 1, 2)
    plt.plot(y2)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Signal II')
    plt.subplot(3, 1, 3)
    sz = np.flipud(s)
    fig = plt.plot(sz)
    im4 = plt.imshow(sz)
    plt.colorbar()
    im4.set_cmap('jet')
    plt.title("Semblance")
    plt.show()
       
# Cargando datos de los espectros 
def Aufladung(Dokumente, columne):
    # --------------------
    # Llamando los archivos que contienen los datos
    e1 = pd.read_csv(Dokumente, header = columne) 
    # Capturando datos
    E1E = pd.DataFrame(e1)
    # Transformadolo a lista
    E1L = E1E.values.tolist()
    #Generando valores zero para que sea un arreglo mas acomodado
    E1A = np.zeros(len(E1L))
    # Convirtiendo en arreglos
    E1 = np.array(E1L)
    for i in range(len(E1)):
        E1A[i] = E1[i]
    return E1A
# --------------------
# Cargando datos de los espectros 
def Ladung(Dokumente, columne):
    # --------------------
    # Llamando los archivos que contienen los datos
    e1 = pd.read_csv(Dokumente, header = 1, usecols = [columne]) 
    # Capturando datos
    E1E = pd.DataFrame(e1)
    # Transformadolo a lista
    E1L = E1E.values.tolist()
    #Generando valores zero para que sea un arreglo mas acomodado
    E1A = np.zeros(len(E1L))
    # Convirtiendo en arreglos
    E1 = np.array(E1L)
    for i in range(len(E1)):
        E1A[i] = E1[i]
    return E1A
# --------------------
# Cargando los datos
s1 = Aufladung('PhiCD-CM-AATB.csv', int(0))
s2 = Aufladung('PhiCF-CM-AATB.csv', int(0)) 
MF = Ladung('Data-CM.csv', 2)
# --------------------
ms1 = Aufladung('PhiCD-MS-AATB.csv', int(0))
ms2 = Aufladung('PhiCF-MS-AATB.csv', int(0)) 
SS = Ladung('Data-MS.csv', 2)
# --------------------
# Haciendo la semblanza entre las dos señales
# Campo magnetico
sCM1 = Semblance(s1, MF, 150)
sCM2 = Semblance(s2, MF, 150)
# --------------------
# Manchas solares
sMS1 = Semblance(ms1, SS, 150)
sMS2 = Semblance(ms2, SS, 150)
