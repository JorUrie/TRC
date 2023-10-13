# Autores: Jorge Uriel Perez Romero
#          Miguel E. Vargas
# Fecha de inicio del programa: 13 de junio de 2022
# Fecha de final:
# Descripcion de lo que hace el prgrama:
# Este programa funciona los diferentes espectros, los cuales son:
# Vos y Marius
# Maurin
# W.R. Webber and J.A. Lockwood
# Burguer
# Della Torre
# Garcia Mugnoz
# Moskalenko
# Exportando bibliotecas

import os
import math
import numpy as np
from ast import Del
import pandas as pd
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
# --------------------------------------------------------------------------------------------------
# Limpiando pantalla
os.system ("cls") 
# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------
# LIS de Maurin Protones 
def LISMP(Ek):
    ciH = [3.4617, -4.1310, -4.6403, -1.4058, -4.7537, 8.5077, 32.637, -28.383, -58.203, 48.129, 33.946, -29.586, 0.61683]
    jmH=0
    if Ek<800:
        for j in range(len(ciH)):        
            jmH += ciH[j]*((math.log(Ek, 10))/(math.log(800, 10)))**j
    else:
        jmH = -3.7995-2.7040*(math.log(Ek/800, 10))
    return (Ek**2)*(10**jmH)
# --------------------------------------------------------------------------------------------------
# LIS de Maurin Helio 
def LISMHe(Ek):
    ciH = (2.2784, -4.5746, -4.865, 0.39567, -1.1578, 4.9893, 16.511, -20.521, -28.367, 31.850, 15, -17.083, 0.60486)
    jmH=0
    if Ek <= 800:
        for j in range(len(ciH)):
            jmH += ciH[j]*((math.log(Ek, 10))/(math.log(800, 10)))**j
        else:
            jmH = -4.9261-2.7140*(math.log(Ek/800, 10))
    return 2*(Ek**2)*(10**jmH)

# -------------------------------------------------------------------------------------------------
# LIS de Burguer Protones 
def LISBP(Ek):
    Jbp = (1.9E+4*Ek**-2.78)/(1+0.4866*Ek**-2.51)
    return Jbp
# --------------------------------------------------------------------------------------------------
# LIS de Burguer Helio
def LISBHe(Ek):
    Jbp=(3.8E+4*Ek**-2.78)/(1+0.9732*Ek**-2.51)
    return Jbp
# -------------------------------------------------------------------------------------------------
# LIS Della Torre para Protones
def LISDTP(P):
    aip = [94.1 -831, 0, 16700, -10200, 0]
    aip.sort(reverse=True)
    b=10800
    c=8590
    d1=-4230000
    d2=3190
    c1=274000
    c2=17.4
    f1=-39400
    f2= 0.464
    g=0
    JDtp = 0.
    if  P <= 1.:
        for j in range(len(aip)):
            JDtp+=(aip[j]*P**j)/(P**-0.7)
    else: 
        JDtp=(b+(c/(P))+(d1/(d2+P))+(c1/(c2+P))+(f1/(f2+P))+(g*P))/P**2.7
    return JDtp
# -----------------------------------------------------------------------------------------------
# LIS Della Torre para Helio
def LISDTHe(P):
    aip = [1.14, 0, -118, 578, 0, -87]
    aip.sort(reverse=True)
    b=3120
    c=-5530
    d1=3370
    d2=1.29
    c1= 134000
    c2= 88.5
    f1=-1170000
    f2= 861
    g=0.03
    JDtp = 0.
    if  P <= 1:
        for j in range(len(aip)):
            JDtp+=(aip[j]*P**j)/(P**-0.7)
    else:
        JDtp=(b+(c/(P))+(d1/(d2+P))+(c1/(c2+P))+(f1/(f2+P))+(g*P))/P**2.7
    return JDtp
# -----------------------------------------------------------------------------------------------
# LIS de Garcia-Mugnoz para protones
def LISJGMP(Ek):
    jgmp = 9.9E+8*(Ek+780*math.exp(-2.5E-4*Ek))**-2.65
    return 1000*jgmp
# -----------------------------------------------------------------------------------------------
# LIS Garcia-Mugnoz para He 
def LISJGMHe(Ek):
    jgmHe = 1.8E+8*(Ek+660*math.exp(-1.4E-4*Ek))**-2.77
    return 1000*jgmHe
# -----------------------------------------------------------------------------------------------
# LIS de Webber & Lookwood de Protones
def LISWLP(Ek):
    jwlp = (2.1*Ek**-2.8)/(1+5.85*Ek**-1.22+1.18*Ek**-2.54)
    return 1000*jwlp
# --------------------------------------------------------------------------------------------------
# LIS de Webber & Lookwood de Helio 
def LISWLHe(Ek):
    jwlh = (1.075*Ek**-2.8)/(1+3.9*Ek**-1.09+0.9*Ek**-2.54)
    return 1000*jwlh

#+-+-+-+-+-+--+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#+-+-+-+-+-+--+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#+-+-+-+-+-+--+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
def ConDifP(Data, phi):
    E0 = 0.938          # Energia de reposo del proton en GeV
    Az = 1.0            # Numero de nucleones/ numero de protones
    # Estableciendo condiciones iniciales
    # Generando vectores
    Pb = Data
    Betab = Pb/(Pb*Pb+(Az*E0)**2)**(0.5)
    Tb = Pb/Betab/Az-E0
    LIS = LISWLP(Pb) #<------- Se cambia el espetro que lo calcula dato por dato
    J = LIS*math.exp(-3.*phi/(Betab*Pb)) 
    return Tb, J, LIS, Betab
# --------------------------------------------------------------------------------------------------
# Algoritmo de Difusion-Conveccion
def ConDifHe(Data, phi):
    E0 = 0.938            # Energia de reposo del proton en GeV
    Az = 2.0              # Numero de nucleones/ numero de protones
    # Estableciendo condiciones iniciales
    # Generando vectores
    Pb = Data
    Betab = Pb/(Pb*Pb+(Az*E0)**2)**(0.5)
    Tb = Pb/Betab/Az-E0
    LIS = LISWLHe(Pb) #<--------------- Se cambia el espectro
    J = LIS*math.exp(-3.*phi/(Betab*Pb)) 
    return Tb, J, LIS, Betab

# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------
# Funcion que calcula las variaciones de modulacion
# Debe recibir la rigidez de corte de la estacion
# Debe recibir los datos del numero de cuentas
# -----------------------------------------------------------------------------------------
import math
import numpy as np
import matplotlib.pyplot as plt
# Funcion que calcula las variaciones de modulacion
# Debe recibir la rigidez de corte de la estacion
# Debe recibir los datos del numero de cuentas
def DphiCD(dat, pc):
    # pc es la funcion de rigidez de corte
    # Numero de datos
    numdat = len(dat)
    # Numero de muestras
    # Es independiente del numero de muestras de los datos. El programa es muy sensible a este valor
    n = numdat
    # Rigidez de corte
    pmin = 0.1
    # Teoricamente debe ser infinito
    rcf = 10000
    # Contstantes
    k = 0.952
    alpha = 10.068
    N01987 = 148.1
    # Número de nucleones para hidrogeno y para helio
    AZh = 1
    AZhe = 2
    # Energía en reposo del electrón [Gev]
    E0 = 0.938
    #-----------Constantes Caballero-Moral-----------
    # Hidrogeno
    F0h = 14000
    p0h = 0.96
    ah = 1.5
    gam1h = -2.7
    gam2h = 2.8

    # Helio
    F0he = 2500
    p0he = 1.5
    ahe = 1.6
    gam1he = -2.7
    gam2he = 3.0
    # Constantes para el cálculo de la razon de las funciones de producción JH y JHe
    F0 = 2.0
    p0 = 0.45
    a = 1.4
    gam1 = 0
    gam2 = 10.0
    
    F0sh = 0.0002
    p0sh = 2.0
    ash = 1.45
    gam1sh = 0.89
    gam2sh = 7.7
    # Contador i correspondiente a la rigidez de corte pc
    # Redondea con una determinada cantidad de decimales
    ipc = round((math.log(pc,10)-math.log(pmin,10))*n/math.log(rcf,10))
    # Se definen arreglos
    P = [0]*n
    dNdP = [0]*n
    FP = [0]*n
    jh = [0]*n
    jhe = [0]*n
    sh = [0]*n
    sh1 = [0]*n
    she = [0]*n
    GPh = [0]*n
    # Del hidrogeno y del helio
    pmh = [0]*n
    pmhe = [0]*n
    jth = [0]*n
    jthe = [0]*n
    dtdph = [0]*n
    dtdphe = [0]*n
    jpph = [0]*n
    jpphe = [0]*n
    deltphi = [0]*n
    Nn = [0]*n
    LL = [0]*n
    Beta = [0]*n
    # Se define un vector de rigidez
    for i in range (1,n+1):
        x = i
        y = n
        P[i-1]=pmin*rcf**(x/y)
    # Es la funcion de Droman. Los valores de 10.068 y 0.952 son validos para cuando el Sol se encuentra
    # en la minima actividad
    datmin = 1-math.e**(-10.068*pc**(-0.952))
    # Se inicializa la variable delta phi
    delphi = 1
    # Condiciones iniciales
    for i in range(1, n+1):
        # dN/dP es la funcion de respuesta tambien conocida como conteo diferencial 
        dNdP[i-1]=-alpha*k*math.e**(-alpha*P[i-1]**(-k))*P[i-1]**(-k-1)
        # Razon de las funciones de produccion de particulas en general de Clem y Dorman
        FP[i-1]=F0*(p0**a+P[i-1]**a)**((gam1-gam2)/a)*P[i-1]**gam2
        # Funcion de produccion de Hidrogeno(Protones)
        jh[i-1]= LISWLP(P[i-1])#F0h*(p0h**ah+P[i-1]**ah)**((gam1h-gam2h)/ah)*P[i-1]**gam2h
        # Funcion de produccion para Helio 
        jhe[i-1]= LISWLHe(P[i-1])#F0he*(p0he**ahe+P[i-1]**ahe)**((gam1he-gam2he)/ahe)*P[i-1]**gam2he
        # Funcion de produccion desarrollada por Caballero-Lopez y Moraal
        sh[i-1]=-dNdP[i-1]/(jh[i-1]+1.584*FP[i-1]*jhe[i-1])
        # Funcion de produccion para el Helio
        sh1[i-1]=F0sh*(p0sh**ash+P[i-1]**ash)**((gam1sh-gam2sh)/ash)*P[i-1]**gam2sh
        # Diferencia entre la duncion de produccion y funcion de produccion desarrollada por Caballero-Lopez
        she[i-1]=FP[i-1]*sh[i-1]
        # Relacion entre la funcion de produccion de Hidrogeno y Helio
        GPh[i-1]=sh[i-1]/sh1[i-1]
    for k in range(1, numdat+1):
        for j in range(1,n+1): #<----------------- OJO!!!!!!!!!!
            jx=j        
            # Condicion para valores muy bajo de energia cinetica        
            if (dat[k-1]<0.1):
                break
            # Para cuando son muy grnades los valores (Los normaliza)
            if (dat[k-1]>100):
                delphi=(-jx)/100
            else:
                delphi=(jx-1)/100
            # Es el algoritmo de Conveccion-Difusion para el Hidrogeno como para el Helio
            for ii in range(n+1):
                pp=[P[ii-1]]
                jth[ii-1]=ConDifP(P[ii-1], delphi)[1]
                jthe[ii-1]=ConDifHe(P[ii-1], delphi)[1]
                jpph=jth#*dtdph[i-1]
                jpphe=jthe#*dtdphe[i-1] 
            # Iniciando vectores dandole valor de cero
            for i in range(1,n+1):
                Nn[i-1]=0
                LL[i-1]=0
                # Parametro de modulacion
            Nn[n-1]=rcf**(1+gam1sh+gam1h)*(F0sh)*(F0h+1.584*F0*F0he)/abs(1+gam1sh+gam1h)
            for i in range(n-2,1,-1):
                LL[i-1]=(sh[i-2]*(jpph[i-2]+1.584*FP[i-2]*jpphe[i-2])+sh[i-1]*(jpph[i-1]+1.584*FP[i-1]*jpphe[i-1]))*(P[i-1]-P[i-2])/2
                Nn[i-1]=Nn[i]+LL[i-1]
            # Variacion en el parametro de modulacion para un determinado tiempo
            print(Nn[ipc-1]/datmin,dat[k-1]/100,  (Nn[ipc-1]/datmin)-(dat[k-1]/100),delphi)
            if (abs((Nn[ipc-1]/datmin)-(dat[k-1]/100))<=0.01):
                break
                # Ya suelta deltaphi
            deltphi[k-1]=delphi
    i = 0
    Deltphi = []
    for i in range(len(deltphi)):
        if deltphi[i] >= 6:
            Deltphi.append(0)
        else:
            Deltphi.append(deltphi[i])
    return Deltphi, jpphe, dat, P, Beta
# jh: Espectro
# dat: Datos del monitor de neutrones
# P: Rigidez
# --------------------------------------------------------------------------------------------------
def Station(Dokumente, Date, Name):
    # Llamando los archivos que contienen los datos
    date = pd.read_csv(Dokumente, header = 0, usecols = [Date])
    e1 = pd.read_csv(Dokumente, header=0, usecols = [Name])
    # Capturando datos
    DateE = pd.DataFrame(date)
    E1E = pd.DataFrame(e1)
    # Transformadolo a lista
    DateL = DateE.values.tolist()
    E1L = E1E.values.tolist()
    #Generando valores zero para que sea un arreglo mas acomodado
    DateA = np.zeros(len(date))
    E1A = np.zeros(len(date))
    # Convirtiendo en arreglos
    Date = np.array(DateL)
    E1 = np.array(E1L)
    for i in range(len(Date)):
        DateA[i] = Date[i]
        E1A[i] = E1[i]
    # Seleccionando valores maximos de cada arreglo para usar valores relativos
    E1max = np.amax(E1)
    # Generando los vectores de valores relativos
    E1R = (E1A/E1max)*100
    return E1R, DateA
# --------------------------------------------------------------------------------------------------
def Plotter(Dokumente, Date):
    # --------------------
    # Llamando los archivos que contienen los datos
    date = pd.read_csv(Date, header = 1, usecols = [0]) # Jahr
    date1 = pd.read_csv(Date, header = 1, usecols = [1]) # Tag
    e1 = pd.read_csv(Dokumente) 
    # Capturando datos
    DateE = pd.DataFrame(date)
    DateE1 = pd.DataFrame(date1)
    E1E = pd.DataFrame(e1)
    # Transformadolo a lista
    DateL = DateE.values.tolist()
    DateL1 = DateE1.values.tolist()
    E1L = E1E.values.tolist()
    #Generando valores zero para que sea un arreglo mas acomodado
    DateA = []
    DateB = []
    E1A = np.zeros(len(date))
    # Convirtiendo en arreglos
    Date = np.array(DateL)
    Date1 = np.array(DateL1)
    E1 = np.array(E1L)
    for i in range(len(Date)):
        DateA.extend([str(Date1[i]) + '/' + str(Date[i])])
        E1A[i] = E1[i]
    # Dando formato a las fechas
    chars = '[, ]'
    DateB = []
    for j in range(len(DateA)):
        DateB.append(DateA[j].translate(str.maketrans('', '', chars)))
    # Convirtiendo a formato de fechas
    x = [dt.datetime.strptime(d,'%m/%Y').date() for d in DateB]
    # ---------------
    plt.rcParams['axes.facecolor'] = 'white'
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m/%Y'))
    plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval = 2550))
    plt.plot(x, E1A, '.', linewidth = 1, color = 'darkred')
    plt.gcf().autofmt_xdate()
    plt.title('Modulation Factor Data with \n Force Field solution', fontsize = 18)
    plt.xlabel('Date', fontsize = 16)
    plt.ylabel('[Counts/s]', fontsize = 16)
    plt.ylim(0, 1.5)
    f1 = dt.datetime(1989, 12, 1)
    f2 = dt.datetime(2016, 12, 1)
    plt.xlim(f1, f2)
    plt.grid()
    plt.savefig('Ausser', bbox_inches = 'tight')
    plt.show()  
# --------------------------------------------------------------------------------------------------
# Programa principal
Dok = 'AATB-MS.csv'
Date = 'Year'
Name = 'AATB'
Ek = Station(Dok, Date, Name)[0]
# Conveccion-Difusion
PhiCD = DphiCD(Ek, 5.90)[0]
# Guardando los datos modulados y las fechas en un txt
np.savetxt('PhiCD-CM-AATB.csv', PhiCD, delimiter = ' ')
# Print spectrum
Plotter('PhiCD-CM-AATB.csv', Dok)