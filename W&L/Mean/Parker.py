import numpy as np
import math
import matplotlib.pyplot as plt
# Este programa solo sirve para graficar
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
Station('PhiCD-CM-HRMS.csv', )