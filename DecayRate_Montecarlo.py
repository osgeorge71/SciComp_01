import matplotlib.pyplot as plt
from scipy import constants # type: ignore
import numpy as np
import random

class Funcion():
    #-Definición de constantes para la función*
    Me = 0.510998950 #Masa del electrón en MeV
    Mp = 938.272088  #Masa del protón en MeV
    Mn = 939.565420  #Masa del neutrón en MeV
    gw = 0.653 #Weak coupling constant
    Mw = 80.379 * 1000 #Masa del bosón W+- por 1000 para MeV
    hBar = 6.582119514e-22
    #-Cota máxima de la función, aumentada arbitrariamente*
    cotaMax = 0.0017

    def __init__(self):
        self.dif_Mnp = self.Mn - self.Mp  #Diferencia Mn - Mp
        self.factorAmplitud = (1/(np.power(np.pi,3)*self.hBar)) * np.power(self.gw/(2*self.Mw),4)

    def f(self, x):
        return self.factorAmplitud * x * np.sqrt(np.power(x,2)-np.power(self.Me,2))*np.power(self.dif_Mnp-x,2)
    
    def integralResueltaGriffiths(self):
        #Constantes que maneja Griffiths
        factor1 = 1/(4*constants.pi**3*6.58212e-22)
        factor2 = (self.gw / (2*80423))**4
        factor3 = 0.5109989**5
        factor4 = 6.54438 #Término del extremo derecho de la integral resuelta
        #Calcular el resultado de la integral
        resultadoGriffiths = factor1*factor2*factor3*factor4
        # print("factor1",factor1)
        # print("factor2",factor2)
        # print("factor3",factor3)
        # print("factor4",factor4)
        return resultadoGriffiths

    def integralResueltaCtesActuales(self):
        #Constantes actualizadas desde el PDG
        factor1 = 1/(4*constants.pi**3*self.hBar)
        factor2 = (self.gw / (2*self.Mw))**4
        factor3 = self.Me**5
        factor4 = 6.54438 #Término del extremo derecho de la integral resuelta
        #Calcular el resultado de la integral (con las constantes actualizadas)
        resultadoSoftware = factor1*factor2*factor3*factor4
        return resultadoSoftware
          
    def graficar(self, x1, x2, tamPaso):
        #Componer muestras en la abscisa
        E = np.arange(x1,x2,tamPaso) #Iterar en tamaño de paso
        #Componer muestras en la ordenada
        dFdE_f = self.f(E)
        #-Graficar*
        plt.plot(E,dFdE_f)
        plt.title("Electron energy distribution from neutron beta decay")
        plt.xlabel("E (MeV)")
        #charGamma = "\u0393"
        plt.ylabel("d\u0393/dE")
        plt.show()

class Montecarlo():
    def __init__(self):
        pass

    def calcIntegral_1(self,f,x1,x2,nSamples):
        sum = 0.0
        xr = np.random.uniform(x1,x2,nSamples)
        for i in range(nSamples):
            sum += f(xr[i])
        return (sum/nSamples)*(x2-x1)
    
    def calcIntegral_2(self,f,x1,x2,yMax,nSamples):
        success = 0.0
        xr = np.random.uniform(x1,x2,nSamples)
        yr = np.random.uniform(0.0,yMax,nSamples)
        for i in range(nSamples):
            if yr[i] <= f(xr[i]):
                success+=1
        return success/nSamples*(x2-x1)*yMax
    
    def graficar(self, cantMuestras, mc1, mc2, mcSegm1, mcSegm2, resSw):
        lst=[resSw]*len(cantMuestras)
        #-Graficar*
        plt.plot(cantMuestras,mc1,label="mc1")
        plt.plot(cantMuestras,mc2,label="mc2")
        plt.plot(cantMuestras,mcSegm1,label="mcSegm1")
        plt.plot(cantMuestras,mcSegm2,label="mcSegm2")
        plt.plot(cantMuestras,lst,label="valor objetivo")
        plt.title("Aproximación de valores con Montecarlo")
        plt.xlabel("Cantidad de muestras")
        plt.ylabel("Valores proporcionados por el método")
        plt.legend()
        plt.show()


def main():
    #Intervalo de trabajo para la función particular
    x1 = 0.511
    x2 = 1.29
    tamPaso = 1e-6
    funcion = Funcion() #Instanciar el objeto función
    # funcion.graficar(x1,x2,tamPaso)
    
    resultadoGriffiths = funcion.integralResueltaGriffiths()
    resultadoSoftware = funcion.integralResueltaCtesActuales()
    print("cálculo final del libro de Griffiths: ", resultadoGriffiths,", tao =",1.0/resultadoGriffiths)
    print("cálculo final constantes actuales: ", resultadoSoftware,", tao =",1.0/resultadoSoftware)
    
    #-Listas para almacenar resultados*
    lstMontecarlo_1 = []
    lstMontecarlo_2 = []
    lstMontecarloSegm_1 = []
    lstMontecarloSegm_2 = []
    #-Generación de las cantidades de muestras para realizar los cálculos*
    lstNumSamples = [i for i in range(10000,100000,20000)]
    for nSamples in lstNumSamples:
        #-----------------------------------------------------------------------
        print("Método de Montecarlo")
        #Llamado al método de la clase para calcular integral con Montecarlo
        mc = Montecarlo() #Instanciación de la clase en el objeto mc
        print("... con el método 1")
        Gamma = mc.calcIntegral_1(funcion.f,x1,x2,nSamples)
        print("Cálculo obtenido con Montecarlo método 1: ", Gamma)
        print("Tiempo de vida = 1/Gamma =", 1.0/Gamma, "contra",1.0/resultadoSoftware,"de la ref. con constantes actualizadas.")
        lstMontecarlo_1.append(1.0/Gamma)
        #-------*
        print("... con el método 2")
        Gamma = mc.calcIntegral_2(funcion.f,x1,x2,funcion.cotaMax,nSamples)
        print("Cálculo obtenido con Montecarlo método 2: ", Gamma)
        print("Tiempo de vida = 1/Gamma =", 1.0/Gamma, "contra",1.0/resultadoSoftware,"de la ref. con constantes actualizadas.")
        lstMontecarlo_2.append(1.0/Gamma)
        #-----------------------------------------------------------------------
        #-----------------------------------------------------------------------
        print("Método de Montecarlo segmentado")
        div = (x2-x1)/10.0  #Se eligen 10 sub-intervalos
        nSamples = int(nSamples/10); #Se divide el numero de muestras
        subInterv = [i for i in np.arange(0.511,1.29,div)]  
        #-------*
        print("... con el método 1")
        segm = [mc.calcIntegral_1(funcion.f,xx,xx+div,nSamples,) for xx in subInterv]
        Gamma = sum(segm)
        print("Tiempo de vida = 1/Gamma =", 1.0/Gamma, "contra",1.0/resultadoSoftware,"de la ref. con constantes actualizadas.")
        lstMontecarloSegm_1.append(1.0/Gamma)
        #-------*
        print("... con el método 2")
        segm = [mc.calcIntegral_2(funcion.f,xx,xx+div,funcion.cotaMax,nSamples,) for xx in subInterv]
        Gamma = sum(segm)
        print("Tiempo de vida = 1/Gamma =", 1.0/Gamma, "contra",1.0/resultadoSoftware,"de la ref. con constantes actualizadas.")
        lstMontecarloSegm_2.append(1.0/Gamma)
    #-Graficar*
    mc.graficar(lstNumSamples,lstMontecarlo_1,lstMontecarlo_2,lstMontecarloSegm_1,lstMontecarloSegm_2,1.0/resultadoSoftware)
    #-Guardar en un archivo para cálculos posteriores*
    # print("MC1")
    # print(lstMontecarlo_1)
    # print("MC2")
    # print(lstMontecarlo_2)
    # print("MCSegm1")
    # print(lstMontecarloSegm_1)
    # print("MCSegm2")
    # print(lstMontecarloSegm_2)
    # print("cantidades de muestras")
    # print(lstNumSamples)

    
if __name__ == "__main__":
    main()