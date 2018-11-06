import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.integrate import quad
from math import exp, log, sqrt, pi
import numpy as np

#Constants from experiment
RESONATOR_LENGTH = 10  # in cm
SPEED_OF_LIGHT = 3e8   # in m/s
MASS_RU_85 = 85        # in u
MASS_RU_87 = 87        # in u
ATOMIC_MASS = 931.494e6  # in eV
BOLTZMANN_CONSTANT = 8.615e-5 # in eV/K
LIGHT_WAVELENGTH = 780e-7   #in cm

# Loads the data in from the oszi format
def loadOsziData(filename, headerskip=18, cellskip=3):
    xdata = []
    ydata = []
    with open(filename) as f:
        for cnt, line in enumerate(f):
            if cnt >= headerskip:
                linesplit = line.split(',')
                xdata.append(float(linesplit[cellskip]))
                ydata.append(float(linesplit[cellskip+1]))
    return xdata, ydata

# At any point in xdata in where ydata is above the theshold, search around
# searchLength in xdata for maximum
def extractReferencePeakPosition(xdata, ydata, threshold=0.03, searchLength=0.001):
    peakTimes = []
    index = 0
    while index < len(xdata):
        if ydata[index] < threshold:
            index+=1
        else:
            maxAmplitude = ydata[index]
            maxIndex = index
            index0=index
            while index < len(xdata) and xdata[index] < xdata[index0] + searchLength:
                if ydata[index] > maxAmplitude:
                    maxAmplitude=ydata[index]
                    maxIndex = index
                index+=1
            peakTimes.append(xdata[maxIndex])
    return peakTimes

def frequenciesAtTime(resonatorFile, fitFunction, threshold=0.03, searchLength=0.001, showGraph=False):
    time, voltage = loadOsziData(resonatorFile)
    peakTimes = extractReferencePeakPosition(time, voltage)
    peakFrequencies = [i*1.0/(4*RESONATOR_LENGTH) for i in range(len(peakTimes))]
    params, conv = curve_fit(fitFunction,peakTimes, peakFrequencies)
    frequency = [fitFunction(t, *params) for t in time]

    if showGraph:
        plt.plot(peakTimes, peakFrequencies, 'bx')
        plt.plot(time, frequency, 'r-')
        plt.xlabel("time [s]")
        plt.ylabel("$\Delta \\nu$ [$\mathrm{cm}^{-1}$]")
        plt.figure()

    return time, frequency

def cutDataAtX(xdata, ydata, xMin, xMax):
    xdataCut = []
    ydataCut = []
    for i in range(len(xdata)):
        if xdata[i] >= xMin and xdata[i] <= xMax:
            xdataCut.append(xdata[i])
            ydataCut.append(ydata[i])
    return xdataCut, ydataCut

def fitAbsorbtionSpectrum(resonatorFile, absorbtionFile):
    time, frequency = frequenciesAtTime(resonatorFile, thirdOrderPoly, showGraph=True)
    _, intensity = loadOsziData(absorbtionFile)

    initialParams = [9e-3/2, 1.6e-1, 7.9e-2/2, 3.7e-1, 2.7e-1/2, 1.3e-1, 4.3e-1/2, 3.5e-2, 350, 65, -36, 27.59]
    initialIntensity = [absorbtionFit(x, *initialParams) for x in frequency]

    frequencyCut, intensityCut = cutDataAtX(frequency, intensity, -0.002, 0.50)

    fittedParams, cov = curve_fit(absorbtionFit, frequencyCut, intensityCut, p0 = initialParams)
    fittedIntensity = [absorbtionFit(x, *fittedParams) for x in frequency]
    err = np.sqrt(np.diag(cov))
    print(err)

    backgroundIntensity = [fittedParams[9]+fittedParams[10]*x+fittedParams[11]*x**2 for x in frequency]

    frequenncyGHz = [f*SPEED_OF_LIGHT*1e-7 for f in frequency]

    plt.plot(frequenncyGHz, intensity, 'bx', markersize = 1)
    #plt.plot(frequency, initialIntensity, 'r-')
    plt.plot(frequenncyGHz, fittedIntensity, 'r-')
    plt.plot(frequenncyGHz, backgroundIntensity, 'g-')
    plt.xlabel("$\Delta \\nu$ [GHz]")
    plt.ylabel("intensity [a.u.]")
    plt.legend(["Messdaten", "Fit", "Background"])

    print("Peak1: %f cm^-1 = %f GHz" % (fittedParams[0], SPEED_OF_LIGHT*fittedParams[0]*1e-7))
    print("Peak2: %f cm^-1 = %f GHz" % (fittedParams[2], SPEED_OF_LIGHT*fittedParams[2]*1e-7))
    print("Peak3: %f cm^-1 = %f GHz" % (fittedParams[4], SPEED_OF_LIGHT*fittedParams[4]*1e-7))
    print("Peak4: %f cm^-1 = %f GHz" % (fittedParams[6], SPEED_OF_LIGHT*fittedParams[6]*1e-7))
    print("Temperatur: %fK" % fittedParams[8])

    plt.show()

def fitSaturationSpectrum(resonatorFile, absorbtionFile):
    time, resonatorVoltage = loadOsziData(resonatorFile)
    peakTimes = extractReferencePeakPosition(time, resonatorVoltage)
    frequency = [(len(peakTimes)-1)*(t-peakTimes[0])/(peakTimes[len(peakTimes)-1]-peakTimes[0])/(4*RESONATOR_LENGTH)*SPEED_OF_LIGHT*1e-7 for t in time]

    _, intensity = loadOsziData(absorbtionFile)
    frequencyCut, intensityCut = cutDataAtX(frequency, intensity, 0.3, 1.2)

    #initialParams = [0.033, 0.08, 5e-4, 0.0435,0.08 ,5e-4, 600, 2e-3, 3.06, -1.33]
    #initialParams = [0.033, 5e-4, 0.043, 5e-4, 0.08, 0.08, 0, 0, 0, 350, 3.0, -1.33]
    #initialParams = [0.0305, 2e-4, 0.0345, 2e-4, 0.0427, 2e-4, 0.05,0.05,0.5, 0,0,0,0,0,0, 600, 2.8, -15]
    #initialParams = [0.65, 1.12, 600, 0.46,0.01,0.50,0.01,0.643,0.01,0,0,0.1,0.1,0.3,0.3, 3.5, 0]  # Peak2
    #initialParams = [0.013/2*SPEED_OF_LIGHT*1e-7, 2, 1500, 0.010/2*SPEED_OF_LIGHT*1e-7,1e-4,0.012/2*SPEED_OF_LIGHT*1e-7,1e-4,0.0156/2*SPEED_OF_LIGHT*1e-7,1e-4,0.01,0.01,0.1,0.1,0.3,0.3]  #Peak3
    #initialParams = [0.0305,4e-4,0.0345,4e-4, 0.0427,2e-3, 0.1, 0.1, 1.0,0.05,0.05,0.2,0.15,0.3,0.4,700, 2.3, -4]

    #initialParams = [0.0342/2*SPEED_OF_LIGHT*1e-7, 2.349, 1599, 0.17,1.14e-3,0.0172/2*SPEED_OF_LIGHT*1e-7,1.14e-3,0.0348/2*SPEED_OF_LIGHT*1e-7,1.14e-3,1.02e-6,2.09e-2,5.129e-2,1.92e-2,1.414e-1,2.451e-1, 3.9, -2.6*2/(SPEED_OF_LIGHT*1e-7)]  #Peak1
    #initialParams=[0.52, 2.4, 350,0.17, 0.02, 0.255, 0.01, 0.52, 0.01, 0.05, 0.02, 0.1, 0,0.1,0.1, 4, 0]  #Peak 1 new


    #initialParams = [0.0366, 2.55, 778, 0.033,7e-4,0.0354,1.27e-3,0.0428,1.04e-3,6e-2,2e-2,6.1e-2,2.1e-2,3.68e-2,9.38e-2, 3.0889, -2.16]  #Peak4
    initialParams = [0.736, 2.54, 334, 0.628, 0.0154, 0.75, 0.015, 0.86, 0.015, 0, 0, 0.05, 0.05, 0, 0.05, 3.07, -0.1]# Peak 4 new

    #initialParams = [1.80938000e-01,1.91527758e+00,350,1.47426891e-01,3.02025318e-03,1.76128042e-01,1.06632172e-02,2.37209167e-01,1.74246140e-02,1.75540675e-15,9.49667214e-02,1.04641174e-01,1.13356670e-01,8.51282898e-02,2.22552926e-01, 3, 0] #Peak3new

    initialIntensity = [saturationEmpiricGaussian(f, *initialParams) for f in frequencyCut]

    paramBounds = ([-np.inf]*17, [np.inf]*17)
    for i in range(9,15):
        paramBounds[0][i]=0

    #paramBounds[0][16] = -0.1
    #paramBounds[1][16] = 0.1

    fittedParams, cov = curve_fit(saturationEmpiricGaussian, frequencyCut, intensityCut, p0=initialParams, maxfev=20000, bounds=paramBounds)
    fittedIntensity = [saturationEmpiricGaussian(f, *fittedParams) for f in frequencyCut]

    print("params: ")
    print(fittedParams)

    print("error: ")
    err = np.sqrt(np.diag(cov))
    print(err)

    plt.plot(frequencyCut, intensityCut, 'b-')
    #plt.plot(frequencyCut, initialIntensity, 'r-')
    plt.plot(frequencyCut, fittedIntensity, 'r-')
    plt.xlabel("$\Delta\\nu$[GHz]")
    plt.ylabel("Intensitaet [a.u]")


    print("Peak1:  %f GHz" % fittedParams[3])
    print("Peak2: %f GHz" % fittedParams[5])
    print("Peak3: %f GHz" % fittedParams[7])

    plt.figure()
    plt.plot(frequencyCut, correctDopplerIntensities(frequencyCut, intensityCut, fittedParams), 'b-')
    plt.plot(frequencyCut, correctDopplerIntensities(frequencyCut, fittedIntensity, fittedParams), 'r-')
    plt.xlabel("$\Delta\\nu$[GHz]")
    plt.ylabel("Intensitaet [a.u]")

    plt.show()

# fit curves:

def correctDopplerIntensities(frequency,intensities, gaussianParams):
    intensityResult = []
    for i in range(len(intensities)):
        x = frequency[i]
        bg0 = gaussianParams[15]
        bg1 = gaussianParams[16]
        minVoltage = gaussianParams[1]
        xMid = gaussianParams[0]
        T = gaussianParams[2]
        newInt = intensities[i]-(bg0+bg1*x+(minVoltage-bg0-bg1*xMid)*dopplerVerbreiterung(x-xMid, T, MASS_RU_85, SPEED_OF_LIGHT/LIGHT_WAVELENGTH*1e-7))
        intensityResult.append(newInt)
    return intensityResult

def firstOrderPoly(x, a0, a1):
    return a0+a1*x

def thirdOrderPoly(x, a0, a1, a2, a3):
    return a0+a1*x+a2*x**2+a3*x**3

def dopplerVerbreiterung(x, T,m, nu0):
    return np.exp(-(x**2/(nu0**2*2*BOLTZMANN_CONSTANT*T/(m*ATOMIC_MASS))))

def absorbtionFit(x, x1, amp1, x2, amp2, x3, amp3, x4, amp4, T, bg0, bg1, bg2):
    intensity = 1
    intensity -= amp1*dopplerVerbreiterung(x-x1, T, MASS_RU_87, 1.0/LIGHT_WAVELENGTH)
    intensity -= amp2*dopplerVerbreiterung(x-x2, T, MASS_RU_85, 1.0/LIGHT_WAVELENGTH)
    intensity -= amp3*dopplerVerbreiterung(x-x3, T, MASS_RU_85, 1.0/LIGHT_WAVELENGTH)
    intensity -= amp4*dopplerVerbreiterung(x-x4, T, MASS_RU_87, 1.0/LIGHT_WAVELENGTH)
    intensity*=(bg0+bg1*x+bg2*x**2)
    return intensity

def saturationFit(x, x1, amp1, fwhm1, x2, amp2, fwhm2, x3, amp3, fwhm3, T, powerRatio, bg0, bg1):
    intensity = 1

    #Dopper
    intensity -= amp1*quad(absorbtionIntegrand, -np.inf, np.inf, (x, x1, fwhm1, T, MASS_RU_87*ATOMIC_MASS))[0]
    intensity -= amp2*quad(absorbtionIntegrand, -np.inf, np.inf, (x, x2, fwhm2, T, MASS_RU_87*ATOMIC_MASS))[0]
    intensity -= amp3*quad(absorbtionIntegrand, -np.inf, np.inf, (x, x3, fwhm3, T, MASS_RU_87*ATOMIC_MASS))[0]


    #Lamb peaks:
    intensity += powerRatio*amp1*amp1*quad(saturationIntegrand, -np.inf, np.inf, (x, x1, fwhm1, x1, fwhm1, T, MASS_RU_87*ATOMIC_MASS))[0]
    intensity += powerRatio*amp2*amp2*quad(saturationIntegrand, -np.inf, np.inf, (x, x2, fwhm2, x2, fwhm2, T, MASS_RU_87*ATOMIC_MASS))[0]
    intensity += powerRatio*amp3*amp3*quad(saturationIntegrand, -np.inf, np.inf, (x, x3, fwhm3, x3, fwhm3, T, MASS_RU_87*ATOMIC_MASS))[0]

    #Lamb peaks:
    intensity += 2*powerRatio*amp1*amp2*quad(saturationIntegrand, -np.inf, np.inf, (x, x1, fwhm1, x2, fwhm2, T, MASS_RU_87*ATOMIC_MASS))[0]
    intensity += 2*powerRatio*amp2*amp3*quad(saturationIntegrand, -np.inf, np.inf, (x, x2, fwhm2, x3, fwhm3, T, MASS_RU_87*ATOMIC_MASS))[0]
    intensity += 2*powerRatio*amp3*amp1*quad(saturationIntegrand, -np.inf, np.inf, (x, x3, fwhm3, x1, fwhm1, T, MASS_RU_87*ATOMIC_MASS))[0]

    intensity *= (bg0+x*bg1)
    return intensity

def saturationFitApprox(x, x1, amp1, fwhm1, x2, amp2, fwhm2, T, powerRatio, bg0, bg1):
    intensity = 1

    #Doppler
    intensity -= amp1*dopplerVerbreiterung(x-x1, T, MASS_RU_87, 1.0/LIGHT_WAVELENGTH)
    intensity -= amp2*dopplerVerbreiterung(x-x2, T, MASS_RU_87, 1.0/LIGHT_WAVELENGTH)
    #intensity -= amp3*dopplerVerbreiterung(x-x3, T, MASS_RU_87, 1.0/LIGHT_WAVELENGTH)

    #Lamb Peaks
    intensity += powerRatio*amp1*amp1*lorentzian(2*(x-x1), 2*fwhm1)
    intensity += powerRatio*amp2*amp2*lorentzian(2*(x-x2), 2*fwhm2)
    #intensity += powerRatio*amp3*amp3*lorentzian(2*(x-x3), 2*fwhm3)

    #crossover peaks
    intensity += 2*powerRatio*amp1*amp2*maxwellBoltzmann((x2-x1)*LIGHT_WAVELENGTH/2,T,MASS_RU_87)*lorentzian(2*x-x1-x2, fwhm1+fwhm2)
    #intensity += 2*powerRatio*amp2*amp3*maxwellBoltzmann((x3-x2)*LIGHT_WAVELENGTH/2,T,MASS_RU_87)*lorentzian(2*x-x2-x3, fwhm2+fwhm3)
    #intensity += 2*powerRatio*amp3*amp1*maxwellBoltzmann((x1-x3)*LIGHT_WAVELENGTH/2,T,MASS_RU_87)*lorentzian(2*x-x3-x1, fwhm3+fwhm1)

    intensity *= (bg0+x*bg1)
    return intensity

def saturationFitApprox2(x, x1 ,fwhm1, x2,fwhm2,x3,fwhm3,amp1, amp2,amp3, amp11, amp22, amp33, amp12, amp13, amp23,T, bg0, bg1):
    intensity = bg0+x*bg1

    #Doppler
    intensity -= amp1*dopplerVerbreiterung(x-x1, T, MASS_RU_87, 1.0/LIGHT_WAVELENGTH)
    intensity -= amp2*dopplerVerbreiterung(x-x2, T, MASS_RU_87, 1.0/LIGHT_WAVELENGTH)
    intensity -= amp3*dopplerVerbreiterung(x-x3, T, MASS_RU_87, 1.0/LIGHT_WAVELENGTH)

    #Lamb Peaks
    intensity += amp11*lorentzian(2*(x-x1), 2*fwhm1)
    intensity += amp22*lorentzian(2*(x-x2), 2*fwhm2)
    intensity += amp33*lorentzian(2*(x-x3), 2*fwhm3)

    #Crossover peaks
    intensity += amp12*lorentzian(2*x-x1-x2, fwhm1+fwhm2)
    intensity += amp13*lorentzian(2*x-x1-x3, fwhm1+fwhm3)
    intensity += amp23*lorentzian(2*x-x2-x3, fwhm2+fwhm3)

    return intensity

def saturationEmpiricParabolic(x, xMid, minVoltage,secondOrder, x1, fwhm1, x2, fwhm2, x3, fwhm3, amp11, amp22, amp33, amp12, amp13, amp23):
    intensity = minVoltage+secondOrder*(x-xMid)**2
    intensity += amp11*lorentzian(2*(x-x1), 2*fwhm1)
    intensity += amp22*lorentzian(2*(x-x2), 2*fwhm2)
    intensity += amp33*lorentzian(2*(x-x3), 2*fwhm3)

    intensity += amp12*lorentzian(2*x-x1-x2, fwhm1+fwhm2)
    intensity += amp13*lorentzian(2*x-x1-x3, fwhm1+fwhm3)
    intensity += amp23*lorentzian(2*x-x2-x3, fwhm2+fwhm3)

    return intensity

def saturationEmpiricGaussian(x, xMid, minVoltage,T, x1, fwhm1, x2, fwhm2, x3, fwhm3, amp11, amp22, amp33, amp12, amp13, amp23, bg0, bg1):
    intensity = bg0+bg1*x+(minVoltage-bg0-bg1*xMid)*dopplerVerbreiterung(x-xMid, T, MASS_RU_85, SPEED_OF_LIGHT/LIGHT_WAVELENGTH*1e-7)

    intensity += amp11*lorentzian(2*(x-x1), 2*fwhm1)
    intensity += amp22*lorentzian(2*(x-x2), 2*fwhm2)
    intensity += amp33*lorentzian(2*(x-x3), 2*fwhm3)

    intensity += amp12*lorentzian(2*x-x1-x2, fwhm1+fwhm2)
    intensity += amp13*lorentzian(2*x-x1-x3, fwhm1+fwhm3)
    intensity += amp23*lorentzian(2*x-x2-x3, fwhm2+fwhm3)

    return intensity

def lorentzian(x, fwhm):
    return fwhm**2/((x**2+fwhm**2))

def maxwellBoltzmann(v, T, m):
    val = -v**2*m/(2*BOLTZMANN_CONSTANT*T*SPEED_OF_LIGHT**2)
    #return sqrt(m/(2*pi*BOLTZMANN_CONSTANT*T*SPEED_OF_LIGHT**2))*np.exp(val)
    return np.exp(val)

def saturationIntegrand(v, x, x1, fwhm1, x2, fwhm2, T, m):
    return maxwellBoltzmann(v, T, m)*lorentzian(x-v/(2*SPEED_OF_LIGHT*LIGHT_WAVELENGTH)-x1, fwhm1)*lorentzian(x-v/(2*SPEED_OF_LIGHT*LIGHT_WAVELENGTH)-x2, fwhm2)

def absorbtionIntegrand(v, x, x1, fwhm1, T, m):
    return maxwellBoltzmann(v, T, m)*lorentzian(x-v/(2*SPEED_OF_LIGHT*LIGHT_WAVELENGTH)-x1, fwhm1)

def main():
    #fitAbsorbtionSpectrum("./absorbtion/F0007CH1.CSV", "./absorbtion/F0007CH2.CSV")
    fitSaturationSpectrum("./saturation/peak4/F0013CH1.CSV", "./saturation/peak4/F0013CH2.CSV")

main()
