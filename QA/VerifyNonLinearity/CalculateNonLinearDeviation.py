from ROOT import *
import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
from scipy.interpolate import interp1d
from scipy import integrate


def gaus(x,a,x0,sigma):
        return a*exp(-(x-x0)**2/(2*sigma**2))

def pdfgaus(x,x0,sigma):
        return (1./np.sqrt(2*np.pi*sigma**2))*exp(-(x-x0)**2/(2*sigma**2))
    
def GetHistArray2D(hist):
    n     = [hist.GetNbinsY()+2,hist.GetNbinsX()+2]
    n2    = (n[0])*(n[1])
    d     = hist.GetArray()
    dX    = hist.GetXaxis().GetBinWidth(0)
    dY    = hist.GetYaxis().GetBinWidth(0)
    
    Xmin  = hist.GetXaxis().GetXmin()
    Xmax  = hist.GetXaxis().GetXmax()
    Xlow  = np.linspace(Xmin, Xmax+dX, n[1]-2)
    Xhigh = np.linspace(Xmin+dX, Xmax, n[1]-2)
    
    Ymin  = hist.GetYaxis().GetXmin()
    Ymax  = hist.GetYaxis().GetXmax()
    Ylow  = np.linspace(Ymin, Ymax+dY, n[0]-2)
    Yhigh = np.linspace(Ymin+dY, Ymax, n[0]-2)
    
    
    d.SetSize(n2)
    histA = np.array(d)
    histA = np.reshape(histA,n)[1:-1,1:-1]
    return histA,Xlow,Ylow
                                                                                    

def GetCalibrationCurve(hist,Yaxis):
    NY, NX = hist.shape
    Max = []
    Std = []
    X   = []
    for n in range(NX):
        try:
                y= hist[:,n]
                maxval = np.max(y)
                maxid  = np.argmax(y)
                X.append(Yaxis[n])
                popt,pcov = curve_fit(gaus,Yaxis,y,p0=[maxval,Yaxis[maxid],10], maxfev = 10000)
                Max.append(popt[1])
                Std.append(popt[2])
        except:
                Max.append(0)
                Std.append(0)
            
    return np.array(Max),np.array(Std)
        

## Load the various data
f          = TFile(sys.argv[1],"update")



for stage in range(0,5):
        ## Non-linear Function from my file
        calibGraph = f.Get("calWEPL/RangeVsEnergy_"+str(stage)) 
        calibY     = np.array(calibGraph.GetY())
        calibX     = np.array(calibGraph.GetX())
        NLF        = interp1d(calibX,calibY, fill_value = "extrapolate")

        ## Non-linear Function from the stage combined calibration
        #calibGraph = f.Get("calWEPL/NewRange_YX")
        #calibY     = np.array(calibGraph.GetY())
        #calibX     = np.array(calibGraph.GetX())
        #NLF        = interp1d(calibX,calibY, fill_value = "extrapolate")

        ## From previous calibration
        #calibList = np.loadtxt("WcalibHe.txt").reshape(5,340)
        #XRange    = np.arange(0,len(calibList[0]),1)
        #NLF        = interp1d(XRange,calibList[stage], fill_value = "extrapolate")

        
        ## Standard deviation
        calibhist  = f.Get("ProfileE")
        calibhist,Xaxis, Yaxis = GetHistArray2D(calibhist)
        Emax, StdEmax = GetCalibrationCurve(calibhist,Yaxis) ## Standard deviation of the Gaussian
        ids = np.where((Xaxis>-50)& (Xaxis<50))
        Std = np.mean(StdEmax[ids])*2## For now we fix it
        WETMean = []
        ## Define the region of integration
        Elist = np.arange(100,700,30)
        for E in Elist:
                Region = np.arange( E-3*Std, E+3*Std, 1) ## Region of interest in energy
                #Region = Region[Region<270] 
                PDF    = pdfgaus(Region, E,Std) ## Transformation with a Gaussian PDF
                WETMean.append(integrate.trapz(NLF(Region)*PDF,Region,dx=1))

                plt.plot(Region,PDF*1000)
                Transformed = NLF(Region)#*PDF        
                plt.plot(Region,Transformed,'k-') ## Calibration curve
                plt.plot([0,250],[NLF(E), NLF(E)],'r--')
                plt.plot([E, E], [0, 250],'r--')
                
                ## New Distribution
                WETPDF = Transformed*PDF/np.sum(Transformed*PDF)
                Mean   = np.average(Transformed, weights=PDF)
                MeanID = np.argmin(np.abs(Transformed-Mean))
                print Mean, NLF(E)
                plt.plot(WETPDF*1000,Transformed)
                plt.plot(WETPDF[MeanID]*1000,Transformed[MeanID],'ro')
        plt.title("Stage "+str(stage))
        plt.xlim([0,250])
        plt.ylim([NLF(Elist[-1]),NLF(Elist[0])])
        plt.xlabel("Energy")
        plt.ylabel("WET")
        plt.show()
