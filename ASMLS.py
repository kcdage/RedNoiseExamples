###a recent version of astropy will be needed to run the LombScargle periodogram
### need to install DElightcurve simulation
##see README here: https://github.com/samconnolly/DELightcurveSimulation/blob/master/README.md
#https://github.com/samconnolly/DELightcurveSimulation.git

###run plotdata() to view the superorbital period variation with errorbars. default is to use .pickle file that is included in this folder
#to generate your own, run plotdata(gendata=True). WARNING: this takes ~2 hours because of all the steps involved. If I had more time, I would have parallelized this.

from DELCgen import *
import numpy as np
import matplotlib.pyplot as plt
import pickle
from astropy.stats import LombScargle
from scipy.optimize import curve_fit
import time
from sklearn import preprocessing
#import astropy.stats




def loaddata(verbose=True):

    '''loading data from text file and returning pertinent items as a dictionary'''
    #verbose toggle True will print out debugging statements.
    adata=np.genfromtxt('smcx1_831.dat.txt')

    
    if verbose:
        print "constructing ASCII data dictionary"

    ddict = {}
    ddict['TIME']  = adata[:,0]
    ddict['chisq'] = adata[:,4]
    ddict['RATE']  = adata[:,7]
    ddict['ERROR'] = adata[:,8]
    ddict['BKG'] = adata[:,12]

    return ddict


def examinelc(vrb=True):
    '''plots full lightcurve'''
    ##a look at the full light curve. vrb parameter toggles  print statements
    
    mydata=loaddata(verbose=vrb)

    vtime=mydata['TIME']
    vcounts=mydata['RATE']
    verr=mydata['ERROR']

    plt.errorbar(vtime, vcounts, yerr=verr, fmt='ko')
    plt.savefig('fulllc.png')
    plt.show()


    return
    

def sliverLS(nSliver=25, vrb=False, MaxChi=1.5, nBoots=1000, nMinPerStrip=50, DoPlot=False):

    '''I want to split my full dataset into T/nsliver stripes to better observe the superorbital period modulation'''

    ###nSliver splits the data into 25 stripes
    ###vrb is toggle for verbosity
    ##MaxChi is maximum chisq to consider in data
    ##nBoots is the number of bootstrap trials to run
    ###DoPlot is a toggle for plotting or not.

    anfang=time.time()
 

    
    mydata=loaddata(verbose=vrb)
    
    vtime=mydata['TIME']
    vrate=mydata['RATE']
    verr=mydata['ERROR']
    vchi=mydata['chisq']
    vbkg=mydata['BKG']

    ##filter data according to chisq and background but have at least 50 points per strip
    g = np.where((abs(vchi-1) < MaxChi) & (vbkg < 10))[0]

    tMin=np.min(vtime[g])
    tMax=np.max(vtime[g])

    delta=(tMax-tMin)/nSliver

    ##initalizing array to collect the slivered powers, as well as split times into.
    
    vArr=np.array([])
    vBegin=np.array([])
    vEnd=np.array([])


    ##now to add in parts for the bootstrap
    sim_err=np.zeros(nSliver)
    PeakPeriods=np.zeros(nSliver)
    bootstrappow=np.zeros(nSliver)
    PeaksStDev=np.zeros(nSliver)
    peakpowers=np.zeros(nSliver)
    jackknife=np.zeros(nSliver)

    burn=np.zeros(nSliver)
    ##the logic to sliver the dataset and plot the lightcurve and LS
    
    for j in range(nSliver):

        MinStrip=tMin+j*delta

        MaxStrip=MinStrip + delta


        gStrip=np.where( (vtime[g]>=MinStrip) & (vtime[g] <= MaxStrip) )[0]

        if np.size(gStrip) < nMinPerStrip:
            print("too few points in strip %f-%f, %i points" % (MinStrip, MaxStrip, np.size(gStrip)))
            continue

        else:
            print("working on time strip %f-%f, %i points" % (MinStrip, MaxStrip, np.size(gStrip)))

        ##okay so using autopower doesn't keep it with the same frequency space each time
        ##need to specify a frequency, then pass as LombScargle().power(freqs)

        vFreqs=np.linspace(0.001, 3, 10000)
        vPowers=LombScargle(vtime[g][gStrip], vrate[g][gStrip], verr[g][gStrip]).power(vFreqs)

        
        vPer=1.0/vFreqs
        #print vFreqs.size
    
        ##gathering up all LS powers
        
        #vArr.append(vPowers)
        #vBegin.append(MinStrip)
        #vEnd.append(MaxStrip)
        print('DBG--does it get this far?')
        if np.size(vArr)<1:
            vArr=np.copy(vPowers)
            vBegin=np.copy(MinStrip)
            vEnd=np.copy(MaxStrip)
            
        else:
            
            vArr=np.vstack((vArr, vPowers))
            vBegin=np.vstack((vBegin, MinStrip))
            vEnd=np.vstack((vEnd, MaxStrip))
        
            
            ##finding and gathering peak periods
        ThisPeak, PeakPow, jknife=FindPeakPeriod(vPer, vPowers, vrb=vrb, DoPlot=DoPlot)
        PeakPeriods[j]=ThisPeak
        peakpowers[j]=PeakPow

        jackknife[j]=jknife

        #print jknife, 'jackknife'

        

        #print ThisPeak, PeakPow
        #print peakpowers, 'blah blah blah'
        #print PeakPeriods, 'argelfraster'
        #print ThisPeak
                
        
        
        
        
        PeaksStDev[j]=BootstrapTrials(vtime[g][gStrip], vrate[g][gStrip], nBoots)

        outerr=np.zeros(5)
        p=0
        while p <5:
            print(j, ThisPeak)
            outerr[p]=tklcsim(vtime[g][gStrip], vrate[g][gStrip], verr[g][gStrip], j,  ThisPeak, DoPlot=True)
            p+=1
        sim_err[j]=np.mean(outerr)
     
           
           
      
        
        
        if DoPlot:
            ##trying to come up with a sensible file name (It didn't work)
            
            strmin=int(np.min(MinStrip))
            strmax=int(np.min(MaxStrip))
            
            ##plotting a sliver with errorbars
            plt.figure(1)
            plt.errorbar(vtime[g][gStrip], vrate[g][gStrip], yerr=verr[g][gStrip], fmt='bo')
            Fout="LC_%f-%f.png" %(strmin, strmax)
            plt.xlabel('Time, MJD')
            plt.ylabel('counts')
            plt.savefig(Fout)
            plt.show()

            
            plt.loglog(vPer, vPowers)
            plt.xlabel('Period (Days)')
            plt.ylabel('Powers')
            Fname='LS_%f-%f.png' %(strmin, strmax)
            plt.savefig(Fname)
            plt.show()


        ##returning and dumping as a dictionary so I don't have to run this over and over.
    #print peakpowers.shape, 'blah'
    SliverDict={'Powers':vArr, 'Periods':1.0/vFreqs,\
                        'TStart':vBegin, 'TEnd':vEnd,'boots':PeaksStDev, 'peaks':PeakPeriods, 'peakpowers':PeakPow, 'jknifeerr':jackknife,'simerr': sim_err}

    pickle.dump(SliverDict, open('StackedPowers.pickle', 'w'))

    ende=time.time()
    elap=(ende-anfang)/3600
    print(elap, "time elapsed for all parts of this code")
    return SliverDict

def GaussFunc(x, *P):

    '''returns a gaussian for a given peak'''
    ##takes in x range and three parameters and returns a gaussian
    

    return P[0]*np.exp(-(x-P[1])**2/(2.0*P[2]**2.0))
        
def FindPeakPeriod(vPeriods, vPowers, vrb=False, DoPlot=False, jack_knife=True):
    

    '''Need to return superorbital period value for the given sliver'''
    ##takes in periods and powers from the LS and finds the peak.
    #vrb is toggle for printing statements, DoPlot will toggle plots
    #returns peak period from a gaussian fit.

    #hardcoding in this range, since 70 is where the plateau starts
    #and below 3 is the orbital period
    gWhere=np.where((vPeriods<70.0) & (vPeriods>3.0))


    ##finding the largest peak
    vSig= np.argmax(vPowers[gWhere])

    if vrb:
        print(vSig)
        #print vPeriods[gWhere]

    ##First guess for peak value based on the above
    vGuess_Init=vPeriods[gWhere][vSig]

    if vrb:
        print(vGuess_Init, "guess peak period")

    
    vGuess=[vPowers[gWhere][vSig], vGuess_Init, 5.0]
    
    if vrb:
        print(vGuess)
        
    ##since I'm really only interested in the superorbital period,
    #restrict range to periods close to it.
    #axed this since curve_fit did not like it  and it did a pretty good job of fitting
    #gSup=np.where(np.abs(vPeriods-vGuess_Init)<2.0)[0]

    ##using scipy to fit to a gaussian
    
    coeff, fit_cov=curve_fit(GaussFunc, vPeriods, vPowers, vGuess, maxfev=2000000)

    
    ##evaluating the function
    MyGauss=GaussFunc(vPeriods, coeff[0], coeff[1], coeff[2])
    #print coeff[1]
    maxper_ind=np.argmax(MyGauss)

    maxper=vPeriods[maxper_ind]
    maxpow=MyGauss[maxper_ind]

    #print maxper, maxpow
    jknife_est=jackknife(vPeriods, GaussFunc, coeff[0], coeff[1], coeff[2])

    #print jknife_est
    
    if DoPlot:
        plt.figure(13)
        plt.clf()
        plt.semilogx(vPeriods, vPowers, 'r-.')
        plt.plot(maxper, maxpow, 'ko', zorder=10, markersize=5)
        plt.semilogx(vPeriods, MyGauss, 'b--')
        plt.errorbar(maxper, maxpow, xerr=jknife_est[maxper_ind])
        plt.ylim(0, 0.5)
        plt.show()

    if jack_knife==True:
        return maxper, maxpow, jknife_est[maxper_ind]
    else:
        return maxper, maxpow
    
def BootstrapTrials(vTimes, vRates, nTrials):

    '''bootstrap montecarlo, finding stDev of  peak period for LS of random subsets per trial'''

    ##takes in times and rates, as well as the number of trials over which to run the bootstrap
    ##returns standard deviation of peak values 
    
    ##Using this link to determine number of bootstrap trials to run
    ## https://stats.stackexchange.com/questions/86040/
    #rule-of-thumb-for-number-of-bootstrap-samples

    fPer=np.zeros(nTrials)
    burn=np.zeros(nTrials)

    ##Take a random sample of the given times/rates and LS them
    #find the peak period of the result
    #return standard deviation of results. This will be the error.
    
    for i in range(nTrials):

        #initalizing variables/arrays
        vPeak=0.0
        vSamplPow=np.zeros(vTimes.size)

        ##drawing random points for the bootstrap!

        #here are the indices I am using this time.
        sampl_ind=np.random.random_integers(0, (vTimes.size-1), vTimes.size)
        #print sampl_ind.size
        #print vTimes
        
        vSamplTime=vTimes[sampl_ind]
        vSamplRates=vRates[sampl_ind]

        vFreqs=np.linspace(0.001, 3, 2000)
        vPowers=LombScargle(vSamplTime, vSamplRates).power(vFreqs)

        vPer=1.0/vFreqs

        thisPeak, apow=FindPeakPeriod(vPer, vPowers, jack_knife=False)

        fPer[i]=thisPeak

    ##remove any nans
    #fPer[np.isnan(fPer)]=0
    ##remove any zero points
    ind_good=np.where(fPer >0)[0]

    stDev=np.std(fPer[ind_good])

    #print stDev
    return stDev
        
        
    
def plotdata(gendata=False):
    '''loads pickle dumped output of sliverLS and plots colormap of LS periodograms with errorbars on superorbital period. gendata toggle will run sliverLS for you if the file is not already in the directory'''
    if gendata==True:
        myPgrams=sliverLS()
    else:
        myPgrams=pickle.load(open('StackedPowers.pickle', 'r'))

    #pulling out the pieces I need for this plot
    vPer=myPgrams['Periods']
    zeit=0.5*(myPgrams['TStart']+myPgrams['TEnd'])
    LSPowers=myPgrams['Powers']
    berr=myPgrams['boots']
    pp=myPgrams['peaks']
    peakpowers=myPgrams['peakpowers']
    sim_err=myPgrams['simerr']
    j_err=myPgrams['jknifeerr']

    


    #print pp.shape, peakpowers.shape, berr.shape
    ##renormalizing powers to be at the maximum

    myMask=np.ma.masked_where(LSPowers<1000.0, LSPowers)
    g=np.where((vPer<70.0) & (vPer>3.0))
    
    #print myMask.shape

    for l in range(5,8):
	#print np.max(LSPowers[l][g])
        #for j in range(np.shape(LSPowers)[0]):
        plt.figure(1)
        plt.clf()
        plt.loglog(vPer, LSPowers[l])
        plt.errorbar(pp[l], np.max(LSPowers[l][g]), xerr=berr[l], color='m')
        plt.errorbar(pp[l], np.max(LSPowers[l][g]), xerr=sim_err[l], color='b')
        plt.errorbar(pp[l], np.max(LSPowers[l][g]), xerr=j_err[l], color='k')
        plt.xlim(0, 500)
        outname='powers_with_errors_%i.png'  %l
        plt.savefig(outname)
    #return
    
    #maxind=(np.shape(LSPowers))[0]
    normval=np.copy(LSPowers)
    #for k in range(maxind):
        #maxPow=np.max(LSPowers[k])
        #print zeit[k], maxPow
        #normval[k]/=maxPow
    logper=np.log10(vPer)
    
    
    ##specify extent of graph
    myExtent=[np.max(logper), np.min(logper), np.max(zeit), np.min(zeit)]

    #print zeit.shape, 'meh'
    
    #normval=preprocessing.normalize(LSPowers, norm='l2')
    ##setting up code for contour plot

    ##good contour examples here: https://matplotlib.org/examples/pylab_examples/contour_demo.html
  


    #plt.imshow(LSPowers,extent=myExtent, cmap='jet')

    #plt.xlabel('log10 period (Days)')
    #plt.ylabel('MJD (Oct1995-Mar2012)')
    #reverse axis
    #plt.xlim(2.0, 1.4)

    x_err=berr/pp
    #print pp.shape
    vzeit= np.transpose(zeit)[0]
     
    #return
    #plt.errorbar(pp[0], zeit[0], xerr=x_err[0], color='m', marker='o', lw=2, zorder=15)
    #plt.plot(pp, zeit, color='m', alpha=0.75)
   


    plt.figure(19)
    plt.clf()
    plt.plot(pp, zeit, alpha=0.25)
    plt.errorbar(pp, zeit, xerr=berr, color='g', fmt='.')
    
    plt.errorbar(pp, zeit, xerr=j_err, ecolor='b', fmt='.')
    plt.errorbar(pp, zeit, xerr=sim_err, color='r', fmt='.')
    plt.xlim(30,70)
    plt.xlabel('period (days)')
    plt.ylabel('time, MJD')
    plt.show()

            

        



def tklcsim(time, flux, flux_err, inte, peakperiod,DoPlot=False, InjectPeak=True, nSliver=25, MaxChi=1.5, nMinPerStrip=50):


    '''simulates red noise light curve according to Timmer and Koenig 1995, using DELCgen. takes in time, rates, the error, simulates a new light curve due only to red noise. Takes a toggle to inject a fixed peak in the data to see how well that is detected. Any variance in detection from the standard deviation of the detected peak versus the injected signal is the errorbar'''


    ## adding a break in the combined powerlaws at the frequency the signal is at
    #this is to test whether or not red noise could cause this specific freq
    nu_break=1.0/peakperiod
    pos_flux=np.copy(flux)
    for i in range(flux.size):
        if pos_flux[i]<0:
            
            pos_flux[i]=0
    #print np.where(flux<0), 'argelfraster'
    
    
    #plt.plot(time, pos_flux)
    #plt.show()
    #print nu_break.size, 'nu'
    ##simulation parameters
    RedNoiseL,RandomSeed,aliasTbin, tbin = 100,12,1,1000000

    ##make a lightcurve out of the time and fluxs
    
    mylc=Lightcurve(time, pos_flux, tbin, flux_err)

    #print time.size, flux.size, flux_err.size
    ##broken powerlaw parameters
    A=3.7
    a_low=1.1
    a_high=2.2
    c =0
    PSDpar=[A, nu_break, a_low, a_high, c]
    
    tklc=Simulate_TK_Lightcurve(BendingPL, (A, nu_break, a_low, a_high, c), lightcurve=mylc)
    
    if InjectPeak:
        #peakperiod=63.0
        peaksignal=A*np.sin((2.0*np.pi*tklc.time)/peakperiod)
        # Defaults for SMC X-1: OrbAmp = 2.5, OrbPeak = 3.89, OrbPhase=0.
        OrbitalSine = 2.5*np.sin(2.0*np.pi*(tklc.time)/3.89)
        tklc.flux=tklc.flux+peaksignal+OrbitalSine

    if DoPlot:
        plt.figure(5)
        plt.clf()
        plt.title('Light Curves')

        #plt.errorbar(mylc.time, flux, mylc.errors, fmt='ko', ecolor='k', alpha=0.5)
        plt.errorbar(tklc.time, tklc.flux, mylc.errors, fmt='mo',  ecolor='m')

        plt.xlabel('time(days)')
        plt.ylabel('rate (counts/s)')

      
        fname='real_and_simulated_%i.png' %inte 
        #plt.show()
        plt.savefig(fname)

    #print tklc.flux, 'argh'
    ##LombScargle real and simulated data for comparison
    Freq=np.linspace(0.001, 3, 10000)
    Data_LS=LombScargle(time, flux, flux_err).power(Freq)
    Sim_LS=LombScargle(tklc.time, tklc.flux, mylc.errors).power(Freq)
    
    if DoPlot:
        plt.figure(2)
        plt.clf()
        plt.title('LS Periodogram')
        plt.loglog(1.0/Freq, Data_LS,'k-.', label='Data', alpha=0.25)
        plt.loglog(1.0/Freq, Sim_LS, 'm--', label='Simulation')

        plt.xlabel('period (days)')
        plt.ylabel('power')
        f_out='LS_%i.png' %inte
        plt.savefig(f_out)


  



    if InjectPeak:

        ThisPeak, PeakPow, jknife=FindPeakPeriod(1.0/Freq, Sim_LS)

        diff=ThisPeak-peakperiod
        #print ThisPeak
        #print peakperiod
        #print diff
        return diff
 

        

      

    ###using inspiration from here:
    #https://people.duke.edu/~ccc14/sta-663/ResamplingAndMonteCarloSimulations.html

def jackknife(x, func, par1, par2, par3):
    """Jackknife estimate of the estimator func, takes in dataset, function and parameters, returns errors"""
    n = len(x)
    idx = np.arange(n)
    return np.sum(func(x[idx!=i], par1, par2, par3) for i in range(n))/float(n)

    

