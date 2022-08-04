This file code for making red noise simulations of data ('ASMLS.py'), as well as the writeup and figures for it. 

The DElightcurve package needs to be installed. See here: https://github.com/samconnolly/DELightcurveSimulation.git

ASMLS.py contains several modules. The first, "loaddata" loads the data from the text file and returns time, chisq, rate, error and background as a dictionary. 

The next module is "examinelc", which plots the full light curve.

The main part to this code is "sliverLS", which splits the dataset into 25 slivers, plots them, runs the LombScargle periodogram on them, and does error analysis on the peak periods. 

The following modules are called on from "sliverLS": 
"GaussFunc", which returns a gaussian mode for given parameters, "FindPeakPeriod", which fits the Gaussian to the given data and returns the maximum value. "BootstrapTrials", which implements the bootstrap Monte Carlo, tklcsim, which simulates red noise data in the form of the RXTE data, adds in a peak and returns the difference from the detected peak and the injected peak. Lastly, there is "jackknife", which returns a jackknife Monte Carlo error estimation. 

Finally, there is "plotdata", which takes in the output of sliverLS and plots it. 
# RedNoiseSimulations



