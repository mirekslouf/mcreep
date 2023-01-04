'''
mcreep.fit
-----------
Fit various models to creep curves + calculate statistics of fitting.  
'''

import sys
import numpy as np
from scipy import optimize

def fit(MODEL, data, t_fstart, t_fend):
    '''
    Fit {MODEL} to {data}.

    Parameters
    ----------
    MODEL : mcreep.model.Model object
        This object contains all parameters (including fitting function)
        needed to fit model to data, calculate statistics, and show results.
    
    data : 2D numpy array
        Numpy array containing two columns with creep data
        (time and deformation).
        We note that this procedure is universal and fits data from
        both indentation and tensile creep experiments
        (see code below and comments inside it).
    
    t_fstart,t_fend : float,float
        The creep data are fitted to model function in interval
        [t_fstart; t_fend].
    
    Returns
    -------
    par,cov : list,list
        List of regression parameters {par}
        and their covariances {cov}
        for given fitting function
        (output from function scipy.optimize.curve_fit).
    
    '''
    # Cut data with t < t_fit
    data = data[:,(t_fstart<=data[0])&(data[0]<=t_fend)]
    # Split data into X,Y = t[s],h[um]
    X = data[0]
    Y = data[1]
    # Recalculate deformation data according to experiment type
    Y = recalculate_deformation(MODEL.EPAR, Y)
    # Fit data with given function
    # (Trick: we employ MODEL.iguess as initial guess of parameters
    # (...if iguess is not given, its default value None is used - Ok
    par,cov = optimize.curve_fit(MODEL.func,X,Y, p0=MODEL.iguess)
    # Return result
    return(par,cov)

def recalculate_deformation(EPAR,Y):
    '''
    Recalculation of deformation
    so that it could be fitted with model function in a correct way.

    Parameters
    ----------
    EPAR : Experiment object
        Experiment object contains information about creep experiment.
        In this function, we need parameters (EPAR.m,EPAR.K},
        which are used for h recalculation.
        
    Y : numpy array
        Deformation data for fitting.
        This data must be recalculated, as explained below.

    Returns
    -------
    Y : numpy array
        Y-data (= deformation) for fitting.
        The data are recalculated as necessary.

    Notes
    -----
    * Why the recalculation is needed?
        * The model function fits different types of deformation,
          depending on the type of creep experiment.
        * Model.EPAR.etype == 'Tensile'
          => deformation = strain = `epsilon(t)`
        * Model.EPAR.etype == 'Vickers' or 'Berkovich' or 'Spherical'
          => deformation = `f[h(t)] = [h(t)**m]/K`
        * The constants `m,K` are determined automatically within EPAR object
          (where EPAR = mcreep.const.Experiment).
    * What if we fitted just h(t) for indentatin creep experiments?
        * The fitting would work somehow,
          but the fitting/regression coefficients from idnentation creep
          would be incompatible with the coefficients from tensile creep.
        * Moreover, the fitting coefficients
          from different types of indentation experiments
          (for example Vickers vs. spherical) would not be comparable either.

    '''
    if EPAR.etype == 'Tensile':
        pass
    elif EPAR.etype in ('Vickers','Berkovich','Spherical'):
        Y = (Y**(EPAR.m))/(EPAR.K)
    else:
        sys.exit('Unknown model!')
    return(Y)

def coefficient_of_determination(MODEL, par, data):
    '''
    Calculate R2 = coefficient of determination ~ goodness of fit.
    
    Parameters
    ----------
    MODEL : mcreep.model.Model object
        This object contains all properties for calculation.
        
    par : list of floats
        Parameters of the fitting function
        (output from mcreep.fit.fit function).

    data : 2D numpy array
        XY data.
        Here: the creep data used for fitting (X = time, Y = deformation).
    
    Returns
    -------
    R2 : float
        The calculated coefficient of determination.

    Notes
    -----
    * R2 characterizes how well the data are predicted by fitting function.
    * It takes values:
      from -oo
      (extremely bad prediction, the mean of the data provides better fit)
      through 0
      (poor prediction, the mean of the data provides an equivalent fit)
      to +1
      (perfect prediction).
    * More info: <https://en.wikipedia.org/wiki/Coefficient_of_determination>
    * Here, R2 is calculated from the original {data} and {fitting function}.
    * The {fitting function} is saved in {MODEL.func}
      and its parameters are supplied in argument {par}.

    '''
    # Coefficient of determination = R2
    # R2 values: 1 = perfect, 0 = estimate ~ average(Y), negative = very bad. 
    # https://en.wikipedia.org/wiki/Coefficient_of_determination
    X = data[0]
    Y = data[1]
    # Recalculate deformation data according to experiment type
    Y = recalculate_deformation(MODEL.EPAR, Y)
    # Calculate R2 according.
    Yave = np.average(Y)
    Yfit = MODEL.func(X,*par)
    SSres = np.sum((Y-Yfit)**2)
    SStot = np.sum((Y-Yave)**2)
    R2 = 1 - SSres/SStot
    return(R2)
