'''
mcreep.const
------------
Global constants for package mcreep.

* Global constants are saved in the form of objects.
* Reason: Object parameters can be [re]defined during initialization.
'''

from math import pi,sqrt,tan
import matplotlib.pyplot as plt

class Experiment:
    '''
    Experiment = object describing the creep experiment.
    
    Parameters
    ----------
    etype : str; which takes one o the following values
        'Tensile'   = (macroscale) tensile creep experiment.
        'Vickers'   = indentation creep with Vickers tip.
        'Berkovich' = indentation creep with Berkovich tip.
        'Spherical' = indentation creep with spherical tip.
    F : float; optional, default is None
        Loading force during indentatin experiments in [mN].
        Relevant (and obligatory) for etype='Vickers/Berkovich/Spherical'.
    R : float; optional, the default is None
        Radius of spherical tip in [um].
        Relevant (and obligatory) for etype='Spherical'.
    sigma : float; optional, default is None
        Stress during tensile experiments in [GPa].
        Relevant (and obligatory) for etype='Tensile'.
    
    Additional parameters
    ---------------------
    * m,K = additional constants that are used
      for correct recalculations of indentation creep data.
      These parameters are calculated automatically,
      based on the value of parameter `etype`.
    * const = yet additional constant that is used
      for correct fitting of EVP models
      to both indentation and tensile creep data with EVP models.
        - EVP in general:
          `def(t) = const * { A0 + Av*t + Sum[f(Ai,ti)] }`
        - EVP in tensile creep:
          `def(t) = epsilon(t)` and `const = sigma[GPa]`     
        - EVP in indentation creep:
          `def(t) = [h(t)**m]/K` and `const = F[mN]`
        - Note for indentation creep:
            - Mencik_2009: `def(t) = h(t)**m` and `const = F*K`
            - Here (see above): `def(t) = [h(t)**m]/K` and `const = F[mN]`
            - Reason for the difference:
                - more consistent treatment of PL x NL x EVP models
                - for all models we fit: `[h(t)**m]/K`
    
    Sample calculation of sigma
    ---------------------------
    * `sigma[GPa] = (load[kg]*SA*g[m/s2]) / (W[mm]*T[mm]/1e6) / 1e9`
    * `SA` = stress amplifier (optional lever amplifying the load)
    * `g` = net acceleration = gravity of Earth = 9.81[m/s2]
    * `W,T` = width and thickness of testing specimen in [mm]
    * `1e6` = recalculate W*T[mm2] -> [m2]
    * `1e9` = recalculate final result [Pa] -> [GPa]
    
    '''
    def __init__(self, etype, F=None, R=None, sigma=None):
        # Docstring for __init__ is given above in the class definition.
        # Reason: consistent help in Spyder and Pdoc.
        # -----        
        
        # (1) Initialize basic parameters
        self.etype = etype
        self.F = F
        self.R = R
        self.sigma = sigma
        
        # (2) Calculate additional parameters: K,m
        # (parameters K,m calculated according to {Mencik 2011}
        # {Mencik 2011} = Polymer Testing 30 (2011) 101–109; Eq.(9) at p.103.
        if self.etype == 'Tensile':
            self.K = 1
            self.m = 1
        elif self.etype in ['Vickers','Berkovich']:
            alpha = 70.3 * pi/180
            self.K = pi/(2*tan(alpha))
            self.m = 2
        elif self.etype == 'Spherical':
            self.K = 3/(4*sqrt(self.R))
            self.m = 3/2
        
        # (3) Calculate additional parameter: const
        # (const is employed in EVP models,
        # (...in which it is the 1st multiplicative constant
        # (...that differs for tensile and indentation experiments;
        # (this is my generalization, which enables that the package
        # (...can be used for both tensile and indentation experiments        
        if self.etype == 'Tensile':
            self.const = self.sigma
        elif self.etype in ['Vickers','Berkovich','Spherical']:
            # Mencik_2009:      `def(t) = h(t)**m` and `const = F*K`
            # Here (see above): `def(t) = [h(t)**m]/K` and `const = F[mN]`
            # Reason for the difference:
            # consistent treatment of PL x NL x EVP models.
            self.const = self.F


class DataParameters:
    '''
    DataParameters = object defining format of a creep datafile.
    
    Assumptions:
    * creep datafile is a TXT file containing data in columns;
    * one of the columns contains time, some other contains deformation.
    
    Parameters
    ----------
    usecols : list with two integer values; optional, the default is [0,1]
        This parameter is passed to numpy.loadtxt.
        * The 1st value of the list = column with times.
        * The 2nd value of the list = column with deformations.
    
    comments : string or sequence of strings; optional, the default is '#'
        This parameter is passed to numpy.loadtxt.
        The character(s) indicate comment/ignored lines in input datafile.
    
    skiprows : integer; optional, the default is 1
        This parameter is passed to numpy.loadtxt.
        The first {skiprows} lines are skipped.
    
    time_to_seconds : float; optional, the default is 1
        Multiplicative constant, which converts time values to seconds;
        this conversion is necessary for the following calculations.
    
    deformation_to_um : float; optional, the default is 1
        Multiplicative constant, which converts deformation to micrometers;
        this conversion is necessary for the following calculations. 
    '''
         
    def __init__(self, usecols=[0,1], comments='#', skiprows=0,
                 time_to_seconds=1, deformation_to_um=1):
        # Docstring for __init__ is given above in the class definition.
        # Reason: consistent help in Spyder and Pdoc.
        # -----        
        
        self.usecols = usecols
        self.comments = comments
        self.skiprows = skiprows
        self.time_to_seconds = time_to_seconds
        self.deformation_to_um = deformation_to_um
        
class PlotParameters:
    '''
    PlotParameters = object defining local+global parameters for plotting.
    
    Parameters
    ----------
    
    xlabel, ylabel : str, str
        Labels for X and Y axis.
    
    logscale : bool; optional, the default is False
        If logscale==True, both X and Y axes are in logarithmic scale.
    
    e_to_percent : bool; optional, the default is True
        Relevant only to tensile experiments.
        If true, the values of elongation are multiplied by 100,
        i.e. they are converted from epsilon[] to epsilon[%].
        
    rcParams : dict; optional, the default is empty dictionary {}
        The dictionary shoud be formatted for mathplotlib.pyplot.rcParams.
        The argmument is passed to matplotlib.pyplot.
        The initialization procedure creates some default rcParams.
        This argument can override this pre-defined parameters,
        i.e. the default is created anyway
        and then (possibly) supplemented by rcParams argument.
    
    showfigs : bool; optional, the default is True.
        If showfigs==True, the figures are shown + saved in files,
        which is default behavior, suitable for running the script in Spyder.
        If showfigs==False, the figures are just saved in files,
        which is an option, suitable for running the script from CLI.
        
    ax : matplotlib Axes object, the default is None
        * If ax == None, create and save results as a single plot,
          which is a typical usage.
        * If ax is defined, create the plot within given ax object,
          which can combined with `fig,ax = plt.subplots()`
          in order to create multile figures.
    '''
    
    def __init__(self, xlabel, ylabel,
                 logscale=False, e_to_percent=True,
                 rcParams={}, showfigs=True, ax=None):
        # Docstring for __init__ is given above in the class definition.
        # Reason: consistent help in Spyder and Pdoc.
        # -----        
        
        # (1) Initialize basic parameters
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.logscale = logscale
        self.e_to_percent = e_to_percent
        self.rcParams = rcParams
        self.showfigs = showfigs
        self.ax = ax
        
        # (2) Set global plot settings using rcParams
        PlotParameters.set_default_rcParams(rcParams)
        
    @classmethod
    def set_default_rcParams(cls, my_rcParams={}):
        '''
        A class method defining global plot parameters (plt.rcParams).
        
        Parameters
        ----------
        my_rcParams : dict
            The dictionary re-defines selected plt.rcParams keys.
            
            Example:
                
            >>> PlotParameters.set_default_rcParams({'figure.dpi':500})
            
        Returns
        -------
        None
            The function does not return anything,
            BUT it re-defines the global variable rcParams.

        Notes
        -----
        * This is a @classmethod (because it is used within the whole class)
          but it could be a @staticmethod as well (because it does not use
          cls variable in fact).
        * The method is employed in two ways:
            - Standard usage of MCREEP package: default rcParams are used
              (and possibly modified) in objects of PlotParameters class.
            - Special usage of MCREEP (more figures, multiplots): default
              rcParams are used when definining the axes of (multiple) figures.
        '''
        
        # (1) Set default rcParams
        # (Hardcoded, suitable default for standard plots
        plt.rcParams.update({
            'figure.figsize'     : (8/2.54,6/2.54),
            'figure.dpi'         : 500,
            'font.size'          : 7,
            'lines.linewidth'    : 0.8,
            'axes.linewidth'     : 0.6,
            'xtick.major.width'  : 0.6,
            'ytick.major.width'  : 0.6,
            'grid.linewidth'     : 0.6,
            'grid.linestyle'     : ':'})
        
        # (2) Update default with argument rcParams, if it was given
        # (User-defined in the main program, if necessary
        # (Useful namely for multiplots
        plt.rcParams.update(my_rcParams)
