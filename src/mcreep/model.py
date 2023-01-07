'''
mcreep.model
------------
Definition of Model for fitting the creep data + functions to run it.

The module defines Model class, whose object/instance contains:

* Properties for fitting model to data:
    1. property func: definition of fitting function
    2. property table_of_results: for saving the results of fitting
    3. and a few other properties, which can modify the fitting prodedure.
* Functions/methods for fitting model to data - the most important are:
    1. Model.run, which runs the whole calculation
    2. Model.describe, which briefly describes the model (for outputs, reports)
    3. Model.final_report, which creates final reports from ALL fitting results
* Three objects describing the input/output data:
    1.  mcreep.const.Experiment = description of experiment
    2.  mcreep.const.DataParameters = description of input data
    3.  mcreep.const.PlotParameters = description of the output plot
'''

import sys
import math
import pandas as pd
from pathlib import Path
import mcreep.io, mcreep.fit

class Model:
    '''
    Model for fitting creep data, containing namely:
        
    * fitting function (a function selected from mcreep.func)
    * methods to run the fitting (employing also mcreep.io, mcreep.fit)
    * table/dataframe with results of fitting; updated during calculations
    
    Parameters
    ----------
    EPAR : mcreep.const.Experiment object
        Description of experimental parameters (such as experiment type...).
    
    DPAR : mcreep.const.DataParameters object
        Description of data (such as which columns to read, rows to skip...)
    
    PPAR : mcreep.const.PlotParameters object
        Definition of plot parameters (xlabel, ylabel, global properties...).
    
    name : str
        Name of the fitting function in human-readable form.
    
    func : function object
        One of the fitting functions defined in module mcreep.func,
        such as `mcreep.func.power_law`, `mcreep.func.evp_s_d_1kv`, ...
    
    rtimes : list of floats; optional, default is None
        List of retardation times.
        If defined, the retardation times of EVP fitting functions
        are fixed at given values.
        
    iguess : list of floats; default is None
        Initial guess = list of starting values for fitting.
        The number of values must correspond to variables in given model.
        If defined, iguess parameter is passed to scipy.optimize.curve_fit
    
    print_covariances : bool; default is False
        If true, print covariance matrix.
        The diagonal elements of covariance matrix are `sigma**2(par1,par1)`
        => if they are close to zero, sigma of par1 is very low.
        The off-diagonal elements of covariance matrix = `sigma**2(par1,par2)`
        => if they are zero, par1 and par2 are independent on each other.
        Additional theory concerning the covariation matrix
        can be found elsewhere: https://stats.stackexchange.com/q/151018
        
    output_dir : str or path-like object; optional, default is '.'
        Directory for output files.
        If not given, we use current directory,
        i.e. the directory of the script that runs the program.

    Returns
    -------
    Model object :  
        * The Model.__init__ function
          returns the object with all properties and methods
        * The Model.run function
          fits the selected model function to experimental creep data
          and saves the results to {Model.table_of_results} property.
    
    Additional parameters
    ---------------------
    * table_of_results : pandas.DataFrame
        - Table of fitting results, which contains
          datafile name + values of fitted variables.
    * table_of_results_evp : pandas.DataFrame
        - Table of recalculated fitting results for EVP models, containing:
            * final compliances: `C0,Cv,C1,C2...`
            * final retardation times: `tau1,tau2...`
        - These values can be used to predict E(t):
            * `E(t) = f(C0,Cv,C1...,tau1...)`
    '''
    
    def __init__(self, EPAR, DPAR, PPAR, name, func, 
                 rtimes=None, iguess=None, print_covariances=False,
                 output_dir='.'):
        # Docstring for __init__ are given above in class description.
        # Reason: In this way, the parameters are visible in Spyder/Ctrl+I.

        # (1) Read basic properties.
        # (read during initialization
        self.EPAR = EPAR
        self.DPAR = DPAR
        self.PPAR = PPAR
        self.name = name
        self.func = func
        self.rtimes = rtimes
        self.iguess = iguess
        self.print_covariances = print_covariances
        self.output_dir = Path(output_dir)
        # (2) Modify/initialize additional properties
        # `(a) fname = name of the function; needed in future modifications
        self.fname = func.__name__
        # (b) results = pandas.Dataframes that keeps the results of fitting
        self._initialize_tables_of_results()
        # (3) Fix constants in EVP functions 
        # (a) fix basic constant (multiplicative constant in EVP models)
        if self.fname.startswith('evp'):
            self.func = self._evp_fix_basic_constants(EPAR)
        # (b) fix retardation times if requested => if rtimes != None
        if self.fname.startswith('evp') and rtimes != None:
            self.func = self._evp_fix_retardation_times(rtimes)
    
    def _initialize_tables_of_results(self):
        '''
        Initilize two tables,
        into which we are going to save fitting results.

        Returns
        -------
        None
            The result is correct initialization of two tables:
            table_of_results = for all models,
            table_of_results_evp = extra table for EVP models.

        '''
        # (1) Prepare columns = names of parameters for different models
        if self.fname == 'power_law':
            table1 = ['C','n']
            table2 = None
        elif self.fname == 'nutting_law':
            table1 = ['e0','C','n']
            table2 = None
        elif self.fname == 'evp_s_d_1kv':
            table1 = ['const','B0','Cv','D1','tau1']
            table2 = ['C0','Cv','C1','tau1']
        elif self.fname == 'evp_s_d_2kv':
            table1 = ['const','B0','Cv','D1','D2','tau1','tau2']
            table2 = ['C0','Cv','C1','C2','tau1','tau2']
        elif self.fname == 'evp_s_d_3kv':
            table1 = ['const','B0','Cv','D1','D2','D3','tau1','tau2','tau3']
            table2 = ['C0','Cv','C1','C2','C3','tau1','tau2','tau3']
        else:
            sys.exit('Unknown model!')
        # (2) Two additional columns for statistics
        # R2fit = coefficient of determination for fitted = [t_fstart; t_fend]
        # R2all = coefficient of determination for all data = [t_start; t_end]
        table1.extend(['R2fit','R2all']) 
        # (3) Create empty tables with prepared columns
        self.table_of_results = pd.DataFrame(columns=table1)
        self.table_of_results_evp = pd.DataFrame(columns=table2)            
            
    def _evp_fix_basic_constants(self, EPAR):
        '''
        Fix basic constants in EVP functions.

        Parameters
        ----------
        self : Model object
            Object with model fitting function.
        
        EPAR : Experiment object.
            Object with experimental data, including constants (F,K).
        
        Returns
        -------
        func_with_fixed_constants : function object
            The original function, but with fixed/constant parameters (F,K).
        
        Notes
        -----
        * SciPy algorithm for curve fitting `scipy.optimize.curve_fit`
          does not take into account keyword/fixed/constant arguments,
          i.e. it always fits all arguments regardless if they have value.
        * Therefore, the function has to be modified
          so that the keyword arguments are replaced by the values
          and removed from the fitting function definition.
        * https://stackoverflow.com/q/12208634
        
        '''
        # (1) Initialize
        # The following assignment is necessary.
        # If we use self.func directly, it causes mysterious errors.
        func = self.func
        # (2) Change functions = fix constants = eliminate constant arguments.
        # (here we just change the initial multiplicative constant
        # (the following assignement is just for convenience
        const = self.EPAR.const
        # Now we eliminate the const argument from the function definition...
        if self.fname == 'evp_s_d_1kv':
            func_with_fixed_constants = \
                lambda t,B0,Cv,D1,tau1 : \
                    func(t,const,B0,Cv,D1,tau1)
        elif self.fname == 'evp_s_d_2kv':
            func_with_fixed_constants = \
                lambda t,B0,Cv,D1,D2,tau1,tau2 : \
                    func(t,const,B0,Cv,D1,D2,tau1,tau2)
        elif self.fname == 'evp_s_d_3kv':
            func_with_fixed_constants = \
                lambda t,B0,Cv,D1,D2,D3,tau1,tau2,tau3 : \
                    func(t,const,B0,Cv,D1,D2,D3,tau1,tau2,tau3)
        else:
            # Non-EVP function => nothing to fix...
            func_with_fixed_constants = self.func
        # (3) Keep the name of original function
        # (this is not necessary here, fname is saved in self.fname
        # (however, it is "more correct" and may be necessary in non-OO code
        func_with_fixed_constants.__name__ = self.fname
        # (4) Return final function with fixed constants
        return(func_with_fixed_constants)
   
    def _evp_fix_retardation_times(self, rtimes):
        '''
        Fix retardation times in EVP functions
        if argument rtimes is not empty.

        Parameters
        ----------
        rtimes : list of floats
            Fixed values of retardation times of EVP models = tau1,tau2...
            The number of values must correspond to given model, of course.

        Returns
        -------
        None
        
        Notes
        -----
        * SciPy algorithm for curve fitting `scipy.optimize.curve_fit`
          does not take into account keyword/fixed/constant arguments,
          i.e. it always fits all arguments regardless if they have value.
        * Therefore, the function has to be modified
          so that the retardation times are replaced by constants
          and removed from the fitting function definition
          (if they are fixed by means of argument rtimes).
        * https://stackoverflow.com/q/12208634
        
        '''
        # (1) Initialize
        # The following assignment is necessary.
        # If we use self.func directly, it causes mysterious errors.
        func = self.func
        # (2) Change functions = fix constants = eliminate constant arguments.
        if self.fname == 'evp_s_d_1kv':
            func_with_fixed_rtimes = \
                lambda t,B0,Cv,D1 : \
                    func(t,B0,Cv,D1,rtimes[0])
        elif self.fname == 'evp_s_d_2kv':
            func_with_fixed_rtimes = \
                lambda t,B0,Cv,D1,D2 : \
                    func(t,B0,Cv,D1,D2,rtimes[0],rtimes[1])
        elif self.fname == 'evp_s_d_3kv':
            func_with_fixed_rtimes = \
                lambda t,B0,Cv,D1,D2,D3 : \
                    func(t,B0,Cv,D1,D2,D3,rtimes[0],rtimes[1],rtimes[2])
        else:
            # Non-EVP function => nothing to fix...
            func_with_fixed_rtimes = self.func
        # (3) Keep the name of original function
        # (this is not necessary here, fname is saved in self.fname
        # (however, it is "more correct" and may be necessary in non-OO code
        func_with_fixed_rtimes.__name__ = self.fname
        # (4) Return final function with fixed constants
        return(func_with_fixed_rtimes)
    
    def describe(self, fh=sys.stdout):
        '''
        Print brief description of the fitting model and its outputs.

        Parameters
        ----------
        fh : filehandle; optional, default is sys.stdout
            By default, the model description is printed on stdout.
            If {fh} is given, the output is redirected to {fh} filehandle.
            Reason: the description can be printed both to stdout
            and to text file.

        Returns
        -------
        None
            The output is the text printed to stdout or text file.
        '''
        # (0) If {fh} is given redirect standard output to filehandle.
        # (this is a dirty trick, but it simplifies/shortens code a bit
        if fh != None:
            original_standard_output = sys.stdout
            sys.stdout = fh
        # (1) Print description
        # (either to sys.stdout or to sys.stdout = fh
        # ...brief info about experiment/measurement type.
        if self.EPAR.etype == 'Tensile':
            print('Tensile creep experiment.')
        elif self.EPAR.etype == 'Vickers':
            print('Indentation creep with Vickers tip.')
        elif self.EPAR.etype == 'Berkovich':
            print('Indentation creep with Berkovich tip.')
        elif self.EPAR.etype == 'Spherical':
            print('Indentation creep with Spherical tip.')
        # ...information about the model and results.
        if self.fname == 'power_law':
            print('Model function: Power Law => deformation(t) = C * t**n')
            print('Units are relative; n = creep constant ~ creep rate.')
        elif self.fname == 'nutting_law':
            print("Model function: Nutting's Law => def(t) = e0 + C * t**n")
            print("Units are relative; n = creep constant ~ creep rate.")
        elif self.fname.startswith('evp'):
            print('Model function: EVP with S,D and KV components.')
            print('Compliances B,C,D in [GPa], retardation times tau in [s].')
        # ...final end-of-line
        print()
        # (2) If {fh} was used, restore original standard output...
        if fh != None:
            sys.stdout = original_standard_output
    
    def run(self, datafile, t_start, t_hold, t_fstart=None, t_fend=None):
        '''
        Run the model = fit model functin to experimental data.
        
        * This includes reading data + fitting + plotting + saving results.
        * The method employs other functions defined mcreep.io + mcreep.fit.

        Parameters
        ----------
        datafile : str or path-like object
            Full name of the datafile containing creep data.
        
        t_start : float
            The creep data are read for the interval
            [t_start; t_start+t_hold];
            t_start is the initial time for the data reading.
                For indentation experiments, t_start ~ when Fmax is reached.
                For tensile experiments, t_start ~ the first detected time.
                Alternatively, t_start can be a bit higher
                in order to ignore the initial period.
        
        t_hold : float
            The creep data are read for the interval
            [t_start; t_start+t_hold];
            t_hold is used to calculate the final time for the data reading.
                For both indentation and tensile experiments,
                it is the time for which the Fmax is held.
                It is reasonable to insert a slightly lower time
                (to be safe + to consider possible increase in t_start).
        
        t_fstart : float; optional
            The model function is fitted to creep data in the interval
            [t_fstart; t_fend].
            If t_fstart is not given, it is set equal to t_start;
            i.e. the function is fitted to all data that were read.
            
        t_fend : float; optional
            The model function is fitted to creep data in the interval
            [t_fstart; t_fend].
            If t_fstart is not given, it is set equal to (t_start + t_hold);
            i.e. the function is fitted to all data that were read.

        Returns
        -------
        None
            Nevertheless, this is the key method of {Model} object which...
            
            * reads creep data
            * fits model function to creep data
            * saves results to self.table_of_results
              and self.table_of_results_evp
            * Note: the printing and saving data to file
              is performed after running the model (self.run)
              by means of another method (self.final_report).
        '''
        # (0) Set t_fstart,t_fend to defaults,
        # if they were not given as arguments.
        if t_fstart == None: t_fstart = t_start
        if t_fend == None: t_fend = t_start+t_hold
        # (1) Read datafile = experimental data.
        data = self.read_datafile(datafile, t_start, t_hold)
        # (2) Fit datafile/experimental data with model function.
        par,cov = self.fit_function_to_data(data, t_fstart, t_fend)
        # (3) Print the result of fitting.
        # (Note: we convert datafile to pure/string name without path
        # (Reason: datafile might have been gi ven as pathlib object...
        # (...and the pathlib object cannot be printed easily as string
        datafile = Path(datafile).name
        self.print_fitting_result(datafile, par, cov)
        # (3) Plot the result of fitting.
        self.plot_fitting_result(datafile, data, par)
        # (4) Calculate statistics
        R2fit, R2all = self.calculate_statistics(par, data, t_fstart, t_fend)
        # (5) Save the result of fitting to MODEL object for final processing.
        self.save_fitting_result(datafile, par, R2fit, R2all, t_start)

    def read_datafile(self, datafile, t_start, t_hold):
        '''
        Read datafile with creep data
        (just a wrapper for function mcreep.io.read_datafile).
        '''
        data = mcreep.io.read_datafile(self, datafile, t_start, t_hold)
        return(data)
    
    def fit_function_to_data(self, data, t_fstart, t_fend):
        '''
        Fit model function to creep data
        (just a wrapper for function mcreep.fit.fit).
        '''
        par,cov = mcreep.fit.fit(self, data, t_fstart, t_fend)
        return(par, cov)
    
    def print_fitting_result(self, datafile, par, cov):
        '''
        Print the results of creep data fitting with given model.
        (just a wrapper for function mcreep.io.print_fitting_result).
        '''
        mcreep.io.print_fitting_result(self, datafile, par, cov)
    
    def plot_fitting_result(self, datafile, data, par):
        '''
        Save the results of creep data fitting with given model.
        (just a wrapper for function mcreep.io.plot_fitting_result).
        '''
        mcreep.io.plot_fitting_result(self, datafile, data, par)
    
    def calculate_statistics(self, par, data, t_fstart, t_fend):
        '''
        Calculate statistics for model function fitted to creep data.
        It calls mcreep.fit.coefficient_of_determination twice:
        
        * once for all data (= data in interval [t_start, t_start+t_hold])
        * once for fitted data (= data in interval [t_fstart, t_fend])
        
        Parameters
        ----------
        par : list of floats
            Result of fitting procedure;
            usually the result of procedure mcreep.fit.fit.
        
        data : 2D-numpy array
            Crep data read from input datafile;
            usually the result of procedure mcreep.io.read_datafile.
            
        t_fstart, t_fend : float, float
            Times determining the fitting interval = [t_fstart; t_fend].
        
        Returns
        -------
        R2fit, R2all : float, float
            Coefficients of determination for data in fitting interval
            [t_fstart, t_fend] and for the whole dataset in interval
            [t_start, t_start+t_hold].
        
        '''
        data_fit = data[:,(t_fstart<=data[0])&(data[0]<=t_fend)]
        R2fit = mcreep.fit.coefficient_of_determination(self, par, data_fit)
        R2all = mcreep.fit.coefficient_of_determination(self, par, data)
        return(R2fit, R2all)
    
    def save_fitting_result(self, datafile, par, R2fit, R2all, t_start):
        '''
        Save the results of fitting to Model object.
        The fitting results are saved in two object properties:
        
        * self.table_of_results = the results of fitting
        * self.table_of_results_evp = recalculated results for EVP models

        The saved results in Model object can be printed and saved using
        another property mcreep.model.Model.final_report.
        
        Parameters
        ----------
        
        datafile : str
            Datafile containing creep data.
            In this function it is used just for output;
            the (shortened) datafile name denotes the processed dataset.
        
        par : list of floats
            Result of fitting procedure;
            usually the result of procedure mcreep.fit.fit.
        
        R2fit, R2all : float, float
            Coefficients of determination for fitting interval and all data;
            output of the function mcreep.model.Model.calculate_statistics.
        
        t_start : float
            Starting time of fitting.
            In this procedure it is used for calculation of
            RCF = Ramp Correction Factors of EVP models
            (see the notes inside the code).
        
        Returns
        -------
        None
            * Formally it returns self, but this is not necessary.
            * The result are fitting results saved in self = in Model object.

        '''
        # (0) Define function that calculates RCF = ramp correction factors
        # (RCV = rho = constants employed in EVP recalculations below
        # {Mencik 2011} = Polymer Testing 30 (2011) 101–109; Eq.(9) at p.103.
        # (a) tR = time of ramping = until the maximum F is reached
        # (tR is used in rho definition and EVP calculations below
        tR = t_start
        # (b) rho = ramp correction factor, calculated from tau and tR
        def rho(tau): return( (tau/tR) * (math.exp(tR/tau) - 1) )
        # (1) Models
        # (A) Empirical models
        # (Aa) PL = Power Law
        if self.fname == 'power_law':
            results1 = {'C':par[0], 'n':par[1]}
        # (Ab) NL = Nutting's Law
        elif self.fname == 'nutting_law':
            results1 = {'e0':par[0],'C':par[1], 'n':par[2]}
        # (B) EVP models
        elif self.fname.startswith('evp'):
            # Add parameters that are common to all EVP models
            const = self.EPAR.const
            B0,Cv = par[0:2]
            results1 = {'const':const, 'B0':B0, 'Cv':Cv}
            # Prepare parameter rtimes (it will be used many times)
            rtimes = self.rtimes
            # (Ba) EVP: S+D+1KV
            if self.fname == 'evp_s_d_1kv':
                # Get parameters from par and/or rtimes
                # Two specific tricks for EVP1 (unlike EVP2,EVP3)
                # a) we use tau1=rtimes[0] ..to get just one float, not array
                # b) we use D1 = par[-1] ....to get just one float, not array
                # * Note: for EVP2,EVP3 (see below) it works like:
                #   ... D1,D2 = par[-2:] ....list of floats (from array slice)
                #   ... tau1,tau2 = rtimes ..list of floats (from whole array)
                # * Generalization (important to avoid sneaky errors):
                #   To get single float we cannot use whole array/array slice
                #   because this returns array of floats!
                if rtimes: D1 = par[-1]; tau1=rtimes[0]
                else: D1,tau1 = par[-2:]
                # Save FITTED compliances to results1
                results1 = {**results1, 'D1':D1}
                # Calculate and save FINAL compliances to results2
                if self.EPAR.etype == 'Tensile':
                    # calculation for tensile experiments => straightforward
                    C0 = results1['B0'] - D1
                    results2 = {'C0':C0, 'Cv':Cv, 'C1':D1}
                else:
                    # calculation for indentation experiments => uses RCF=rho
                    C1 = results1['D1']/rho(tau1)
                    C0 = results1['B0'] - results1['Cv']*tR/2 - C1
                    results2 = {'C0':C0, 'Cv':Cv, 'C1':C1}
                # Add retardation times to both results1 and results2
                results1 = {**results1, 'tau1':tau1}
                results2 = {**results2, 'tau1':tau1}
            # (Bb) EVP: S+D+2KV
            elif self.fname == 'evp_s_d_2kv':
                # Get parameters from par and/or rtimes
                if rtimes: D1,D2 = par[-2:]; tau1,tau2 = rtimes
                else: D1,D2,tau1,tau2 = par[-4:]
                # Save FITTED compliances to results1
                results1 = {**results1, 'D1':D1, 'D2':D2}
                # Calculate and save FINAL compliances to results2
                if self.EPAR.etype == 'Tensile':
                    # calculation for tensile experiments => straightforward
                    C0 = results1['B0'] - D1 - D2
                    results2 = {'C0':C0, 'Cv':Cv, 'C1':D1, 'C2':D2}
                else:
                    # calculation for indentation experiments => uses RCF=rho
                    C1,C2 = D1/rho(tau1),D2/rho(tau2)
                    C0 = results1['B0'] - results1['Cv']*tR/2 - C1 - C2
                    results2 = {'C0':C0, 'Cv':Cv, 'C1':C1, 'C2':C2}
                # Add retardation times to both results1 and results2
                results1 = {**results1, 'tau1':tau1, 'tau2':tau2}
                results2 = {**results2, 'tau1':tau1, 'tau2':tau2}
            # (Bc) EVP: S+D+3KV
            elif self.fname == 'evp_s_d_3kv':
                # Get parameters from par and/or rtimes
                if rtimes: D1,D2,D3 = par[-3:]; tau1,tau2,tau3 = rtimes
                else: D1,D2,D3,tau1,tau2,tau3 = par[-6:]
                # Save FITTED compliances to results1
                results1 = {**results1, 'D1':D1, 'D2':D2, 'D3':D3}
                # Calculate and save FINAL compliances to results2
                if self.EPAR.etype == 'Tensile':
                    # calculation for tensile experiments => straightforward
                    C0 = results1['B0'] - D1 - D2 - D3
                    results2 = {'C0':C0, 'Cv':Cv, 'C1':D1, 'C2':D2, 'C3':D3}
                else:
                    # calculation for indentation experiments => uses RCF=rho
                    C1,C2,C3 = D1/rho(tau1),D2/rho(tau2),D3/rho(tau3)
                    C0 = results1['B0'] - results1['Cv']*tR/2 - C1 - C2 - C3 
                    results2 = {'C0':C0, 'Cv':Cv, 'C1':C1, 'C2':C2, 'C3':C3}
                # Add retardation times to both results1 and results2
                results1 = {**results1, 'tau1':tau1, 'tau2':tau2, 'tau3':tau3}
                results2 = {**results2, 'tau1':tau1, 'tau2':tau2, 'tau3':tau3}
            # (Bd) Some other/unknown EVP model...
            else: sys.exit('Unknown EVP model!')
        # (D) Completely unknown model...
        else: sys.exit('Unknown model!')
        # (2) Save {results1} to self.table_of_results
        # (a) Add statistics
        results1 = {**results1, 'R2fit':R2fit, 'R2all':R2all}
        # (b) Add/save parameters to  self.table_of_results
        # onvert datafile to pure/string name without path
        # (Reason: datafile will be used as a (string) index of the DataFrame
        datafile = Path(datafile).name
        # (c) Add the new row (created above) to the DataFrame
        # (df.loc[new_index] => new row in df; here [new_index] = [datafile]
        # (other methods of adding rows to df - more difficult, bad indexing
        self.table_of_results.loc[datafile] = results1
        # (3) Save {results2} to self.table_of_results_evp (EVP models only)
        if self.fname.startswith('evp'):
            self.table_of_results_evp.loc[datafile] = results2
        # (4) Return is not necessary (data have been stored in self directly)
        return(self)
    
    def final_report(self, output_file=None):
        '''
        Print and save final report = summarize the results of fitting.
        The report is calls the following methods
        of Model object:
        
        * Model.describe = basic text description of the model.
        * Model.print_table_of_results = print results of fitting.
        * Model.print_table_of_results_evp = print addidional results
          of fitting for EVP models.
          
        The results are both printed and saved in the text file.

        Parameters
        ----------
        output_file : str, optional
            Name of the output text file, into which the results are saved.
            If the parameter is given

        Returns
        -------
        None
            Output is the printed report
            and saved file with the results of fitting.

        '''
        # (0) Get name of output file and open it
        # ...if output_file was not given, create default = input_file.txt
        if output_file == None: output_file = sys.argv[0] + '.txt'
        fh = open(output_file, 'w')
        # (1) Save basic description of the model to file
        self.describe(fh)
        # (2) Print + save basic table of results
        # (for all models, this table contains fitting + statistics
        self.print_table_of_results()
        self.save_table_of_results(fh)
        # (3) Print + save additional table of results for EVP models
        # (only for EVP: recalculation of parameters to predict J(t) = f(t)
        if self.fname.startswith('evp'):
            self.print_table_of_results_evp()
            self.save_table_of_results_evp(fh)
        # (4) End of report, close output file
        fh.close()
        
    def print_table_of_results(self):
        '''
        Print table of results as text
        (this is done by a small trick employing pandas.DataFrame).
        '''
        print('\nFitting results and statistics:\n')
        s = self.table_of_results.to_string(float_format = '%.4f')
        print(s)
        
    def print_table_of_results_evp(self):
        '''
        Print table of additional result of EVP models as text
        (this is done by a small trick employing pandas.DataFrame).
        '''
        print('\nFinal compliances and retardation times of EVP model:\n')
        s = self.table_of_results_evp.to_string(float_format = '%.4f')
        print(s)
        
    def save_table_of_results(self, fh):
        '''
        Save table of results as text
        (this is done by a small trick employing pandas.DataFrame).
        '''
        s = 'Fitting results & statistics:\n\n'
        s = s + self.table_of_results.to_string(float_format = '%.6f')
        fh.write(s)
        
    def save_table_of_results_evp(self, fh):
        '''
        Save table of additional results of EVP models as text
        (this is done by a small trick employing pandas.DataFrame).
        '''
        s = '\n\nFinal compliances and retardation times of EVP model:\n\n'
        s = s + self.table_of_results_evp.to_string(float_format = '%.6f')
        fh.write(s)
