'''
mcreep.io
-----------
Input/output functions for package mcreep.    
'''

import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def read_datafile(MODEL, datafile, t_start, t_hold):
    '''
    Read datafile containing creep data.
    
    Parameters
    ----------
    MODEL : MODEL object
        The object keeps experiment parameters in MODEL.EPAR property.
        We use the object to convert time and deformation to correct units.
    datafile : string or pathlib object
        Name of datafile = a text file containing
        at least two columns (time,deformation);
        see section {Expected format of the datafile}
        below for more information.
    t_start : float
        Start of step II = holding step = period at which F = Fmax.
    t_end : float
        End of step II = holding step = period at which F = Fmax.
        
    Returns
    -------
    2D numpy array.
    
    Expected format of the datafile
    -------------------------------
    * Two colums (time,deformation) separated by whitespace.
        - 1st column = t[s] = time in seconds
        - 2nd column = deformation[um] = deformation in micrometers
            - in tensile experiments: deformation = strain
            - in indentation experiments: deforamtion = penetration depth
    * The datafile is read by means of numpy.loadtxt function
      and some parameters can be adjusted by means of MODEL argument
      (which contains DataParameters object describing the input data):
        - see source code of this function
        - and mcreep.const.DataParameters object
        - for example, it is possible to specify column numbers
          if the datafile contains more columns than required
    '''
    # Read file to numpy array and try to catch possible errors/exceptions
    try:
        # read h-t file
        data = np.loadtxt(
            datafile, unpack=True,
            usecols = MODEL.DPAR.usecols,
            skiprows = MODEL.DPAR.skiprows,
            comments = MODEL.DPAR.comments)
        # take just section II of h-t curve
        # section I = loading  = (0..t_start] ...without 0 due to logarithms
        # section II = holding = [t_start..t_hold] ...time of maximal loading
        data = data[:,(data[0]>=t_start) & (data[0]<=(t_start+t_hold))]
        # recalculate time and deformation to [s] and [um], respectively
        data[0] = data[0] * MODEL.DPAR.time_to_seconds
        data[1] = data[1] * MODEL.DPAR.deformation_to_um
    except OSError as err:
        print('OSError:', err)
        sys.exit()
    except ValueError:
        print('ValueError: probably a wrong format of the file.')
        print('Expected: TXT file with two columns t[s], def[length_units]')
        print('...where def[length_units] is conveted to def[um]')
        print('...with a user-defined constant MODEL.DPAR.deformation_to_um')
        sys.exit()
    # Return final 2D numpy array
    # (1st column = t[s], 2nd col = def[um], range = [t_start;t_start+t_hold]
    return(data)

def print_fitting_result(MODEL, datafile, par, cov):
    '''
    Print ONE result of fitting
    (i.e. print fitting results for one/currently processed datafile).

    We note that the complete results of each fitting (for each datafile)
    are kept in mcreep.model.Model object and usually reported/printed
    together at the very end of the whole processing.
    
    Exception: Covariance matrixes of each fitting are just printed here,
    NOT kept in mcreep.model.Model object and NOT saved to file.
    Moreover, they are printed just on request (MODEL.print_covariances=True).
    They are needed just occasionally, and if really needed, the values
    can be copy+pasted from stdout and saved to file.

    Parameters
    ----------
    MODEL : mcreep.model.Model object
        This object contains all properties needed to print the result.
        
    datafile : str
        Name of the datafile;
        in this procedure, it is just printed to stdout
        in order to denote the source of the fitting results.
        
    par : list of floats
        List of optimized parameters of given fitting function;
        output from the procedure mcreep.fit.fit.
        
    cov : 2D-array of floats
        Covariance matrix for all fitted/regression parameters;
        supplementary output from the procedure mcreep.fit.fit.

    Returns
    -------
    None
        The output are results printed on stdout.
    '''
    # (1) Convert parameters array to tuple => suitable for formatted printing 
    par = tuple(par)
    # (2) Print datafile name (this is the same for all models)
    print(f'{datafile} ', end='')
    # (3) Print fitting/regression parameters...
    # * The parameters are printed in the format/order specific to given model.
    if MODEL.fname == 'power_law':
        print('[C,n]: %8.4f %8.4f' % par)
    elif MODEL.fname == 'nutting_law':
        print('[e0,C,n]: %8.4f %8.4f %8.4f' % par)
    elif MODEL.fname == 'evp_s_d_1kv':
        if MODEL.rtimes:
            print('[B0,Cv,D1]: %8.4f %8.4f %8.4f' % par)
        else:
            print('[B0,Cv,D1,tau1]: %8.4f %8.4f %8.4f %6.2f' % par) 
    elif MODEL.fname == 'evp_s_d_2kv':
        if MODEL.rtimes:
            print('[B0,Cv,D1,D2]: ', end='')
            print('%8.4f %8.4f %8.4f %8.4f' % par)
        else:
            print('[B0,Cv,D1,D2,tau1,tau2]: ', end='')
            print('%8.4f %8.4f %8.4f %8.4f %6.2f %6.2f' % par) 
    elif MODEL.fname == 'evp_s_d_3kv':
        if MODEL.rtimes:
            print('[B0,Cv,D1,D2,D3]: ',end='')
            print('%8.4f %8.4f %8.4f %8.4f %8.4f' % par)
        else:
            print('[B0,Cv,D1,D2,D3,tau1,tau2,tau3]: ', end='')
            print('%8.4f %8.4f %8.4f %8.4f %8.4f %6.2f %6.2f %6.2f' % par)
    else:
        print('Warning: unknown model during printing!')
        print(par)
    # (4) Print covariance matrix showing independence of parameters...
    # * The cov.matrix is printed on request (MODEL.print_covariances=True).
    # * The order of rows/columns of the cov.matrix the same like in item (3).
    #   => therefore, we do not re-print the parameter names like in item (3).
    # * Fitting/regression parameters are usually saved to file
    #   (by means of MODEL.final_report() method).
    # * Covariance matrixes is just printed to stdout here.
    #   (they can be copy+pasted from stdout and saved to file manually).
    if MODEL.print_covariances == True:
        print('\nCovariance matrix of all parameters after fitting:')
        for row in cov:
            for val in row:
                print(f'{val:10.6f}', end='')
            print()

def plot_fitting_result(MODEL, datafile, data, par):
    '''
    Plot the result of fitting (single plot OR axes defined within MODEL).
        
    Parameters
    ----------
    MODEL : mcreep.model.Model object
        This object contains all properties needed to print the result.
        
    datafile : str or path-like object
        Name of the datafile containing input data.
        * Important: We use this argument just for creating
          the name of the output graph => {datafile}.png.
        * If MODEL.PPAR.ax == None: the output goes to single PNG file,
          whose name is created as `datafile.PNG`.
        * If MODEL.PPAR.ax != None: the output goes to pre-prepared axes
          object, and datafile parameter is ignored.
    
    data : 2D numpy array
        Creep data read from the currently processed datafile.
        It is expectted 
        
    par : list of floats
        List of optimized parameters of given fitting function;
        output from the procedure mcreep.fit.fit.

    Returns
    -------
    None
        The output is either single plot OR axes in a multiplot.
        * If MODEL.PPAR.ax == None:
          a single output graph is shown (option) + saved (always).
          The name of the output graph is {datafile}.png.
        * If MODEL.PPAR.ax == axes_object
          a graphs is created within the pre-defined axes_object.
          This can be employed in creating user-defined multiplots.
    '''
    # (1) Read data
    X,Y = data
    # (2) Calculate fitted data
    # Trick: *(list(par)) = convert saved parameters to list and expand
    Yfit = MODEL.func(X, *(list(par)))
    # (3) For tensile experiments, convert elongation to % if required.
    if MODEL.EPAR.etype == 'Tensile' and MODEL.PPAR.e_to_percent == True:
        Y = Y*100
        Yfit = Yfit*100
    # (4) Recalculate fitting function
    # Reason: Tensile x Indentation experiments fit different deformations...
    # More details => see the description of the recalculating fucntion.
    Yfit = recalculate_fitted_data(MODEL.EPAR, Yfit)
    # (5) If required, convert data to logscale
    # AND change x/ylabels accordingly.
    if MODEL.PPAR.logscale == True:
        # (a) Convert data to logscale
        X = np.log10(X)
        Y = Y = np.log10(Y)
        Yfit = np.log10(Yfit)
        # (b) change x/ylabels
        # Trick: this function is called repeatedly
        # => check if the change of x/ylabes has not been already done!
        if not(MODEL.PPAR.xlabel.startswith('log')):
            MODEL.PPAR.xlabel = 'log('+MODEL.PPAR.xlabel+')'
        if not(MODEL.PPAR.ylabel.startswith('log')):
            MODEL.PPAR.ylabel = 'log('+MODEL.PPAR.ylabel+')'
    # (6) Create plot
    # (user-adjustable + global plot settings are saved in PPAR object
    # (MODEL.PPAR object is an instance of mcreep.const.PlotParameters class
    # (6a) No axes_object was given as argument => create+save single plot
    if MODEL.PPAR.ax == None:
        # (a) Create name of output graph => name of the plot to save.
        # (output graph = MODEL.output_dir/datafile.png
        output_filename = Path(datafile).name + '.png'
        output_graph = Path(MODEL.output_dir, output_filename)
        # (b) Create the plot
        plt.plot(X, Y, color='orange', label='Experiment')
        plt.plot(X, Yfit, 'k:', label=MODEL.name)
        plt.xlabel(MODEL.PPAR.xlabel)
        plt.ylabel(MODEL.PPAR.ylabel)
        plt.grid()
        plt.legend()
        plt.tight_layout()
        # (c) Save the plot
        plt.savefig(output_graph)
        # (d) Show and close the plot
        # (Default is showfigs=True => figs shown+closed - suitable for Spyder
        # (Option is to set showfigs==False => just close - suitable for CLI
        # (Technical notes:
        # ( * plt.show ...show (and close) plot in Spyder
        # ( * plt.close ..close plot explicitly
        # (   needed for multiple plots + showfigs==False
        # (   otherwise all plots would be drawn into just one figure
        if MODEL.PPAR.showfigs == True:
            plt.show()
        plt.close()
    # (6b) axes object was given as argument => create plot within the axes
    # (this is employed when creating multiple plots or multiplots
    else:
        ax = MODEL.PPAR.ax
        ax.plot(X, Y, color='orange', label='Experiment')
        ax.plot(X, Yfit, 'k:', label=MODEL.name)
        ax.set_xlabel(MODEL.PPAR.xlabel)
        ax.set_ylabel(MODEL.PPAR.ylabel)
        ax.grid()
        ax.legend()
        
def recalculate_fitted_data(EPAR, Y_orig):
    '''
    Recalculation of the fitted/calculated data before plotting.
   
    Parameters
    ----------
    EPAR : Experiment object
        Object with experimental parameters.
        In this function, we need (EPAR.m,EPAR.K) for Y-data recalculation.
    Y_orig : 1D numpy array.
        Original Y-data, to which f(t) was fitted.
        The data have to be recalculated, because...
        
        * The models could be fitted to tensile or indentation creep data.
            * Tensile creep data were not modified: deformation(t) = epsilon(t)
            * Indentation creep data were modified: deformation(t) = (h**m)/K)
        * The modification is necessary in order to achieve
          compatibility between tensile and indentation creep results.
            * More precisely, indentation creep data are recalculated
              so that the results
              from (specific) fitting of models to indentation creep data
              are compatible with the results
              from (standard) fitting of models to tensile creep data.
        * For INDENTATION CREEP, the modification means that...
            * Model was not fitted to [h(t)] but to [h(t)**m]/K
            * Fitted data calculated from model are [h(t)**m]/K
            * These data have to be calculated to [h(t)] before plotting:
                * `Y_orig = [h(t)**m]/K`
                * `h(t) = (Y_orig * K) ** (1/m)`

    Returns
    -------
    Y_recalc : 1D numpy array
        Recalculated Y-data for plotting.
        The recalculated data are:
            * epsilon(t) for tensile creep
            * h(t) for indentation creep
    
    '''
    if EPAR.etype == 'Tensile':
        Y_recalc = Y_orig
    else:
        Y_recalc = (Y_orig * EPAR.K) ** (1/EPAR.m)
    return(Y_recalc)