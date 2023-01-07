'''
mcreep.func
-----------
Definition of functions for fitting.

* The functions are defined in general form:`deformation(t) = f(t)`
* The package fits the functions to both tensile and indentation creep data.

Pre-defined functions available in current version:

* power_law = empirical Power Law model:   `deformation(t) = C * t**n`
* Nutting's law = empirical Nutting's law: `deformation(t) = e0 * C * t**n`
* evp_s_d_1kv = EVP model with Spring + Dashpot + 1 Kelvin/Voigt element
* evp_s_d_2kv = EVP model with Spring + Dashpot + 2 Kelvin/Voigt elements
* evp_s_d_3kv = EVP model with Spring + Dashpot + 3 Kelvin/Voigt elements

The modifications for tensile and indentation creep data:
    
* In both tensile and indentation creep, we have: `deformation(t) = f(t)`.
* In tensile creep: `deformation(t) = strain(t) = epsilon(t)`.
* In indentation creep: `deformation(t) = [h(t)**m]/K`.
    - the constants (m,K) are taken from publication {Mencik 2011}
    - this ensures the compatibility between tensile and indentation creep
    - more details can be found in our own publication {Slouf 2013}
'''

import numpy as np

def power_law(t, c, n):
    '''
    Function defining PL model (Power Law model).
        
    Parameters
    ----------
    t : float
        Time = X-variable for fitting.
    
    c, n : float, float
        Parameters of power law model.

    Returns
    -------
    Function
        Expression for fitting procedure.
    
    '''
    # Define fitting function - here - power law
    func = c * t**n
    # Return final function
    return(func)

def nutting_law(t, e0, c, n):
    '''
    Function defining NL model (Nutting's law/model).
        
    Parameters
    ----------
    t : float
        Time = X-variable for fitting.
    
    e0, c, n : float, float, float
        Parameters of Nutting's model.

    Returns
    -------
    Function
        Expression for fitting procedure.
    
    '''
    # Define fitting function - here - power law
    func = e0 + c * t**n
    # Return final function
    return(func)
    

def evp_s_d_1kv(t, const,B0,Cv, D1, tau1):
    '''
    Function defining EVP model with elements [S + D + 1*KV].

    Parameters
    ----------
    t : float
        Time = X-variable for fitting.

    const : float
        Multiplicative constant for EVP models.
        * Tensile experiments: const = applied stress = sigma[GPa]
            - `sigma` = constant, calculated+saved in mcreep.const.Experiment
        * Indentation experiments: const = `F*K`
            - `F[mN]` = loading force
            - `K` = constant, calculated+saved in mcreep.const.Experiment

    B0 : float
        Parameter corresponding to S-element of EVP model.

    Cv : float
        Parameter corresponding to D-element of EVP model.

    D1 : float
        Parameter of KV-element of EVP model ~ compliance.

    tau1 : float
        Parameter of KV-element of EVP model = retardation time.
    
    Returns
    -------
    Function
        Expression for fitting procedure.
    
    '''
    # Define model function (constants will be fixed later)
    func = const * (B0 + Cv*t - D1*np.exp(-t/tau1))
    # Return final function
    return(func)

def evp_s_d_2kv(t, const,B0,Cv, D1,D2, tau1,tau2):
    '''
    Function defining EVP model with elements [S + D + 2*KV].

    Parameters
    ----------
    t : float
        Time = X-variable for fitting.

    const : float
        Multiplicative constant for EVP models.
        * Tensile experiments: const = applied stress = sigma[GPa]
            - `sigma` = constant, calculated+saved in mcreep.const.Experiment
        * Indentation experiments: const = `F*K`
            - `F[mN]` = loading force
            - `K` = constant, calculated+saved in mcreep.const.Experiment

    B0 : float
        Parameter corresponding to S-element of EVP model.

    Cv : float
        Parameter corresponding to D-element of EVP model.

    D1, D2 : float, float
        Parameters of KV-elements of EVP model ~ compliance.

    tau1, tau2 : float, float
        Parameters of KV-elements of EVP model = retardation times.
    
    Returns
    -------
    Function
        Expression for fitting procedure.
    
    '''
    # Define model function (constants will be fixed later)
    func = const * (B0 + Cv*t - (D1*np.exp(-t/tau1) + D2*np.exp(-t/tau2)))
    # Return final function
    return(func)


def evp_s_d_3kv(t, const,B0,Cv, D1,D2,D3, tau1,tau2,tau3):
    '''
    Function defining EVP model with elements [S + D + 3*KV].

    Parameters
    ----------
    t : float
        Time = X-variable for fitting.

    const : float
        Multiplicative constant for EVP models.
        * Tensile experiments: const = applied stress = sigma[GPa]
            - `sigma` = constant, calculated+saved in mcreep.const.Experiment
        * Indentation experiments: const = `F*K`
            - `F[mN]` = loading force
            - `K` = constant, calculated+saved in mcreep.const.Experiment

    B0 : float
        Parameter corresponding to S-element of EVP model.

    Cv : float
        Parameter corresponding to D-element of EVP model.

    D1, D2, D3 : float, float, float
        Parameters of KV-elements of EVP model ~ compliance.

    tau1, tau2, tau3 : float, float, float
        Parameters of KV-elements of EVP model = retardation times.
    
    Returns
    -------
    Function
        Expression for fitting procedure.
    
    '''
    # Define model function (constants will be fixed later)
    func = const * (B0 + Cv*t - (
        D1*np.exp(-t/tau1) + D2*np.exp(-t/tau2) + D3*np.exp(-t/tau3)))
    # Return final function
    return(func)
