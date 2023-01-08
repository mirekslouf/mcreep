# Python package initialization file.

'''
MCREEP package
-------------
Fitting of various pre-defined models to creep data.

* The creep data:
    - a text file with two columns
    - the 1st column = time, the 2nd column = deformation
* The creep data may come from various experiments:
    - tensile creep => deformation(t) = strain(t) = epsilon(t)
    - indentation creep => deformation(t) = penetration_depth(t) = h(t)
* The pre-defined models within the pagkage are:
    - simple empirical models = Power law, Nutting's law
    - phenomenological models = elasto-visco-plastic models = EVP models
* Technical notes:
    - tensile creep is usually measured in macroscale
    - indentation creep is usually measured in micro/nanoscale
    - the models used in this package work for all three scales

Typical (simple/minimalistic but complete/working) example is as follows:
    
>>> # (0) Import modules
>>> import mcreep.const, mcreep.func, mcreep.model
>>>
>>> # (1) Define experiment, data and plot parameters
>>> EPAR = mcreep.const.Experiment(etype='Vickers', F=100)
>>> DPAR = mcreep.const.DataParameters(deformation_to_um=0.001)
>>> PPAR = mcreep.const.PlotParameters(xlabel='t [s]', ylabel='h [um]')
>>>
>>> # (2) Define model (using functions prepared in mcreep.func)
>>> MODEL = mcreep.model.Model(EPAR, DPAR, PPAR,
>>>     name='S+D+2KV', func=mcreep.func.evp_s_d_2kv, rtimes=[5,20])
>>>
>>> # (3) Run the whole calculation
>>> MODEL.run('pe1.txt', t_start=2.2, t_hold=99)
>>>
>>> # (4) Create report = save and print fitting results
>>> MODEL.final_report()
'''

__version__ = "1.1.2"
