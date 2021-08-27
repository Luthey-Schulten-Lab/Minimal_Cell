"""
A file to integrate out the ODE Simulation Model over a timestep

Author: David Bianchi
"""

from pycvodes import integrate_predefined
from pycvodes import integrate_adaptive
from scipy import integrate
import odecell
import numpy as np
import time as timer

### Constants
step = 0.1 # s
atol = 1e-6 # tolerance
rtol = 1e-6 # tolerance

def setSolver(model):
    """
    Set the solver for the model

    Parameters:
         
        (odecell.model) - the model object

    Returns:

        solvFunctor - a functor for the solver

    """

    ## We are NOT building for odeint (gives us more room to chose between CVODES and SciPy-ODE).
    ## We are NOT using a jacobian, since we do not have the partial derivatives for all rate forms.
    ## We are building with Cython for speed, this is a big model.

    # Builds the solver using a *Functor* interface
    solvFunctor = odecell.solver.ModelSolver(model) 
    solvFunctor.prepareFunctor() #OG uncomment

    # Set verbosity to 0 for now, below uncomment OG
    rxnIdList = solvFunctor.buildCall(odeint=False, useJac=False, cythonBuild=True, functor=True, verbose=0)
    #rxnIdList = solvFunctor.buildCall(odeint=True, useJac=False, cythonBuild=False, functor=False, transpJac=False, verbose=0, noBuild=True)

    # Sets up the actual solver, with updated parameter values
    #modelOptSpace, initParamVals=model.getOptSpace()
    initParamVals = model.getInitVals()
    solvFunctor = solvFunctor.functor( np.asarray(initParamVals, dtype=np.double) )

    return solvFunctor

### NOTE: Have to get a callable f(y,t) for scipy.ode without creating the functor
def noCythonSetSolver(model):
    """
    Set the solver without compiling via Cython

    Parameters:

    model (odecell Model object): The model object

    Returns:

    solver (odecell Solver object): The Solver object, to solve the system of ODEs representing metabolic reactions
    """

    # Construct a Model Solver Object
    solver = odecell.solver.ModelSolver(model)

    rxnIdList =solver.buildCall(verbose=0, useJac=False, transpJac=0, nocheck=False, odeint=False, cythonBuild=False, functor=False, noBuild=True)

    return solver 

def f_wrap(solv, t, y, dydt):

    #solv = setSolver(model)
    dydt[:] = solv(0,np.asarray(y))[:]

def getInitVals(model):
    y0=model.getInitVals()
    return y0

def runODE(y0,time,oldTime,ts,solv,model):
    """
    Run the ODE Model after getting initial conditions

    Parameters:

    y0 (seems non-necessary) - can remove
    time (float): the current hybrid simulation time
    oldTime (float): the last simulation time at which we stepped into the "hookSolver"
    ts (float): the timestep for the adaptive ODE Solver
    solv (odecell Solver object): The solver object, with call built
    model (odecell Model object): The model object

    Returns:

    results (np.array): the array containing ODE Simulation Results (Maybe only the last time should be passed?)
    """

    integrator = integrate.ode(solv)#, solv.calcJac)

    # TODO: Set which integrator to use (Appears vode is default), lsoda adaptive
    lsodaBool = True
    if (lsodaBool):
        integrator.set_integrator("lsoda")

    integrator.set_initial_value(model.getInitVals())

    ### With fixed timestepping
    step = ts
    delt = round(time-oldTime,2) # Get the deltaT we have accumulated between communication steps (May be greater than tau - user defined)
    totalTime = max(delt,ts) #time+
    results = np.empty((0,len(model.getInitVals())), float)

    startIntTime = timer.time()
    #print("Integration started at: ",startIntTime)

    while integrator.successful() and (integrator.t < totalTime) and (np.isclose([integrator.t],[totalTime])==([False])):
        #print("ode step time is: ", integrator.t)
        currConcentration = integrator.integrate(integrator.t + step)
        # Silence integrator output for now
        #print(integrator.t, currConcentration)
        results = np.append(results, [np.asarray(currConcentration)], axis=0 )

    endIntTime = timer.time()
    #print("Total Integration time was: ",endIntTime-startIntTime, " per step")
#     results = results[-1,:]

    # Return only the last timestep for the results, this is all thats really needed
    return results
