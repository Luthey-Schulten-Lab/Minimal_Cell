"""
A file to integrate out the ODE Simulation Model over a timestep

Author: David Bianchi
"""

from pycvodes import integrate_predefined
from pycvodes import integrate_adaptive
from scipy import integrate
from scipy.integrate import odeint as odeint
#from scipy.integrate import ode as integrate
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

def f_wrap(y, t, solv, dydt):

    #solv = setSolver(model)
    dydt[:] = solv(0,np.asarray(y))[:]
    #dydt[:] = solv(np.asarray(y))[:]
    return dydt

def getInitVals(model):
    y0=model.getInitVals()
#     print(y0)
    return y0

# def noCythonRunODE():
#     """
#     Run the ODE Model without compiling via Cython

#     Parameters:

#     Returns:
#     None
#     """
#     return 0

def runODE(y0,time,delt,ts,solv,model):
    """
    Run the ODE Model after getting initial conditions

    Parameters:

    y0 (seems non-necessary) - can remove
    time (float): the current hybrid simulation time
    delt (float): the communication timestep between stochastic and deterministic simulation
    ts (float): the timestep for the adaptive ODE Solver
    solv (odecell Solver object): The solver object, with call built
    model (odecell Model object): The model object

    Returns:

    results (np.array): the array containing ODE Simulation Results (Maybe only the last time should be passed?)
    """

    #y0 = model.getInitVals()
    #print("shape: ",len(y0))
    #y0 = np.asarray(y0,dtype=np.double)

    #modelOptSpace,initParamVals=model.getOptSpace()
    #dydt = np.zeros(len(y0))
    #tout, results, info = integrate_adaptive(f_wrap(solv,time+delt,y0,dydt),None,y0,time,time+delt,atol,rtol,dx0=ts,nsteps=10000)
    #tout, results, info = integrate_adaptive(f_wrap(solv,time,y0,dydt),None,y0,time,time+delt,atol,rtol,dx0=ts,nsteps=10000)

    #tout = 0.0
    #info = "place holder"

    #solv = solv.ModelSolver(model)
    #solv.buildCall(odeint=False, useJac=False, verbose=2)
    #integrator = integrate.ode(solv)#, solv.calcJac)

    # TODO: Set which integrator to use (Appears vode is default), lsoda adaptive
    #lsodaBool = True
    #if (lsodaBool):
        #integrator.set_integrator("lsoda")

    #integrator.set_initial_value(model.getInitVals())

    ### With fixed timestepping
    step = float(ts/100.0)
    totalTime = delt #time+
    resultsZero = np.zeros(len(model.getInitVals()))
    print("Len resultsZero is: ", len(resultsZero))

    startIntTime = timer.time()
    print("Integration started at: ",startIntTime)

    # With adaptive integration via LSODA: minStep = 0.1s, maxStep=1.0s
    #results=integrate.ode(integrator,min_step=0.1,max_step=delt)
    print("Solv type is :",type(solv))

    results = odeint(f_wrap,model.getInitVals(),np.linspace(time,time+delt,int(np.ceil(delt/ts)+1)),args = (solv,model.getInitVals()), hmin=1e-3, hmax=delt)
    #while integrator.successful() and integrator.t < totalTime:
        #currConcentration = integrator.integrate(integrator.t + step)
        # Silence integrator output for now
        #print(integrator.t, currConcentration)
        #results = np.append(results, [np.asarray(currConcentration)], axis=0 )
    #integrate(,min_step=0.1,max_step=ts)


    endIntTime = timer.time()
    print("Integration time ended at: ",endIntTime)
#     results = results[-1,:]

    # Return only the last timestep for the results, this is all thats really needed
    return results
