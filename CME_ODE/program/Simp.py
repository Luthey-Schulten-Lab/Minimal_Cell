"""
Author: David Bianchi

The ODE Hook Solver to Integrate with a CME Simulation of JCVI Syn3A Minimal Cell Gene Expression
"""

# SBtab classes source code
from sbtab import SBtab

# SBtab validator
from sbtab import validatorSBtab

# Converter SBtab -> SBML
from sbtab import sbtab2sbml

# Converter SBML -> SBtab
from sbtab import sbml2sbtab

import libsbml
import csv
import pandas as pd
import numpy as np
from collections import defaultdict, OrderedDict
from scipy import integrate

# To access gene and protein data
from Bio import SeqIO
from Bio.Seq import Seq

import importlib

# Imports ODE model manager
import odecell
importlib.reload(odecell)

import Rxns as Rxns

### Related Constants

# Cell radius (meters):
r_cell = 200e-9 # 200 nm radius, 400 nm diameter

# Cell Volume (Liters):
cellVolume = ((4/3)*np.pi*(r_cell)**3)*(1000) # for a spherical cell

# Avogadro:
avgdr = 6.022E+23

# Convert particle counts to (mM) concentrations
countToMiliMol = 1000/(avgdr*cellVolume)

# Default metabolite concentration (mM)
defaultMetConcentration = 0.1

# Default Protein concentration (mM): set to 1 micro-Mol - like 2 proteins approx. 
defaultPtnConcentration = defaultMetConcentration/100

# Default RNA concentration (mM): 
rnaConcentration = round(10*countToMiliMol, 5)  # 10 molecules per cell

minPtnConc = round(10*countToMiliMol, 5)  # 10 molecules per cell

def initModel(pmap):
    """
    Initiate the model object and pass in the CME species counts

    Arguments:

    pmap (particle map): The CME particle Map
    
    Returns:

    upModel (odecell model object): The updated ODE kinetic model for simulation
    """

    # Initialize the ODECell model
    model = odecell.modelbuilder.MetabolicModel()
    
    
    zeroOrderOnOff = '$onoff * $K'
    model.zeroOrderOnOff = odecell.modelbuilder.RateForm(zeroOrderOnOff)
    model.updateAvailableForms()
    

    #print("Have a model object",flush=True)

    # Set verbosity outputs to zero for now to improve performance
    model.setVerbosity(0)


    # Define Rxns and pass in the Particle Map containing enzyme concentrations
    upModel = Rxns.defineRxns(model,pmap)

    return upModel



