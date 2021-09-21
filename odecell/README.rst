
ODECell
================

About
-----

**ODECell** provides the infrastructure to create large systems of ODEs, covering multiple cellular subsystems and 
mechanisms, using arbitrary sets rate forms, parameters, and species. The *modelbuilder* module presents 
an interface to define chemical species and concentrations, global parameters (such as pH), rate laws and reactions. 
Partial derivatives of rate laws can be provided for automatic jacobian calculation for every reaction.

The final model object can then be used by the *solver* module, producing a callable object that can be 
used directly by an ODE solver from SciPy.integrate.ode (see Tutorial notebook). 

If parameters were defined with value ranges (upper and lower bounds), the *paropt* module can be used
to optimize their values using multiple approaches (see Tutorial_Optimization notebook). Targets can be set
for both metabolite concentration and reaction fluxes, and the optimization strategy can be either called
through the *paropt* module (using a Differential Evolution implementation from SciPy) or or by an external
optimizer. In the latter case, the *paropt* module provides a "calcObjective" function that returns the cost
for a given parameter set.

Installation instructions (Linux/Unix)
------------------------

1. Create a Python virtual environment.

If you are not already working in a python or conda environment, to create an environment, one can use Anaconda.

Install Anaconda 3 if not already installed on your computer, download and follow the instructions at
https://www.anaconda.com/

To create a new environment named "envname", use

    conda create -n envname -c conda-forge python=3.7.3
    
Activate the environment

    conda activate envname
    
Change the pip version in the environment

    conda install pip==20.3.1

2. Install SUNDIALS in your system so that both headers and libraries are available, as well as its dependencies ( usually provided by development packages for `blas` and `lapack`).

Under Fedora 26 or newer, one can use the example below (as root):

    dnf install sundials-devel blas-devel lapack-devel

For Ubuntu 18.04.1 LTS, one can use (also as root):

    apt-get install libsundials-dev
    
Install SbTab and PyCVODES

    pip install sbtab==0.9.71
    
    conda install -c conda-forge pycvodes=0.14.0
    
Install other dependencies

    conda install -c conda-forge pydantic=1.8.2
    conda install -c conda-forge cobra
    conda install -c anaconda xlrd==1.2.0
    conda install -c anaconda mpi4py
    
3. Install ODECell.

Move into the repository.

    cd odecell

Use `pip` to install the package in your virtual environment.
    
    python -m pip install .
    
This will install all required python packages, and the install ODECell in your activated python virtual environment.

4. Finalize the installation.

Since jupyter extensions are installed as one of the requirements, run the following command to activate the extensions:
    
    jupyter contrib nbextension install --user

