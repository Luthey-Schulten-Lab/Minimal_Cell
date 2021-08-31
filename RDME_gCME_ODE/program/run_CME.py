"""
Author: Zane Thornburg
"""

import lm
from lm import GillespieDSolver
import sys

CSIMfilename = sys.argv[1]

lm.runSolver(CSIMfilename, 1, solver=GillespieDSolver(), cudaDevices=[0], checkpointInterval=0)
