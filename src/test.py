#!/usr/bin/env python
"""
Testing for Shrimp Project
Rebecca Li 2019
"""

import os

from mathUtil import testMathUtil
from aerodynamics import testAerodynamics
from shrimpClasses import testShrimpClasses
from bladeDynamics import testBladeDynamics
from shrimpController import testShrimpController
from shrimpOde import testShrimpOde
from runSim import testRunSim
from shrimpVisualizer import testShrimpVisualizer
from calculateThrustTorqueProp import testThrustTorqueProp

print('Beginning Shrimp Project Testing')
testMathUtil()
testAerodynamics()
testShrimpClasses()
testBladeDynamics()
testShrimpController()
testShrimpOde()
testRunSim()
testShrimpVisualizer()
testThrustTorqueProp()
print('All Tests Passed.')
# Note: can't use sys exit because vpython does not respond to that. Go figure.
os._exit(0)
