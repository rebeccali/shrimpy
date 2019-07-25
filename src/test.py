#!/usr/bin/env python
"""
Testing for Shrimp Project
Rebecca Li 2019
"""

from mathUtil import testMathUtil
from aerodynamics import testAerodynamics
from shrimpClasses import testShrimpClasses
from bladeDynamics import testBladeDynamics
from shrimpController import testShrimpController
from shrimpOde import testShrimpOde
from runSim import testRunSim

print('Beginning Shrimp Project Testing')
testMathUtil()
testAerodynamics()
testShrimpClasses()
testBladeDynamics()
testShrimpController()
testShrimpOde()
testRunSim()
