#!/usr/bin/env python
"""
Script to calculate the thrust and drag of a propeller. Do this by calculating the aerodynamic forces
and moments at still, then converting these to thrust and drag.
Rebecca Li 2019
"""
from bladeDynamics import getPropForceMoment
from mathUtil import rotPerSec2RadiansPerSecond
from shrimpClasses import PropellerType, zeroShrimpState, PropellerParameters, defaultShrimpParams


def testThrustTorqueProp(propParams=None):
    # TODO: grab from spreadsheet
    # Set up parameters
    rho = 1.225
    rotorFreqHz = 30

    # Calculations
    if propParams is None:
        numBlades = 2
        pitchRootTip = (0.174, 0.174)  # In code this is Beta
        chordRootTip = (0.1, 0.1)
        radiusRootTip = (0, 1)
        height_b2p = 0.1
        propType = PropellerType.SHAFT
        propParams = PropellerParameters.fromRootTipParams(numBlades, pitchRootTip, radiusRootTip, chordRootTip,
                                                           height_b2p, propType)

    shrimpParams = defaultShrimpParams()  # Only rho matters
    shrimpParams.rho = rho

    state = zeroShrimpState()
    state.yawDot_b2p = rotPerSec2RadiansPerSecond(rotorFreqHz)

    (forces_b, moments_b) = getPropForceMoment(propParams, shrimpParams, state)
    print('Thrust [N]: %f' % forces_b[2])
    print('Drag Moment [Nm]: %f' % moments_b[2])
    print('\n')

    print('Forces X: %f' % forces_b[0])
    print('Forces Y: %f' % forces_b[1])
    print('Forces Z: %f' % forces_b[2])
    print('Moments X: %f' % moments_b[0])
    print('Moments Y: %f' % moments_b[1])
    print('Moments Z: %f' % moments_b[2])

if __name__ == "__main__":
    testThrustTorqueProp()
