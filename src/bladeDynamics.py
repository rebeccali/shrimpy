#!/usr/bin/env python
"""
Blade Element Theory for Shrimp Project
Rebecca Li 2019
"""
import numpy as np

from aerodynamics import getAngleOfAttack, getLiftDragElement, getRelativeWind
from shrimpClasses import BladeElementParameters, dummyShrimpState, PropellerType
from shrimpClasses import defaultBodyPropParams, defaultShrimpParams, zeroShrimpState  # TODO make this import prettier
from mathUtil import rpm2RadiansPerSecond, rotPerSec2RadiansPerSecond


def getElementForceMoment(params, state):
    """Returns force and moment vectors in the flyer frame
    Assumes the element radius is the radius at the center of the width.
    The width is only used in a surface area calculation.
    Arguments:
        params (BladeElementParameters)
        state (ShrimpState)
    """
    # TODO: refactor with shorter names
    rho = params.rho
    bladeIndex = params.index
    bladeRadius = params.radius
    vel_w2b_w = state.vel_w2b_w
    angvel_w2b_b = state.angvel_w2b_b
    inflowVelocity = params.inflowVel
    height_b2e = params.height_b2e
    yawDot_b2e = params.yawDot_b2e
    numBlades = params.numBlades
    rot_w2b = state.rot_w2b

    # TODO: make call less ugly. Maybe move r_b2e_b and rot_b2e into the creation of the blade element object
    (relativeWind_e, gamma, r_b2e_b, rot_b2e) = getRelativeWind(bladeIndex, numBlades,
                                                                bladeRadius, vel_w2b_w,
                                                                rot_w2b, angvel_w2b_b,
                                                                inflowVelocity, height_b2e,
                                                                yawDot_b2e)

    alpha = getAngleOfAttack(gamma, params.pitch)
    (lift, drag) = getLiftDragElement(params.width, params.chord, rho, relativeWind_e, alpha)

    # Forces in the propeller frame
    normalForce_e = lift * np.cos(gamma) + drag * np.sin(gamma)
    tangentialForce_e = -lift * np.sin(gamma) + drag * np.cos(gamma)

    # Z-axis is up!
    forces_e = np.array([-tangentialForce_e, 0, normalForce_e])

    rot_e2b = np.transpose(rot_b2e)
    forces_b = rot_e2b.dot(forces_e)
    moments_b = np.cross(r_b2e_b, forces_b)

    return(forces_b, moments_b)


def testBladeElement(radius, chord, pitch, width):
    rho = 1.225
    bladeIndex = 1
    numBlades = 2
    height_b2e = 0.01

    inflowVel = 1
    yawDot_b2e = rpm2RadiansPerSecond(500)
    # Set up this particular propeller element
    bladeParams = BladeElementParameters(rho, pitch, width, radius, chord, bladeIndex,
                                         numBlades, height_b2e, inflowVel, yawDot_b2e)

    state = dummyShrimpState()
    (forces_b, moments_b) = getElementForceMoment(bladeParams, state)
    return (forces_b, moments_b)


def getBladeForceMoment(propParams, shrimpParams, state, bladeIndex):
    """Get the force moment of a single blade.
       Splits up the blade into n sections and evaluates a strip model for each
       section. Also decides inflow velocity and angular velocity based off of whether
       it is the body or shaft propeller
       Arguments:
           propParams: PropellerParameters
           shrimpParams: ShrimpParameters
           state: ShrimpState
           bladeIndex: index of the blade in the propeller
    """
    n = 5  # The number of sections to split the blade into for element theory
    (radiusRoot, radiusTip) = propParams.radiusRootTip
    width = (radiusTip - radiusRoot) / n

    radii = np.linspace(radiusRoot, radiusTip, n)
    pitches = propParams.getPitchFromRadius(radii)
    chords = propParams.getChordFromRadius(radii)

    # Assign propeller specific characteristics
    if propParams.propType == PropellerType.SHAFT:
        inflowVel = state.inflowVel
        yawDot_b2p = state.yawDot_b2p
    elif propParams.propType == PropellerType.BODY:
        inflowVel = 0
        yawDot_b2p = 0
    else:
        raise Exception('Using unimplemented PropellerType!')

    def getThisElemForceMoment(pitch, radius, chord):
        elementParams = BladeElementParameters(shrimpParams.rho, pitch, width, radius,
                                               chord, bladeIndex, propParams.numBlades,
                                               propParams.height_b2p, inflowVel,
                                               yawDot_b2p)
        return getElementForceMoment(elementParams, state)

    forceMomentTuples_b = [getThisElemForceMoment(p, r, c)
                           for (p, r, c) in zip(pitches, radii, chords)]
    forces_b = sum([f for f, _ in forceMomentTuples_b])
    moments_b = sum([m for _, m in forceMomentTuples_b])
    return (forces_b, moments_b)


def getPropForceMoment(propParams, shrimpParams, state):
    """Computes the force across all blades in the propeller
        Arguments:
        propParams: PropellerParameters to be iterated over
        shrimpParams: ShrimpParameters
        state: ShrimpState
    """
    forceMomentTuples_b = [getBladeForceMoment(propParams, shrimpParams, state, i)
                           for i in range(propParams.numBlades)]
    forces_b = sum([f for f, _ in forceMomentTuples_b])
    moments_b = sum([m for _, m in forceMomentTuples_b])
    return forces_b, moments_b


def testPropForceMoment():
    propParams = defaultBodyPropParams()
    shrimpParams = defaultShrimpParams()
    state = dummyShrimpState()

    (f, m) = getPropForceMoment(propParams, shrimpParams, state)
    print('Forces: ', f)
    print('Moments: ', m)


def testBladeDynamics():
    n = 5  # The number of sections to split the blade into for element theory
    radiusRoot = 0
    radiusTip = 0.1  # maximum blade radius / tip to root
    width = (radiusTip - radiusRoot) / n
    chordRoot = 0.1  # Chord at the root of the blade
    chordTip = 0.2  # Chord at the tip of the blade
    pitchRoot = 0 * np.pi / 180
    pitchTip = 15 * np.pi / 180

    radii = np.linspace(radiusRoot, radiusTip, n)
    chords = np.linspace(chordRoot, chordTip, n)
    pitches = np.linspace(pitchRoot, pitchTip, n)

    forceMomentTuples = [testBladeElement(r, c, p, width) for (r, c, p) in zip(radii, chords, pitches)]

    forces = sum([f for f, _ in forceMomentTuples])
    moments = sum([m for _, m in forceMomentTuples])
    print('Forces', forces)
    print('Moments', moments)

    testPropForceMoment()
