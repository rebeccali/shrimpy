#!/usr/bin/env python
"""
Aerodynamics for Shrimp Project
Rebecca Li 2019
"""

import numpy as np
from scipy.spatial.transform import Rotation as R

from mathUtil import euler2Rotm


def getRelativeWind(bladeIndex, numBlades, bladeRadius, vel_w2b_w,
                    rot_w2b, angvel_w2b_b, inflowVelocity, height_b2e,
                    yawDot_b2e):
    """
    Returns relative wind angle as well as the wind angle vector for
    the given quantities

    Return values:
    relativeWind_e : the relative wind in the propellor frame
    relWindAngleGamma: the wind angle relative to the propellor blade disk angle
    r_b2e_b: the vector from the body to the element in the body frame
    rot_b2e: the rotation matrix from the body to element frame

    Arguements:
    bladeIndex: the index of the blade to be used. Assumes blade 0 is aligned with the
    x-axis of the body
    numBlades: total number of blades, assuming equal spacing
    bladeRadius: the radius to do this overall calculation
    vel_w2b_w: velocity of the body relative to the world velocity in the world frame
    rot_w2b: DCM of world to body
    angvel_w2b_b: angular velocity of body relative to world velocity in body axes (pqr)
    inflowVelocity: also called nu - scalar amount of inflow velocity
    height_b2e: the distance from body to element along the body z axis (up is positive)
    yawDot_b2e: speed of element relative to body about the body z-axis (spin rate)

    Notation:
    _e is prop element frame
    _b is the body frame (about the center of gravity)
    _w is world/inertial frame

    """
    # Calculate Rotation Matricies
    propellor_angle = bladeIndex / numBlades * 2 * np.pi
    rot_b2e = R.from_euler('Z', propellor_angle + np.pi / 2).as_dcm()

    # Body Velocity Calculation
    vel_w2b_b = rot_w2b.dot(vel_w2b_w)
    vel_w2b_e = rot_b2e.dot(vel_w2b_b)

    # Inflow velocity Calculation
    # it is always in the negative z direction
    vel_inducedInflow_e = np.array([0, 0, -inflowVelocity])

    # Velocity from rotational component Calculation
    # Calculate the vector from the cg to the element
    pos_x_blade = bladeRadius * np.cos(propellor_angle + np.pi / 2)
    pos_y_blade = bladeRadius * np.sin(propellor_angle + np.pi / 2)
    r_b2e_b = np.array([pos_x_blade, pos_y_blade, height_b2e])

    # Calculate the angular velocity in the body frame
    angvel_b2e_b = np.array([0, 0, yawDot_b2e])
    # Total angular velocity of the propellor is body + prop
    angvel_w2e_b = angvel_b2e_b + angvel_w2b_b

    vel_b2e_b = np.cross(angvel_w2e_b, r_b2e_b)
    vel_b2e_e = rot_b2e.dot(vel_b2e_b)

    # Sum all forces together
    relativeWind_e = -(vel_w2b_e + vel_b2e_e) + vel_inducedInflow_e
    relWindAngleGamma = np.arctan2(relativeWind_e[2], relativeWind_e[0])

    return (relativeWind_e, relWindAngleGamma, r_b2e_b, rot_b2e)


def testRelativeWind():
    """Test for relative wind. """
    eul_w2b = np.pi / 180 * np.array([0, 0, 0])  # ZXY
    rot_w2b = euler2Rotm(eul_w2b)

    numBlades = 2
    bladeRadius = 0.1
    bladeIndex = 1
    vel_w2b_w = np.array([0, 10, 10])
    angvel_w2b_b = np.array([0.0, 0.0, 0.0])
    inflowVelocity = 0
    height_b2e = 2
    yawDot_b2e = 0

    # Using XYZ euler angles such that the last angle Z is the correct psi angle
    (relWind, gamma, _, _) = getRelativeWind(bladeIndex, numBlades, bladeRadius, vel_w2b_w,
                                             rot_w2b, angvel_w2b_b, inflowVelocity, height_b2e,
                                             yawDot_b2e)

    print('gamma', gamma * 180 / np.pi)
    print('relwind', relWind)


def getAngleOfAttack(bladePitch, relWindAngleGamma):
    """ Returns angle of attack between -pi and pi"""
    unwrapped_angle = bladePitch + relWindAngleGamma
    return (unwrapped_angle + np.pi) % (2 * np.pi) - np.pi


def clFlatPlate(alpha):
    """ Coeff of lift for a Flat Plate """
    return 2. * np.sin(alpha) * np.cos(alpha)


def cdFlatPlate(alpha):
    """ Coeff of drag for a Flat Plate """
    return 2. * np.sin(alpha)**2.0


def cmFlatPlate(alpha):
    """ Coeff of pitch moment for a Flat Plate """
    return 0.0


def getLiftDragElement(elementWidth, chord, rho, relWind_e, alpha):
    """ Returns magnitude of lift and drag
    Arguements:
    elementWidth: the width of the blade section
    chord: airfoil chord
    rho: air density
    relWind_e : relative wind in the propellor frame
    alpha: angle of attack of section
    Returns:
    (lift, drag) scalars
    """
    # velocity only incorporates the x and z components in the wind frame
    speed = relWind_e.dot([1, 0, 1])
    coeffLift = clFlatPlate(alpha)
    coeffDrag = cdFlatPlate(alpha)
    quantity = 0.5 * speed * speed * chord * elementWidth * rho
    lift = quantity * coeffLift
    drag = quantity * coeffDrag
    return (lift, drag)


def testAerodynamics():
    """Test Aerodynamics!"""
    testRelativeWind()
