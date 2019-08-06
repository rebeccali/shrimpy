#!/usr/bin/env python
"""
Classes for Shrimp Project
Rebecca Li 2019
"""

from enum import Enum
import numpy as np
from scipy.interpolate import interp1d

from mathUtil import addYaw, euler2Rotm, rotm2Euler, rotmFromYaw, rpm2RadiansPerSecond


class PropellerType(Enum):
    SHAFT = 1
    BODY = 2


class BladeElementParameters:
    def __init__(self, rho, pitch, width, radius, chord, index, numBlades,
                 height_b2e, inflowVel, yawDot_b2e):
        self.rho = rho
        self.pitch = pitch
        self.width = width  # this is in the radial direction (strip model width)
        self.radius = radius  # this should probably "middle" of the overall radius width
        self.chord = chord
        self.index = index  # Blade Index out of n blades
        self.numBlades = numBlades  # total number of blades. A normal prop has 2
        self.height_b2e = height_b2e  # scalar [m]
        self.inflowVel = inflowVel  # Inflow velocity for this element. Scalar [m/s]
        self.yawDot_b2e = yawDot_b2e  # spin velocity of prop about body z axis for this element. Scalar [m/s]


class PropellerParameters:
    """Parameters of all blades in the propeller
       Assumes propeller center is aligned in the z-axis with the body center of mass.
    """

    def __init__(self, numBlades, radiusRootTip, getPitchFromRadius, getChordFromRadius,
                 height_b2p, propType):
        self.numBlades = numBlades
        self.radiusRootTip = radiusRootTip  # Tuple, (root, tip)
        self.height_b2p = height_b2p  # Z-axis distance from cg to propeller, +z is up
        self.propType = propType  # PropellerType Enum
        self.getPitchFromRadius = getPitchFromRadius  # Function that takes radius and returns pitch at that point
        self.getChordFromRadius = getChordFromRadius  # Function that takes radius and returns pitch at that point

    @classmethod
    def fromRootTipParams(cls, numBlades, pitchRootTip, radiusRootTip, chordRootTip,
                          height_b2p, propType):
        getPitchFromRadius = interp1d(radiusRootTip, pitchRootTip)
        getChordFromRadius = interp1d(radiusRootTip, chordRootTip)
        return cls(numBlades, radiusRootTip, getPitchFromRadius, getChordFromRadius,
                   height_b2p, propType)


class BatteryParameters:
    """Parameters for a battery"""

    def __init__(self, voltage, resistance):
        self.voltage = voltage  # nominal voltage
        self.resistance = resistance  # internal resistance, assumed constant heh


class MotorParameters:
    """Parameters for a motor"""

    def __init__(self, K_t, resistance):
        self.K_t = K_t  # Eletromotive forceconstant, or 1/K_v
        self.resistance = resistance  # stator resistance, assumed constant heh


class ShrimpVizParameters:
    """Parameters for a visualization """

    def __init__(self, bodyWidth, bodyHeight, propDiscRadius, propDiscThickness):
        self.propDiscRadius = propDiscRadius
        self.propDiscThickness = propDiscThickness
        self.bodyWidth = bodyWidth
        self.bodyHeight = bodyHeight


class ShrimpParameters:
    """ Parameters of Shrimp Aircraft.
    """

    def __init__(self, bodyPropParams, shaftPropParams, bodyMass, Ixx, Iyy, Izz,
                 propIxx, propIyy, propIzz, rho, motorParams, batteryParams, propMass, vizParams):
        self.bodyPropParams = bodyPropParams
        self.shaftPropParams = shaftPropParams
        self.bodyMass = bodyMass  # mass of everything attached to stator, so not rotor and prop
        self.inertia_b = np.diag([Ixx, Iyy, Izz])  # body Inertia without prop inertia
        self.propInertia_b = np.diag([propIxx, propIyy, propIzz])
        self.rho = rho  # Air density
        self.motorParams = motorParams
        self.batteryParams = batteryParams
        self.propMass = propMass
        self.vizParams = vizParams


def getRotsFromEuler(euler_w2b, yaw_f2b):
    rot_f2b = rotmFromYaw(yaw_f2b)
    rot_w2b = euler2Rotm(euler_w2b)
    rot_b2w = np.transpose(rot_w2b)
    rot_f2w = rot_f2b.dot(rot_b2w)
    return (rot_f2b, rot_w2b, rot_b2w, rot_f2w)


class ShrimpState:
    """ Shrimp state object. Stores state in a clear format, and also converts
        the clear format to a numpy vector for ODE integration.
    """

    def __init__(self, r_w2b_w, vel_w2b_w, euler_w2b, angvel_w2b_b, inflowVel,
                 yaw_b2p, yaw_f2b, yawDot_b2p, yawDot_f2b, rot_f2b, rot_b2w, rot_w2b, rot_f2w,
                 euler_w2f, angvel_w2f_f, angvel_f2b_f):
        self.r_w2b_w = r_w2b_w  # np.array 3x1 position vector from world origin to body/cg
        self.vel_w2b_w = vel_w2b_w  # np.array 3x1 velocity vector from world origin to body/cg
        self.euler_w2b = euler_w2b  # np.array 3x1 euler from world to body frame
        self.angvel_w2b_b = angvel_w2b_b  # np.array 3x1 angular velocity of the body in the body frame
        self.inflowVel = inflowVel  # scalar inflow velocity adjustment due to duct on the shaft propeller [m/s]
        self.yaw_b2p = yaw_b2p  # the angle from the stator/body to the rotor/propeller
        self.yaw_f2b = yaw_f2b  # the angle from the flyer frame to the body frame
        # the speed of the spin of the propeller relative to the body/stator. Scalar [rad/s]
        self.yawDot_b2p = yawDot_b2p
        self.yawDot_f2b = yawDot_f2b  # the speed of the spin of the body relative to the flyer frame. Scalar [rad/s]
        self.rot_f2b = rot_f2b
        self.rot_b2w = rot_b2w
        self.rot_w2b = rot_w2b
        self.rot_f2w = rot_f2w
        self.euler_w2f = euler_w2f
        self.angvel_w2f_f = angvel_w2f_f
        self.angvel_f2b_f = angvel_f2b_f

    @classmethod
    def fromVecState(cls, r_w2b_w, vel_w2b_w, euler_w2b, angvel_w2b_b, inflowVel,
                     yaw_b2p, yaw_f2b, yawDot_b2p, yawDot_f2b):
        (rot_f2b, rot_w2b, rot_b2w, rot_f2w) = getRotsFromEuler(euler_w2b, yaw_f2b)
        euler_w2f = rotm2Euler(np.transpose(rot_f2w))
        rot_b2f = np.transpose(rot_f2b)
        angvel_f2b_f = np.array([0, 0, yawDot_f2b])
        angvel_w2b_f = rot_b2f.dot(angvel_w2b_b)
        angvel_b2f_f = -angvel_f2b_f
        angvel_w2f_f = angvel_w2b_f + angvel_b2f_f
        return cls(r_w2b_w, vel_w2b_w, euler_w2b, angvel_w2b_b, inflowVel,
                   yaw_b2p, yaw_f2b, yawDot_b2p, yawDot_f2b, rot_f2b, rot_b2w, rot_w2b, rot_f2w,
                   euler_w2f, angvel_w2f_f, angvel_f2b_f)

    @classmethod
    def fromOdeState(cls, odeState):
        """ Converts from vector to this object.
            odeState[0] = x body position in world frame
            odeState[1] = y body position in world frame
            odeState[2] = z body position in world frame
            odeState[3] = x body velocity in world frame
            odeState[4] = y body velocity  in world frame
            odeState[5] = z body velocity  in world frame
            odeState[6] = Z ZXY euler angle in flyer frame (yaw)
            odeState[7] = Y ZXY euler angle in flyer frame
            odeState[8] = X ZXY euler angle in flyer frame
            odeState[9] = p is x roll rate in flyer frame
            odeState[10] = q is y pitch rate in flyer frame
            odeState[11] = r is z yaw rate in flyer frame
            odeState[12] = yawDot_f2b is body yaw rate wrt flyer frame
            odeState[13] = yawDot_b2p is the prop speed wrt body frame (stator frame)
            odeState[14] = yaw_f2b from flyer to stator/body frame
            odeState[15] = yaw_b2p from stator/body to rotor/prop frame
        """
        r_w2b_w = odeState[0:3]
        vel_w2b_w = odeState[3:6]
        euler_w2f = odeState[6:9]
        angvel_w2f_f = odeState[9:12]
        yawDot_f2b = odeState[12]
        yawDot_b2p = odeState[13]
        yaw_f2b = odeState[14]
        yaw_b2p = odeState[15]
        euler_w2b = addYaw(euler_w2f, yaw_f2b)
        (rot_f2b, rot_w2b, rot_b2w, rot_f2w) = getRotsFromEuler(euler_w2b, yaw_f2b)
        angvel_f2b_f = np.array([0, 0, yawDot_f2b])  # Angvel of the in flyer frame
        angvel_w2b_f = angvel_w2f_f + angvel_f2b_f
        angvel_w2b_b = rot_f2b.dot(angvel_w2b_f)
        inflowVel = 0
        return cls(r_w2b_w, vel_w2b_w, euler_w2b, angvel_w2b_b, inflowVel,
                   yaw_b2p, yaw_f2b, yawDot_b2p, yawDot_f2b, rot_f2b, rot_b2w, rot_w2b, rot_f2w,
                   euler_w2f, angvel_w2f_f, angvel_f2b_f)


def defaultShaftPropParams():
    """ Parameters for Cheerson
        Data from propeller_cheerson_cx10.m
    """
    numBlades = 2
    # TODO: do pitches/radius with interpolable sections
    pitchRootTip = (0.175, 0.18)  # In code this is Beta
    chordRootTip = (3.22e-3, 4.12e-3)
    radiusRootTip = (0, 0.0146)
    height_b2p = -1.55 / 1000  # Piccolissimo_V11, positive is UP
    propType = PropellerType.SHAFT
    return PropellerParameters.fromRootTipParams(numBlades, pitchRootTip, radiusRootTip,
                                                 chordRootTip, height_b2p, propType)


def defaultBodyPropParams():
    """ Parameters for body_Piccolissimo_V11"""
    numBlades = 6
    # TODO: do pitches/radius with interpolable sections
    pitchRootTip = (0.15359, 0.7927)
    chordRootTip = (0.0135, 0.0054)
    radiusRootTip = (0, 19.16e-3)
    height_b2p = 2.56 / 1000  # Piccolissimo_V11, positive is UP
    propType = PropellerType.BODY
    return PropellerParameters.fromRootTipParams(numBlades, pitchRootTip, radiusRootTip,
                                                 chordRootTip, height_b2p, propType)


def defaultBatteryParams():
    """ Default battery Parameters
        From PiccolissimoControl.m
    """
    voltage = 3.7
    resistance = 2
    return BatteryParameters(voltage, resistance)


def defaultMotorParams():
    """ From motor_cheerson_cx10.m"""
    resistance = 1.97
    K_t = 1 / (20500 * 2 * np.pi / 60)  # Nm/A
    return MotorParameters(K_t, resistance)


def defaultShrimpParams():
    """ Piccolissimo V11 essentially"""
    bodyPropParams = defaultBodyPropParams()
    shaftPropParams = defaultShaftPropParams()
    massBody = 4.22e-3
    inertiaAdjustment = 2.5
    Ixx = inertiaAdjustment * 383e-9
    Iyy = inertiaAdjustment * 439e-9
    Izz = inertiaAdjustment * 697e-9

    # Taken from propeller_cheerson_cx10.m Though this describes the rotor mass/inertia
    # as well
    propIxx = 3.95e-9
    propIyy = 6.31e-9
    propIzz = 3.21e-9
    propMass = 0.25e-3

    rho = 1.225
    motorParams = defaultMotorParams()
    batteryParams = defaultBatteryParams()

    propDiscRadius = shaftPropParams.radiusRootTip[1]
    c = shaftPropParams.getChordFromRadius(propDiscRadius)
    p = shaftPropParams.getPitchFromRadius(propDiscRadius)
    propDiscThickness = c * np.cos(p)
    bodyWidth = 0.005
    bodyHeight = 0.003
    vizParams = ShrimpVizParameters(bodyWidth, bodyHeight, propDiscRadius, propDiscThickness)

    return ShrimpParameters(bodyPropParams, shaftPropParams, massBody, Ixx, Iyy, Izz,
                            propIxx, propIyy, propIzz, rho, motorParams, batteryParams, propMass, vizParams)

def zeroShrimpState():
    """ Generates Shrimp state where everything is zero/nominal """
    odeState = np.zeros(16)
    state = ShrimpState.fromOdeState(odeState)
    return state

def dummyShrimpState():
    r_w2b_w = np.array([6, 20, 1])
    vel_w2b_w = np.array([1, 2, 3])
    euler_w2b = np.array([2, 3, 5]) * np.pi / 180.
    angvel_w2b_b = np.array([1, 23, 8])
    inflowVel = 0
    yaw_b2p = 0.2
    yaw_f2b = 0.3
    yawDot_f2b = rpm2RadiansPerSecond(10)
    yawDot_b2p = rpm2RadiansPerSecond(500)
    state = ShrimpState.fromVecState(r_w2b_w, vel_w2b_w, euler_w2b, angvel_w2b_b, inflowVel,
                                     yaw_b2p, yaw_f2b, yawDot_b2p, yawDot_f2b)
    return state


def dummyOdeStateVsShrimpState():
    s = dummyShrimpState()
    odeState = np.zeros(16)
    odeState[0:3] = s.r_w2b_w
    odeState[3:6] = s.vel_w2b_w
    odeState[6:9] = s.euler_w2f
    odeState[9:12] = s.angvel_w2f_f
    odeState[12] = s.yawDot_f2b
    odeState[13] = s.yawDot_b2p
    odeState[14] = s.yaw_f2b
    odeState[15] = s.yaw_b2p
    state = ShrimpState.fromOdeState(odeState)
    return (s, state)


def testShrimpClasses():
    temp = dummyShrimpState()
    temp2 = np.zeros(16)
    temp3 = ShrimpState.fromOdeState(temp2)
    temp4 = defaultShrimpParams()
    (s1, s2) = dummyOdeStateVsShrimpState()
    print('The following should be close to zero, if not zero: ')
    print('s.r_w2b_w: ', (s1.r_w2b_w - s2.r_w2b_w))
    print('s.vel_w2b_w: ', (s1.vel_w2b_w - s2.vel_w2b_w))
    print('s.euler_w2f: ', (s1.euler_w2f - s2.euler_w2f))
    print('s.angvel_w2f_f: ', (s1.angvel_w2f_f - s2.angvel_w2f_f))
    print('s.yawDot_f2b: ', (s1.yawDot_f2b - s2.yawDot_f2b))
    print('s.yawDot_b2p: ', (s1.yawDot_b2p - s2.yawDot_b2p))
    print('s.yaw_f2b: ', (s1.yaw_f2b - s2.yaw_f2b))
    print('s.yaw_b2p: ', (s1.yaw_b2p - s2.yaw_b2p))
    return (temp, temp2, temp3, temp4)
