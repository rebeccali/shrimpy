#!/usr/bin/env python
"""
Shrimp Ode/equations of motion for Shrimp Project
Rebecca Li 2019
"""
import numpy as np

from bladeDynamics import getPropForceMoment
from shrimpController import shrimpController
from shrimpClasses import ShrimpState, defaultShrimpParams
from mathUtil import angVel2EulerAngleVel, rotateTensor, rotmFromYaw

# TODO:rebecca make this betterer
_odeOutputs = {
    "times": [],
    "thrustsProp": [],
    "forcesAero_w": [],
    "momentsAero_w": [],
    "forcesGravity_w": [],
    "forcesAngVel_w": [],
    "momentBodyProp_f": [],
    "momentShaftProp_f": [],
    "momentNegMotor_f": [],
    "momentGyro1_f": [],
    "momentGyro2_f": [],
    "momentBodyInertial_f": [],
    "momentPropInertial_f": [],
}

def flyerOde(odeState, t, p):
    """ This assembles everything for the ode integration and returns dx/dt in the Flyer Frame
        Arguments:
            odeState: the 1x16 np array containing the state vector.
            t: a time instant. Unused, since this is a time-invariant system.
            p: ShrimpParameters
            The ODE state vector is as follows:
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
        Returns the derivative of the ODE state vector.
        TODO: rename this function
    """
    _odeOutputs["times"].append(t)

    # Set up useful quantities
    s = ShrimpState.fromOdeState(odeState)  # ShrimpState
    totalMass = p.bodyMass + p.propMass
    # Rotations
    rot_f2b = s.rot_f2b
    rot_f2w = s.rot_f2w
    rot_w2f = np.transpose(rot_f2w)
    rot_b2f = np.transpose(rot_f2b)
    rot_b2p = rotmFromYaw(s.yaw_b2p)
    rot_p2b = np.transpose(rot_b2p)
    rot_p2f = rot_p2b.dot(rot_b2f)

    # Rotate some quantities
    vel_w2b_f = rot_w2f.dot(s.vel_w2b_w)
    # velocity of cg flyer is the same as velocity of cg body
    vel_b2f_f = np.zeros(3)
    vel_w2f_f = vel_w2b_f + vel_b2f_f
    propInertia_f = rotateTensor(p.propInertia_b, rot_b2f)
    inertia_f = rotateTensor(p.inertia_b, rot_b2f)

    # CALCULATIONS BELOW ###### ------------------------------------------------------------
    # Calculate controller
    pwm = shrimpController(p, s)

    # Calculate Aerodynamic Forces
    (bodyPropForce_b, bodyPropMoment_b) = getPropForceMoment(p.bodyPropParams, p, s)
    (shaftPropForce_b, shaftPropMoment_b) = getPropForceMoment(p.shaftPropParams, p, s)
    bodyPropMoment_f = rot_b2f.dot(bodyPropMoment_b)
    shaftPropMoment_f = rot_b2f.dot(shaftPropMoment_b)
    forcesAero_f = rot_b2f.dot(bodyPropForce_b + shaftPropForce_b)
    momentsAero_f = bodyPropMoment_f + shaftPropMoment_f
    _odeOutputs["thrustsProp"].append(forcesAero_f[2])
    _odeOutputs["forcesAero_w"].append(rot_f2w.dot(forcesAero_f))
    _odeOutputs["momentsAero_w"].append(rot_f2w.dot(momentsAero_f))


    # Calculate Gravity
    gravity_w = np.array([0, 0, -9.81])*totalMass
    gravity_f = rot_w2f.dot(gravity_w)
    _odeOutputs["forcesGravity_w"].append(gravity_w)

    # Motor equations
    K_t = p.motorParams.K_t
    motorR = p.motorParams.resistance
    battR = p.batteryParams.resistance
    battV = p.batteryParams.voltage
    totalResistance = max(1e-6, pwm * pwm * battR + motorR)  # Don't divide by zero
    motorCurrent = (battV * pwm - K_t * s.yawDot_b2p) / totalResistance
    motorMoment_f = np.array([0, 0, K_t * motorCurrent])  # TODO: double check this angle

    # Rotational Aspects
    # omg_rDot - prop z-rotational acc in prop frame
    angacc_b2p_f = np.linalg.inv(propInertia_f).dot((shaftPropMoment_f + motorMoment_f) * np.array([0, 0, 1]))
    # omg_bDot - body z-rotational acc in flyer frame
    angacc_f2b_f = np.linalg.inv(inertia_f).dot((bodyPropMoment_f - motorMoment_f) * np.array([0, 0, 1]))
    # MGyro1
    angvel_w2b_f = rot_b2f.dot(s.angvel_w2b_b)
    moment_gyro1_f = -np.cross(s.angvel_w2f_f, inertia_f.dot(angvel_w2b_f))
    # MGyro2
    angvel_f2p_f = rot_p2f.dot(np.array([0, 0, s.yawDot_f2b + s.yawDot_b2p]))
    angvel_w2p_f = s.angvel_w2f_f + angvel_f2p_f
    moment_gyro2_f = -np.cross(s.angvel_w2f_f, propInertia_f.dot(angvel_w2p_f))
    bodyInertialMoment_f = - inertia_f.dot(angacc_f2b_f)
    propInertialMoment_f = - propInertia_f.dot(angacc_b2p_f)
    _odeOutputs["momentBodyProp_f"].append(bodyPropMoment_f)
    _odeOutputs["momentShaftProp_f"].append(shaftPropMoment_f)
    _odeOutputs["momentNegMotor_f"].append(-motorMoment_f)
    _odeOutputs["momentGyro1_f"].append(moment_gyro1_f)
    _odeOutputs["momentGyro2_f"].append(moment_gyro2_f)
    _odeOutputs["momentBodyInertial_f"].append(bodyInertialMoment_f)
    _odeOutputs["momentPropInertial_f"].append(propInertialMoment_f)

    allMoments_f = (bodyPropMoment_f + shaftPropMoment_f - motorMoment_f +
                    moment_gyro1_f + moment_gyro2_f + bodyInertialMoment_f + propInertialMoment_f)


    allInertiaInv_f = np.linalg.inv(inertia_f + propInertia_f)
    # Piccoli's Euler equations
    angacc_w2f_f = allInertiaInv_f.dot(allMoments_f)
    velEuler_w2f = angVel2EulerAngleVel(s.angvel_w2f_f, s.euler_w2f)

    # Newton's equations
    forcesAngVel_f = - np.cross(s.angvel_w2f_f, vel_w2f_f)
    acc_w2f_f = (1 / totalMass) * (forcesAero_f + gravity_f) + forcesAngVel_f
    _odeOutputs["forcesAngVel_w"].append(rot_f2w.dot(forcesAngVel_f))

    acc_w2f_w = rot_f2w.dot(acc_w2f_f)
    acc_f2b_w = np.zeros(3)  # flyer and body occupy the same point bruh
    acc_w2b_w = acc_w2f_w + acc_f2b_w

    dOdeState = np.zeros(16)
    dOdeState[0] = s.vel_w2b_w[0]  # x body position in world frame
    dOdeState[1] = s.vel_w2b_w[1]  # y body position in world frame
    dOdeState[2] = s.vel_w2b_w[2]  # z body position in world frame
    dOdeState[3] = acc_w2b_w[0]  # x body velocity in world frame
    dOdeState[4] = acc_w2b_w[1]  # y body velocity  in world frame
    dOdeState[5] = acc_w2b_w[2]  # z body velocity  in world frame
    dOdeState[6] = velEuler_w2f[0]  # Z ZXY euler angle in flyer frame (yaw)
    dOdeState[7] = velEuler_w2f[1]  # Y ZXY euler angle in flyer frame
    dOdeState[8] = velEuler_w2f[2]  # X ZXY euler angle in flyer frame
    dOdeState[9] = angacc_w2f_f[0]  # p is x roll rate in flyer frame
    dOdeState[10] = angacc_w2f_f[1]  # q is y pitch rate in flyer frame
    dOdeState[11] = 0  # r is z yaw rate in flyer frame, would be duplicative if not zero
    dOdeState[12] = angacc_f2b_f[2]  # yawDot_f2b is body yaw rate wrt flyer frame
    dOdeState[13] = angacc_b2p_f[2]  # yawDot_b2p is the prop speed wrt body frame (stator frame)
    dOdeState[14] = s.yawDot_f2b  # yaw_f2b from flyer to stator/body frame
    dOdeState[15] = s.yawDot_b2p  # yaw_b2p from stator/body to rotor/prop frame
    return dOdeState


def testShrimpOde():
    params = defaultShrimpParams()
    t = 0
    x = np.zeros(16)
    dX = flyerOde(x, t, params)
    print(dX)
