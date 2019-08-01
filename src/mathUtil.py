#!/usr/bin/env python
"""
Math utilities for Shrimp Project
Rebecca Li 2019
"""

import numpy as np
from scipy.spatial.transform import Rotation as R

_eulerConvention = 'ZXY'


def rotm2Euler(x):
    """Rotation matrix to Euler Angles ZXY intrinsic vector for Shrimp"""
    r = R.from_dcm(x)
    return r.as_euler(_eulerConvention)


def euler2Rotm(x):
    """Euler angles ZXY intrinsic vector to Rotation Matrix for Shrimp"""
    r = R.from_euler(_eulerConvention, x)
    return r.as_dcm()


def quat2Euler(x):
    """Quaterion to Euler Angles ZXY intrinsic vector for Shrimp"""
    r = R.from_quat(x)
    return r.as_euler(_eulerConvention)


def euler2Quat(x):
    """Euler angles ZXY intrinsic vector to Quaternion for Shrimp"""
    r = R.from_quat(_eulerConvention, x)
    return r.as_euler()


def addYaw(eulerAngles, yaw):
    """ Add yaw in a safe way since we want to control ZXY from here"""
    return eulerAngles + np.array([yaw, 0, 0])


def rotmFromYaw(yaw_a2b):
    """ Returns rot_a2b rotation matrix about z-axis from a to b """
    return R.from_euler('z', -yaw_a2b).as_dcm()


def getYaw(eulerAngle):
    return eulerAngle[0]


def getPitch(eulerAngle):
    return eulerAngle[1]


def getRoll(eulerAngle):
    return eulerAngle[2]


def eulerExtXYZfromEulerShrimp(x):
    """Converts euler extrinsic from eulerShrimp"""
    eulerExtrinsic = R.from_euler(_eulerConvention, x).as_euler('xyz')
    return eulerExtrinsic


def angVel2EulerAngleVel(pqr, eulerZXY):
    """ Angular velocity in the body/flyerish frame to the Euler ZXY angular velocity.
        The reason it's body/flyerish depends on whether you're keeping track of the yaw
        angle in the Euler frame or not, so be careful!
        Stolen from MEAM 620
    """
    phi = eulerZXY[0]
    theta = eulerZXY[1]
    # psi = eulerZXY[2]
    rot_eul2pqr = np.array([[np.cos(theta), 0, -np.cos(phi) * np.sin(theta)],
                            [0, 1, np.sin(phi)],
                            [np.sin(theta), 0, np.cos(phi) * np.cos(theta)]])
    eulVel = np.linalg.inv(rot_eul2pqr).dot(pqr)
    return eulVel


def rotateTensor(tensor, rotationMatrix):
    """ Rotate those tensors!
        Tensors (like inertia tensors) rotate as I_B = R_A^B I_B (R_A^B)^T
        Where
        I_A is the tensor in frame A
        I_B is the tensor in frame B
        R_A^B is the rotation matrix from frame A to frame B

        Note this doesn't do the parallel axis theorem, so still be careful!
    """
    invRot = np.transpose(rotationMatrix)
    firstMultiplication = rotationMatrix.dot(tensor)
    return firstMultiplication.dot(invRot)


def rpm2RadiansPerSecond(x: float) -> float:
    return 0.10472 * x


def rotPerSec2RadiansPerSecond(x: float) -> float:
    return np.pi * 2 * x


def testMathUtil():
    eul = np.array([np.pi / 2, 0, 0])
    print('eul', eul)
    k = [0, 2, 0]
    r = euler2Rotm(eul)
    print(r.dot(k))

    inertia = np.diag([1, 2, 3])

    rot = R.from_euler('z', 45, degrees='true')
    rotm = R.as_dcm(rot)
    print(r.dot(k))
    print(rotm)
    rotateTensor(inertia, rotm)
    yaw = np.pi / 2
    print('eul', eul)
    print("Should be pi, 0, 0", addYaw(eul, yaw))

    yaw_a2b = 90 * np.pi / 180
    rot_a2b = rotmFromYaw(yaw)
    expected90deg = np.array([[0, 1, 0], [-1, 0, 0], [0, 0, 1]])

    if not np.allclose(rot_a2b, expected90deg):
        print('ERR: rotmFromYaw')
        print('Should be 90 degree rotation!')
        print(yaw_a2b)
        print('!= ')
        print(expected90deg)
