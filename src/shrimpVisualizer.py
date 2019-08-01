#!/usr/bin/env python
"""
3D visualizer for Shrimp Project
Rebecca Li + Allen Yu 2019

All units are SI units, using meters.

A Note on VPython Axes:
X is to the right of the screen
Y is towards the top of the screen
Z is out of the screen
"""

import numpy as np
import vpython as vp

from mathUtil import addYaw, eulerExtXYZfromEulerShrimp
from shrimpClasses import defaultShrimpParams


class SceneParameters:
    """ Parameters of the Vpython Shrimp Visualization
    """

    def __init__(self, groundWidth, groundLength, groundHeight, canvasHeight, canvasWidth):
        self.groundWidth = groundWidth
        self.groundLength = groundLength
        self.groundHeight = groundHeight
        self.canvasWidth = canvasWidth
        self.canvasHeight = canvasHeight


def vpVecFromArr(x):
    """ Cast a 3x1 array to Vp vector """
    assert np.shape(x) == (3,), "Error: trying to cast array that is not 3x1 to a Vpython Vector"
    return vp.vector(x[0], x[1], x[2])


def odeStateToOrientationVecs(odeState):
    """ Gets the position and orientation from ODE vector
        Returns:
            r_w2b_w (vp.vector)
            euler_w2b_xyz (vp.vector) the orientation  of the body in the Euler Extrinsic XYZ order
    """
    r_w2b_w = odeState[0:3]
    # In the VPython frame, rotation around x-axis
    rot_w2v = np.array([[1, 0, 0], [0, 0, -1], [0, -1, 0]])
    r_w2b_v = rot_w2v.dot(r_w2b_w)
    euler_w2f = odeState[6:9]
    yaw_f2b = odeState[14]
    # yaw_b2p = odeState[15]
    euler_w2b = addYaw(euler_w2f, yaw_f2b)
    euler_w2b_w_xyz = eulerExtXYZfromEulerShrimp(euler_w2b)
    # Reorder euler angles around the correct axes
    euler_w2b_v_xyz = rot_w2v.dot(euler_w2b_w_xyz)
    return vpVecFromArr(r_w2b_v), vpVecFromArr(euler_w2b_v_xyz)


def defaultSceneParameters():
    groundWidth = 1e-3
    groundLength = 1e-3
    groundHeight = 0.0001
    canvasHeight = 400
    canvasWidth = 600
    return SceneParameters(groundWidth, groundLength, groundHeight, canvasHeight, canvasWidth)


def drawGround(p):
    """ Generates default ground
        Arguments:
            p (SceneParameters)
    """
    _ = vp.box(pos=vp.vector(0, -p.groundHeight, 0), length=p.groundLength, width=p.groundWidth, height=p.groundHeight,
               color=vp.color.cyan)


def rotateShrimpBody(shrimpBody, oldEuler_w2b, newEuler_w2b):
    """ Rotates the body about the world axis centered at the body
    """
    deltaEuler = oldEuler_w2b - newEuler_w2b
    # The axes for the world are different than our
    shrimpBody.rotate(angle=deltaEuler.x, axis=vp.vector(1, 0, 0))
    shrimpBody.rotate(angle=deltaEuler.y, axis=vp.vector(0, 1, 0))
    shrimpBody.rotate(angle=deltaEuler.z, axis=vp.vector(0, 0, 1))


def drawShrimpBody(shrimpParams, r_w2b_w, euler_w2b):
    """ Generates the Shrimp object to draw, so the vp.cylinder and propeller
        Arguments:
            shrimpParams (ShrimpParameters)
            r_w2b_w (vp.vector)
            euler_w2b (vp.vector)
    """

    # Assume Body position is the center of mass, or center of b frame
    p = shrimpParams.vizParams
    height_b2p = shrimpParams.shaftPropParams.height_b2p

    shrimpStator = vp.box(pos=r_w2b_w, length=p.bodyWidth, width=p.bodyWidth, height=p.bodyHeight,
                          color=vp.color.orange)
    propPosOffset = vp.vector(0, height_b2p, 0)
    propPos = propPosOffset + r_w2b_w
    shrimpProp = vp.cylinder(pos=propPos, axis=vp.vector(0, p.propDiscThickness, 0),
                             radius=p.propDiscRadius, color=vp.color.magenta)
    shrimpBody = vp.compound([shrimpStator, shrimpProp])

    rotateShrimpBody(shrimpBody, vp.vector(0, 0, 0), euler_w2b)
    return shrimpBody


def setupShrimpScene(vizParams, shrimpParams, r_w2b_w, euler_w2b, autoplay):
    """ Set up scene with Shrimp. Returns Shrimp Drawing object.
        Arguments:
            vizParams (SceneParameters)
            shrimpParams (ShrimpParameters)
            r_w2b_w (vp.vector)
            euler_w2b (vp.vector)
    """
    scene = vp.canvas(width=vizParams.canvasWidth, height=vizParams.canvasHeight,
                      title='Shrimp World. Press any key to begin')
    # scene.autoscale = False
    drawGround(vizParams)  # draws ground
    shrimpBody = drawShrimpBody(shrimpParams, r_w2b_w, euler_w2b)
    sceneText = initializeShrimpText(shrimpBody.pos)
    if not autoplay:
        scene.waitfor('keydown')
    scene.camera.follow(shrimpBody)
    scene.autoscale = True
    return scene, shrimpBody, sceneText


def initializeShrimpText(shrimpPosition):
    """ Initializes Text in the scene dialog
        Arguments:
            shrimpPosition: vector
    """
    L = vp.label(pos=shrimpPosition,
                 text=shrimpPosition, xoffset=20,
                 yoffset=50, space=30,
                 height=16, border=4,
                 font='sans')
    return L


def updateText(textObject, shrimpPosition):
    """ Updates text in the scene dialog
        Arguments:
             textObject: TODO
             shrimpPosition : vector
    """
    textObject.text = shrimpPosition
    textObject.pos = shrimpPosition


def updateShrimpBody(shrimpBody, oldEuler_w2b, newOdeState):
    """ Updates shrimp position and orientation
    """
    (pos, newEuler_w2b) = odeStateToOrientationVecs(newOdeState)
    shrimpBody.pos = pos
    rotateShrimpBody(shrimpBody, oldEuler_w2b, newEuler_w2b)
    return (newEuler_w2b)


def drawShrimp(shrimpParams, timeStamps, odeStates, autoplay=False):
    """ Draws and animates 3d shrimp from odestates and timestamp
        Arguments:
            shrimpParams (ShrimpParameters)
            timeStamps (np.ndarray): nx1 array of times [s]
            odeStates (np.ndarray): nx16 array of ODE states, fully described in shrimpClasses
    """
    vizParams = defaultSceneParameters()

    initialOdeState = odeStates[0]
    odeStatesTail = odeStates[1:]

    (initial_r_w2b_w, initial_euler_w2b) = odeStateToOrientationVecs(initialOdeState)
    (scene, shrimpBody, shrimpText) = setupShrimpScene(vizParams, shrimpParams, initial_r_w2b_w,
                                                       initial_euler_w2b, autoplay)

    oldEuler_w2b = initial_euler_w2b
    # dt = timeStamps[1] - timeStamps[0]  # calculate what dt is
    # TODO: maybe calculate dt based off of actual time, or have a slowdown factor
    dt = 0.05  # 20 frames per second
    timeStampsTail = timeStamps[1:]
    positions = [shrimpBody.pos]
    # For each timestamp, we're going to draw the body
    print('Animating Shrimp...')
    for (timeStamp, newOdeState) in zip(timeStampsTail, odeStatesTail):  # animation loop here

        # update body position and orientation
        oldEuler_w2b = updateShrimpBody(shrimpBody, oldEuler_w2b, newOdeState)
        # call update function, pass in the text object and position of shrimp
        updateText(shrimpText, shrimpBody.pos)
        scene.title = 'time: %f seconds' % timeStamp
        # sleep for appropriate amount of time
        positions.append(shrimpBody.pos)
        vp.sleep(dt)
    scene.title = 'Simulation complete. Time: %f seconds' % timeStamps[-1]
    print('Visualizer finished. You may close the program')
    return scene


def testShrimpVisualizer():
    print('Testing shrimpVisualizer')
    p = defaultShrimpParams()
    n = 30
    dt = 0.1
    t = np.linspace(0, n * dt, n)
    deltaPos = np.linspace(0, 1, n)
    stateSize = 16
    states = np.zeros((n, stateSize))
    states[:, 8] = deltaPos
    scene = drawShrimp(p, t, states, autoplay=True)
    scene.delete()
