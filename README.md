TransformUtilities
==================

*TransformUtilities* is an Objective C class built on Apple's GLKit to facilitate breaking a 3D transformation matrix into a series of simple component transformations.  It is packaged as a class, rather than functions, to allow future configuration of matrix decomposition behavior.

Background
----------

Transformation matrices used in 3D rendering are generally created by combinig a series of basic transformations, falling into the categories of translation, rotation, scaling, shear and perspective.  These may be combined in any sequence, and since matrix multiplication (composition) is not commutative, it is not generally possible to fully reconstruct the original series of operations that went in to creating a given transformation.

It is, however, possible to generate a canonical sequence of basic transformations that will reproduce a given transformation when recombined.  Furthermore, by incorporating certain assumptions, it is possible to create such a series of transformations that reflects standard usage of available APIs, and which may therefore match or come close to the original sequence of transformations.

Decomposition of transformation matrices can be useful whenever the sequence of operations leading to a transformation is not readily accesible and a precise breakdown is not needed, such as in animation, or does not appear to be behaving as expected, where the decomposition can aid in debugging.  Some decomposition is also part of some 3D APIs, such as Apple's Scene Kit, so an open implementation can aid in porting.

Sources
-------

A version of matrix decomposition was apparently provided in Graphics Gems II, and the source can be found at http://tog.acm.org/resources/GraphicsGems/gemsii/unmatrix.c .  However, there are issues with using this code 'as is', not the least of which is that it was written before the adoption of current graphics APIs.

A more recent adaptation of the basic algorithm can be found at http://dev.w3.org/csswg/css-transforms/#matrix-decomposing .  The most notable contribution of this version is the switch from *Euler Angles* to *Quaternions* for representing rotations, which helps avoid some issues related to gimbal lock.

Contents
--------

This repository includes the decomposition class in ./Utility Classes/TUTransformUtilities.h and .m, packaged as an Xcode project using Apple's GLKit template code and incorporating a SenTest unit test suite.  All content is available for use under the MIT license.
