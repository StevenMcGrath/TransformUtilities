//
//  TransformUtilities.h
//  TransformUtilities
//
//  Created by Steven McGrath on 5/8/13 and made available under the MIT license.
/*
 The MIT License (MIT)
 
 Copyright (c) 2013 Steven McGrath
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
 */

/*
 * SPMathUtils.h is adapted from unmatrix.h:
 *
 * unmatrix.h - Definitions for using unmatrix
 *
 * Author:	Spencer W. Thomas
 *		University of Michigan
 *
 * http://users.soe.ucsc.edu/~pang/160/f98/Gems/GemsII/unmatrix.h
 *
 */

#import <GLKit/GLKit.h>

/* The unmatrix subroutine fills in a vector of floating point
 * values.  These symbols make it easier to get the data back out.
 */


@interface TUTransformUtilities : NSObject

/* decomposeTransform is based on unmatrix:
 *
 * unmatrix - Decompose a non-degenerate 4x4 transformation matrix into
 * 	the sequence of transformations that produced it.
 * [Sx][Sy][Sz][Shearx/y][Sx/z][Sz/y][Rx][Ry][Rz][Tx][Ty][Tz][P(x,y,z,w)]
 *
 * The coefficient of each transformation is returned in the corresponding
 * element of the vector arguments.
 *
 * Returns true upon success, false if the matrix is singular.
 */

// Comprehensive decompositions:

//
// Decompose a combined model/view/projection transformation matrix into components such that transformation = projection *
// translation * rotation * shear * scale.  This routine returns the projection component as a vector of four values
// representing field of view, aspect ratio, near z and far z, which are the inputs for GLKMatrix4MakePerspective. Rotation
// can be retrieved as either a three component vector of Euler angles, or as a quaternion. The quaternion output is
// prefered as it does not introduce problems with gimbal lock, and can easily be converted to axis / angle form if desired.
//
// Note that if a projection element of the transformation is found, any scale factors from the input may be incorporated
// in the perspective calculation, and the scale element may be empty, even when scaling as part of the original model/view
// portion of the input transform.  A further limitation at present is that any shear present in the transform can interfere
// with returning results that match to original composition of component transforms when a projection is present.
//
+ (bool)decomposeTransform:(GLKMatrix4)transformIn perspective:(GLKVector4*)perspectiveOut translation:(GLKVector3*)translationOut rotation:(GLKVector3*)rotationOut quaternion: (GLKQuaternion*)quaternionOut shearRatios:(GLKVector3*)shearRatiosOut scale:(GLKVector3*)scaleOut;

//
// Decompose a combined model/view/projection transformation matrix into components such that transformation = projection *
// translation * rotation * shear * scale.  This routine returns the projection component as two vectors of three values
// each, representing the lower and upper bounds of a frustum, which are the inputs for GLKMatrix4MakeFrustum. Rotation
// can be retrieved as either a three component vector of Euler angles, or as a quaternion. The quaternion output is
// prefered as it does not introduce problems with gimbal lock, and can easily be converted to axis / angle form if desired.
//
// Note that if a projection element of the transformation is found, any scale and shear factors from the input may be
// incorporated in the frustum calculation, and the scale and shear elements may be empty, even when scaling or shearing
// are part of the original model/view portion of the input transform.
//
+ (bool)decomposeTransform:(GLKMatrix4)transformIn frustum:(GLKVector3*)lowerBoundsOut :(GLKVector3*)upperBoundsOut translation:(GLKVector3*)translationOut rotation:(GLKVector3*)rotationOut quaternion:(GLKQuaternion*)quaternionOut shearRatios:(GLKVector3*)shearRatiosOut scale:(GLKVector3*)scaleOut;

// Projection Decompositions:

//
// Isolates the elements that comprise a projection matrix as produced by GLKMatrix4MakePerspective, returning those elements,
// field of view, aspect ratio, near z and far z, as the four elements of perspectiveOut.  Any remaining transformations from
// the input are returned in the residualTransform matrix, such that transformIn = projection * residualTransform.  This
// residual comprises a model/view transform.
//
// Limitations: Since part of what makes up a projection matrix is scaling, it is not possible for this routine to distinguish
// scaling that may have been applied to the model from that which was part of the original projection transform.  It will
// therefore only return a perspective output that matches the original input when the model/view transformation does not
// include scaling.
//
// At this time, shear is also not effectively isolated, and results will vary in the presence of shear in the model/view transform.
//
+ (bool)decomposeProjectionTransform:(GLKMatrix4)transformIn intoPerspective:(GLKVector4*)perspectiveOut residualTransform:(GLKMatrix4*)residualOut;

//
// Isolates the elements that comprise a projection matrix as produced by GLKMatrix4MakeFrustum, returning those elements,
// the left, bottom, near z, right top, and far z clip limits, in the vectors lowerBounds and UpperBounds.  Any remaining
// transformations from the input are returned in the residualTransform matrix, such that
// transformIn = projection * residualTransform.  This residual comprises a model/view transform.
//
// Limitations: Since a projection frustum matrix imparts both scaling and shear, it is not possible for this routine to
// distinguish scaling or shearing that may have been applied to the model from that which was part of the original
// projection transform.  It will therefore only return a frustum output that matches the original input when the
// model/view transformation does not include scaling or shear.
//
+ (bool)decomposeProjectionTransform:(GLKMatrix4)transformIn intoFrustum:(GLKVector3*)lowerBounds :(GLKVector3*)upperBounds residualTransform:(GLKMatrix4*)residualOut;

// Model / View Decomposition:

//
// Decompose a model/view transformation matrix into component transformations such that transformation = translation *
// rotation * shear * scale.  Rotation can be retrieved as either a three component vector of Euler angles, or as a
// quaternion.  The quaternion output is prefered as it does not introduce problems with gimbal lock, and can easily be
// converted to axis / angle for if desired.
//
+ (bool)decomposeModelViewTransform:(GLKMatrix4)transformIn intoTranslation:(GLKVector3*)translationOut rotation:(GLKVector3*)rotationOut quaternion: (GLKQuaternion*)quaternionOut shearRatios:(GLKVector3*)shearRatiosOut scale:(GLKVector3*)scaleOut;

@end
