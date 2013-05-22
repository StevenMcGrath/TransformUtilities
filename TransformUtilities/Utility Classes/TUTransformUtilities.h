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
+ (bool)decomposeTransform:(GLKMatrix4)transformIn intoScale:(GLKVector3*)scaleOut  shearRatios:(GLKVector3*)shearRatiosOut rotation:(GLKVector3*)rotationOut quaternion:(GLKQuaternion*)quaternionOut translation:(GLKVector3*)translationOut perspective:(GLKVector4*)perspectiveOut;

+ (bool)decomposeTransform:(GLKMatrix4)transformIn intoScale:(GLKVector3*)scaleOut shearRatios:(GLKVector3*)shearRatiosOut rotation:(GLKVector3*)rotationOut quaternion: (GLKQuaternion*)quaternionOut translation:(GLKVector3*)translationOut frustum:(GLKVector3*)lowerBoundsOut :(GLKVector3*)upperBoundsOut;

@end
