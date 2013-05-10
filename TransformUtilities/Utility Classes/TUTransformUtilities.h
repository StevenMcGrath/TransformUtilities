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
+ (bool)decomposeTransform:(GLKMatrix4)transform intoScale:(GLKVector3*)scale  shearRatios:(GLKVector3*)shearRatios rotation:(GLKVector3*)rotation quaternion:(GLKQuaternion*)quaternion translation:(GLKVector3*)translation perspective:(GLKVector4*)perspective;

@end
