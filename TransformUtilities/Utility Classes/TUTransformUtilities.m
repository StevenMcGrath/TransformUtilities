//
//  TransformUtilities.m
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
 * SPMathUtils.m is adapted from unmatrix.c:
 *
 * unmatrix.c - given a 4x4 matrix, decompose it into standard operations.
 *
 * Author:	Spencer W. Thomas
 * 		University of Michigan
 *
 * http://users.soe.ucsc.edu/~pang/160/f98/Gems/GemsII/unmatrix.c
 *
 */
#include <math.h>
#include "TUTransformUtilities.h"

/* decomposeTransform is based on unmatrix:
 *
 * unmatrix - Decompose a non-degenerate 4x4 transformation matrix into
 * 	the sequence of transformations that produced it.
 * [Sx][Sy][Sz][Shearx/y][Sx/z][Sz/y][Rx][Ry][Rz][Tx][Ty][Tz][P(x,y,z,w)]
 *
 * Returns true upon success, false if the matrix is singular.
 */

@implementation TUTransformUtilities

+ (bool)getRawPerspective:(GLKVector4*)perspectiveOut fromTransform:(GLKMatrix4)transformIn;
{
    GLKMatrix4 pmat;
    int i;
    GLKMatrix4 tinvpmat;
    GLKMatrix4 invpmat;
    GLKVector4 perspective;
 	GLKVector4 prhs;
 	/* pmat is used to solve for perspective, but it also provides
 	 * an easy way to test for singularity of the upper 3x3 component.
 	 */
 	pmat = transformIn;
 	for ( i=0; i<3; i++ )
 		pmat.m[i*4 + 3] = 0;
 	pmat.m33 = 1;
    
 	/* First, isolate perspective. */
 	if ( transformIn.m03 != 0 || transformIn.m13 != 0 ||
 		transformIn.m23 != 0 ) {
 		/* prhs is the right hand side of the equation. */
 		prhs.x = transformIn.m03;
 		prhs.y = transformIn.m13;
 		prhs.z = transformIn.m23;
 		prhs.w = transformIn.m33;
        
 		/* Solve the equation by inverting pmat and multiplying
 		 * prhs by the inverse.  (This is the easiest way, not
 		 * necessarily the best.)
 		 */
        bool invertable;
        invpmat = GLKMatrix4Invert(pmat, &invertable);
        if (!invertable) {
            return false;
        }
		tinvpmat = GLKMatrix4Transpose (invpmat);
 		perspective = GLKMatrix4MultiplyVector4(tinvpmat, prhs);
 	} else	{	/* No perspective. */
 		perspective.x = perspective.y = perspective.z = 0;
        perspective.w = 1; // CSS draft sets w to one here.
    }
    if (perspectiveOut != nil) *perspectiveOut = perspective;
    return true;
}

+ (bool)decomposeProjectionTransform:(GLKMatrix4)transformIn intoPerspective:(GLKVector4*)perspectiveOut residualTransform:(GLKMatrix4*)residualOut;
{
    bool result;
    bool invertable;
    GLKVector4 rawPerspective;
    GLKVector3 scale, shear;
    float fieldOfView, aspectRatio;
    result = [self getRawPerspective:&rawPerspective fromTransform:transformIn];
    if (GLKVector4AllEqualToVector4(rawPerspective, (GLKVector4){0,0,0,1})) {
        // No perspective.
        if (perspectiveOut != nil) {
            *perspectiveOut = (GLKVector4){0,1,0,0};
        }
        if (residualOut != nil) {
            *residualOut = transformIn;
        }
    }
    else {
    //    NSLog(@"Raw perspective: %@", NSStringFromGLKVector4(rawPerspective));

        //
        // Recreate the projection matrix.  This takes the output from the original unmatrix approach to isolating
        // perspective and adjusts it to a representation appropriate to Open GL and GLKit, as is produced when
        // calling GLKMatrixMakePerspective.  For contrast, compare the CSS Transforms Level 1 draft at
        // http://dev.w3.org/csswg/css-transforms/ .  Open GL coordinates have the positive y axis pointing up
        // and the positive z axis pointing away from the viewer, and perspective specified in terms of field of view,
        // while CSS coordinates have the positive y axis pointing down, the positive z axis pointing towards the
        // viewer, and perspective specified in terms of the viewer's presumed distance from the image.
        //
        GLKMatrix4 projection = GLKMatrix4Identity;
        //
        // Convert to inputs to GLKMatrix4MakePerspective
        //
        float nearz, farz, q;
        q = (1 - rawPerspective.z)/(1 + rawPerspective.z);
        nearz = rawPerspective.w *(1 + q)/2.0;
        farz = nearz / q;
        //        NSLog(@"Near z = %f, Far z = %f", nearz, farz);
        {
            GLKMatrix4 projection = GLKMatrix4Identity;
            projection.m03 = rawPerspective.x;
            projection.m13 = rawPerspective.y;
            //
            // Adjust the z scaling to match output of GLKMatrix4MakeFrustrum and GLKMatrix4MakePerspective
            //
            projection.m22 = -1.0/rawPerspective.z;
            projection.m23 = -1.0;
            projection.m32 = -1.0 * rawPerspective.w / rawPerspective.z;
            projection.m33 = 0.0;
            //
            // Invert and multiply to remove the projection from the input transform.
            //
            GLKMatrix4 inverseProjection = GLKMatrix4Invert(projection, &invertable);
            if (!invertable) {
                return false;
            }
            GLKMatrix3 deprojected = GLKMatrix4GetMatrix3(GLKMatrix4Multiply(inverseProjection, transformIn));
            //
            // Get scale and shear, to build the frustum
            //
            GLKMatrix3 residual;
            //
            // Using the transpose allows the same method to be used to factor the scale out from the left (perspective)
            // side of the transformation, rather than the right side, as for canonical model decomposition.
            //
            [self factorTransform:GLKMatrix3Transpose(deprojected) intoScale:&scale shearRatios:&shear residualTransform:&residual];
            //
            // Since we are extracting perspective, and not an uneven frustum, we want to factor out the shear from the
            // right (model side), so we can get the scale out for the perspective,
            //
            // The formula to be factored is Scale * Rotation * Shear.
            //
            // TODO: Develop a formula for decomposing a linear (3 x 3) transformation i the order Scale * Rotation * Shear.
            //
/*
                The following attempts do not work.
 
//            GLKMatrix3 inverseShear = GLKMatrix3Identity;
//            inverseShear.m10 = -shear.v[0];
//            inverseShear.m21 = -shear.v[2];
//            inverseShear.m20 = -shear.v[1] + (shear.v[0] * shear.v[2]);
//            
//            deprojected = GLKMatrix3Multiply(deprojected, inverseShear);
//            
//            [self factorTransform:GLKMatrix3Transpose(deprojected) intoScale:&scale shearRatios:&shear residualTransform:&residual];
            //
            // Since the scale is to be applied to the prespective transform, it need to be transformed to the
            // appropriate space.
            //
            // The residual matrix is a rotation, so its transpose is its inverse.
            //
//            residual = GLKMatrix3Multiply(deprojected, GLKMatrix3Transpose(residual));
//            scale = (GLKVector3){residual.m00, residual.m11, residual.m22};
//            [self factorTransform:residual intoScale:&scale shearRatios:&shear residualTransform:&residual];
//            scale = GLKMatrix3MultiplyVector3(GLKMatrix3Transpose(residual), scale);
            //
*/
            fieldOfView = fabs(2.0 * atanf(1.0/scale.y));
            aspectRatio = fabs(scale.y / scale.x);
        }
        if (perspectiveOut != nil) {
            perspectiveOut->x = fieldOfView;
            perspectiveOut->y = aspectRatio;
            perspectiveOut->z = nearz;
            perspectiveOut->w = farz;
        }
        if (residualOut != nil) {
            projection = GLKMatrix4MakePerspective (fieldOfView, aspectRatio, nearz, farz);
            GLKMatrix4 inverseProjection = GLKMatrix4Invert(projection, &invertable);
            if (!invertable) {
                return false;
            }
            *residualOut = GLKMatrix4Multiply(inverseProjection, transformIn);
        }
    }
    return true;
}

+ (bool)decomposeProjectionTransform:(GLKMatrix4)transformIn intoFrustum:(GLKVector3*)lowerBounds :(GLKVector3*)upperBounds residualTransform:(GLKMatrix4*)residualOut;
{
    bool result;
    GLKVector4 rawPerspective;
    GLKVector3 scale, shear;
    result = [self getRawPerspective:&rawPerspective fromTransform:transformIn];
    //    NSLog(@"Raw perspective: %@", NSStringFromGLKVector4(rawPerspective));
    
    //
    // Recreate the projection matrix.  This takes the output from the original unmatrix approach to isolating
    // perspective and adjusts it to a representation appropriate to Open GL and GLKit, as is produced when
    // calling GLKMatrixMakePerspective.  For contrast, compare the CSS Transforms Level 1 draft at
    // http://dev.w3.org/csswg/css-transforms/ .  Open GL coordinates have the positive y axis pointing up
    // and the positive z axis pointing away from the viewer, and perspective specified in terms of field of view,
    // while CSS coordinates have the positive y axis pointing down, the positive z axis pointing towards the
    // viewer, and perspective specified in terms of the viewer's presumed distance from the image.
    //
    GLKMatrix4 projection = GLKMatrix4Identity;
    bool invertable;
    //
    // Convert to inputs to GLKMatrix4MakeFrustum
    //
    float nearz, farz, q, left, right, top, bottom;
    q = (1 - rawPerspective.z)/(1 + rawPerspective.z);
    nearz = rawPerspective.w *(1 + q)/2.0;
    farz = nearz / q;
    //        NSLog(@"Near z = %f, Far z = %f", nearz, farz);
    {
        GLKMatrix4 projection = GLKMatrix4Identity;
        projection.m03 = rawPerspective.x;
        projection.m13 = rawPerspective.y;
        //
        // Adjust the z scaling to match output of GLKMatrix4MakeFrustrum and GLKMatrix4MakePerspective
        //
        projection.m22 = -1.0/rawPerspective.z;
        projection.m23 = -1.0;
        projection.m32 = -1.0 * rawPerspective.w / rawPerspective.z;
        projection.m33 = 0.0;
        //
        // Invert and multiply to remove the projection from the input transform.
        //
        GLKMatrix4 inverseProjection = GLKMatrix4Invert(projection, &invertable);
        if (!invertable) {
            return false;
        }
        GLKMatrix3 deprojected = GLKMatrix4GetMatrix3(GLKMatrix4Multiply(inverseProjection, transformIn));
        //
        // Get scale and shear, to build the frustum
        //
        GLKMatrix3 residualTransform;
        [self factorTransform:GLKMatrix3Invert(deprojected, &invertable) intoScale:&scale shearRatios:&shear residualTransform:nil];
        GLKMatrix3 inverseScaleMatrix = {scale.x, 0.0, 0.0, 0.0, scale.y, 0.0, 0.0, 0.0, scale.z};
        scale.x = 1.0/scale.x;
        scale.y = 1.0/scale.y;
        scale.z = 1.0/scale.z;
        //
        // This could be cleaner. Simplification needed.
        //
        residualTransform = GLKMatrix3Multiply(inverseScaleMatrix, deprojected);
        [self factorTransform:GLKMatrix3Invert(residualTransform, &invertable) intoScale:nil shearRatios:&shear residualTransform:nil];
        shear.v[0] = -shear.v[0]*scale.x/scale.y;
        shear.v[2] = -shear.v[2]*scale.y/scale.z;
        shear.v[1] = -shear.v[1]*scale.x/scale.z+(shear.v[0]*shear.v[1]);

        float yspan = 2.0 * nearz / scale.y;
        float xspan = 2.0 * nearz /scale.x;
        right = xspan * (shear.v[1]+1.0)/2.0;
        left = xspan * (shear.v[1]-1.0)/2.0;
        top = yspan * (shear.v[2]+1.0)/2.0;
        bottom = yspan * (shear.v[2]-1.0)/2.0;
        //        NSLog (@"Field of View: %f radians (%f degrees), aspect ratio: %f.", fieldOfView, GLKMathRadiansToDegrees(fieldOfView),aspectRatio);
        //        NSLog (@"Frustum bounds: {{ %f, %f, %f } { %f, %f, %f }}", left, bottom, nearz, right, top, farz);
        // NSLog (@"Intermediate transform: %@",NSStringFromGLKMatrix3(residualTransform));
    }
    if ((lowerBounds != nil) && (upperBounds != nil)) {
        lowerBounds->x = left;
        upperBounds->x = right;
        lowerBounds->y = bottom;
        upperBounds->y = top;
        lowerBounds->z = nearz;
        upperBounds->z = farz;
    }
    if (residualOut != nil) {
        projection = GLKMatrix4MakeFrustum (left, right, bottom, top, nearz, farz);
        // NSLog (@"Frustum transform: %@",NSStringFromGLKMatrix4(projection));
        GLKMatrix4 inverseProjection = GLKMatrix4Invert(projection, &invertable);
        if (!invertable) {
            return false;
        }
        *residualOut = GLKMatrix4Multiply(inverseProjection, transformIn);
    }
    return true;
}

+ (void)factorTransform:(GLKMatrix3)transformIn intoScale:(GLKVector3*)scaleOut shearRatios:(GLKVector3*)shearRatiosOut residualTransform:(GLKMatrix3*)residualOut;
{
 	/* Now get scale and shear. */
    register int i;
    GLKVector3 scale, shearRatios;
 	GLKVector3 column[3];
    for ( i=0; i<3; i++ ) {
        column[i] = GLKMatrix3GetColumn(transformIn, i);
 	}
    
    
 	/* Compute X scale factor and normalize first column. */
 	scale.x = GLKVector3Length(column[0]);
 	column[0] = GLKVector3DivideScalar(column[0], scale.x);
    
 	/* Compute XY shear factor and make 2nd column orthogonal to 1st. */
 	shearRatios.v[0] = GLKVector3DotProduct(column[0], column[1]);
 	column[1] = GLKVector3Add(column[1], GLKVector3MultiplyScalar( column[0], -shearRatios.v[0]));
    
 	/* Now, compute Y scale and normalize 2nd column. */
 	scale.y = GLKVector3Length(column[1]);
 	column[1] = GLKVector3DivideScalar(column[1], scale.y);
 	shearRatios.v[0] /= scale.y;
    
 	/* Compute XZ and YZ shears, orthogonalize 3rd column. */
 	shearRatios.v[1] = GLKVector3DotProduct(column[0], column[2]);
 	column[2] = GLKVector3Add(column[2], GLKVector3MultiplyScalar(column[0], -shearRatios.v[1]));
 	shearRatios.v[2] = GLKVector3DotProduct(column[1], column[2]);
 	column[2] = GLKVector3Add(column[2], GLKVector3MultiplyScalar(column[1], -shearRatios.v[2]));
    
 	/* Next, get Z scale and normalize 3rd column. */
 	scale.z = GLKVector3Length(column[2]);
 	column[2] = GLKVector3DivideScalar(column[2], scale.z);
 	shearRatios.v[1] /= scale.z;
 	shearRatios.v[2] /= scale.z;
    
    // Now remove the shear and get the unsheared scale factors
//    
//    GLKMatrix3 inverseShear = GLKMatrix3Identity;
//    inverseShear.m10 = -shearRatios.v[0];
//    inverseShear.m21 = -shearRatios.v[2];
//    inverseShear.m20 = -shearRatios.v[1] + (shearRatios.v[0] * shearRatios.v[2]);
//    
//    transformIn = GLKMatrix3Multiply(inverseShear, transformIn);
//    
//    for ( i=0; i<3; i++ ) {
//        column[i] = GLKMatrix3GetRow(transformIn, i);
// 	}
// 	scale.x = GLKVector3Length(column[0]);
// 	column[0] = GLKVector3DivideScalar(column[0], scale.x);
// 	scale.y = GLKVector3Length(column[1]);
// 	column[1] = GLKVector3DivideScalar(column[1], scale.y);
// 	scale.z = GLKVector3Length(column[2]);
// 	column[2] = GLKVector3DivideScalar(column[2], scale.z);
//
    
    if (scaleOut != nil) *scaleOut = scale;
    if (shearRatiosOut != nil) *shearRatiosOut = shearRatios;
    if (residualOut != nil) *residualOut = GLKMatrix3MakeWithColumns(column[0], column[1], column[2]);
}

+ (bool)decomposeModelViewTransform:(GLKMatrix4)transformIn intoScale:(GLKVector3*)scaleOut shearRatios:(GLKVector3*)shearRatiosOut rotation:(GLKVector3*)rotationOut quaternion: (GLKQuaternion*)quaternionOut translation:(GLKVector3*)translationOut
{
    register int i;
 	GLKVector3 column[3];
    GLKVector3 rotation, translation;
    
 	/* Next take care of translation. */
    translation = (GLKVector3){transformIn.m30, transformIn.m31, transformIn.m32};
    transformIn.m30 = transformIn.m31 = transformIn.m32 =0;
    
    GLKMatrix3 transform3In = GLKMatrix4GetMatrix3(transformIn);
 	/* Now get scale and shear. */
    [self factorTransform:transform3In intoScale:scaleOut shearRatios:shearRatiosOut residualTransform:&transform3In];
    
    for ( i=0; i<3; i++ ) {
        column[i] = GLKMatrix3GetColumn(transform3In, i);
 	}
    
 	/* At this point, the matrix (in columns[]) is orthonormal.
 	 * Check for a coordinate system flip.  If the determinant
 	 * is -1, then negate the matrix and the scaling factors.
 	 */
 	if ( GLKVector3DotProduct( column[0], GLKVector3CrossProduct( column[1], column[2]) ) < 0 )
 		for ( i = 0; i < 3; i++ ) {
            if (scaleOut != nil) {
                scaleOut->v[i] *= -1;
            }
 			column[i].x *= -1;
 			column[i].y *= -1;
 			column[i].z *= -1;
 		}
    
    //
    // Get rotation out as a quaternion
    // What the unmatrix algorithim called rows GLKit calls columns.
    //
    GLKMatrix3 rotationMatrix = GLKMatrix3MakeWithColumns(column[0], column[1], column[2]);
    GLKQuaternion quaternion = GLKQuaternionMakeWithMatrix3(rotationMatrix);
    
 	/* Now, get the rotations out, as described in the gem. */
 	rotation.y = asin(-column[0].z);
    //
    // The following test is designed to address gimbal lock conditions using Euler angles. The test interval should
    // be selected to minimize numerical instability, but is only approximated here.  Note that the value will
    // depend on the floating point precision used, which in this case (GLKit) is single precision.
    //
 	if ( fabsf(cosf(rotation.y)) > 0.001 ) {
 		rotation.x = atan2(column[1].z, column[2].z);
 		rotation.z = atan2(column[0].y, column[0].x);
 	} else {
 		rotation.x = atan2(-column[2].y, column[1].y);
 		rotation.z = 0;
 	}
 	/* All done! */
    if (rotationOut != nil) *rotationOut = rotation;
    if (quaternionOut != nil) *quaternionOut = quaternion;
    if (translationOut != nil) *translationOut = translation;
    return true;
}

+ (bool)decomposeTransform:(GLKMatrix4)transformIn intoScale:(GLKVector3*)scaleOut shearRatios:(GLKVector3*)shearRatiosOut rotation:(GLKVector3*)rotationOut quaternion: (GLKQuaternion*)quaternionOut translation:(GLKVector3*)translationOut perspective:(GLKVector4*)perspectiveOut;
{
    GLKVector4 perspective;
 	bool result;
    
    result = [self decomposeProjectionTransform:transformIn intoPerspective:&perspective residualTransform:&transformIn];
    
    if (result) {
        if (perspectiveOut != nil) *perspectiveOut = perspective;    
        result = [self decomposeModelViewTransform:transformIn intoScale:scaleOut shearRatios:shearRatiosOut rotation:rotationOut quaternion:quaternionOut translation:translationOut];
    }
 	return result;
}

+ (bool)decomposeTransform:(GLKMatrix4)transformIn intoScale:(GLKVector3*)scaleOut shearRatios:(GLKVector3*)shearRatiosOut rotation:(GLKVector3*)rotationOut quaternion: (GLKQuaternion*)quaternionOut translation:(GLKVector3*)translationOut frustum:(GLKVector3*)lowerBoundsOut :(GLKVector3*)upperBoundsOut;
{
    GLKVector3 lowerBounds, upperBounds;
 	bool result;
    
    result = [self decomposeProjectionTransform:transformIn intoFrustum:&lowerBounds :&upperBounds residualTransform:&transformIn];
    
    if (result) {
        if ((lowerBoundsOut != nil) && (upperBoundsOut != nil)) {
            *lowerBoundsOut = lowerBounds;
            *upperBoundsOut = upperBounds;
        }
        result = [self decomposeModelViewTransform:transformIn intoScale:scaleOut shearRatios:shearRatiosOut rotation:rotationOut quaternion:quaternionOut translation:translationOut];
    }
 	return result;
}

@end
