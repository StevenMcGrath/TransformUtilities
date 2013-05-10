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

+ (bool)decomposeTransform:(GLKMatrix4)transform intoScale:(GLKVector3*)scaleOut shearRatios:(GLKVector3*)shearRatiosOut rotation:(GLKVector3*)rotationOut quaternion: (GLKQuaternion*)quaternionOut translation:(GLKVector3*)translationOut perspective:(GLKVector4*)perspectiveOut;
{
 	register int i;
    GLKVector3 scale, shearRatios, rotation, translation;
    GLKVector4 perspective;
 	GLKMatrix4 locmat;
 	GLKMatrix4 pmat, invpmat, tinvpmat;
 	GLKVector4 prhs;
 	GLKVector4 column[4];

 	locmat = transform;
 	
 	if ( locmat.m33 == 0 )
 		return false;
 	/* pmat is used to solve for perspective, but it also provides
 	 * an easy way to test for singularity of the upper 3x3 component.
 	 */
 	pmat = locmat;
 	for ( i=0; i<3; i++ )
 		pmat.m[i*4 + 3] = 0;
 	pmat.m33 = 1;

 	/* First, isolate perspective. */
 	if ( locmat.m03 != 0 || locmat.m13 != 0 ||
 		locmat.m23 != 0 ) {
 		/* prhs is the right hand side of the equation. */
 		prhs.x = locmat.m03;
 		prhs.y = locmat.m13;
 		prhs.z = locmat.m23;
 		prhs.w = locmat.m33;

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
        projection.m03 = perspective.x;
        projection.m13 = perspective.y;
        //
        // Adjust the z scaling to match output of GLKMatrix4MakeFrustrum and GLKMatrix4MakePerspective
        //
        projection.m22 = -1.0/perspective.z;
        projection.m23 = -1.0;
        projection.m32 = -1.0 * perspective.w / perspective.z;
        projection.m33 = 0.0;
        //
        // Invert and multiply to remove the projection from the input transform.
        //
        GLKMatrix4 inverseProjection = GLKMatrix4Invert(projection, &invertable);
        if (!invertable) {
            return false;
        }
        locmat = GLKMatrix4Multiply(inverseProjection, locmat);
        
 		/* Clear the perspective partition. */
        GLKMatrix4SetRow(locmat, 3, (GLKVector4){0, 0, 0, 1});
        
 	} else	{	/* No perspective. */
 		perspective.x = perspective.y = perspective.z = 0;
        perspective.w = 1; // CSS draft sets w to one here.
    }
 	/* Now get scale and shear. */
 	for ( i=0; i<4; i++ ) {
        column[i] = GLKMatrix4GetColumn(locmat, i);
 	}


 	/* Compute X scale factor and normalize first column. */
 	scale.x = GLKVector4Length(column[0]);
 	column[0] = GLKVector4DivideScalar(column[0], scale.x);

 	/* Compute XY shear factor and make 2nd column orthogonal to 1st. */
 	shearRatios.x = GLKVector4DotProduct(column[0], column[1]);
 	column[1] = GLKVector4Add(column[1], GLKVector4MultiplyScalar( column[0], -shearRatios.x));

 	/* Now, compute Y scale and normalize 2nd column. */
 	scale.y = GLKVector4Length(column[1]);
 	column[1] = GLKVector4DivideScalar(column[1], scale.y);
 	shearRatios.x /= scale.y;

 	/* Compute XZ and YZ shears, orthogonalize 3rd column. */
 	shearRatios.y = GLKVector4DotProduct(column[0], column[2]);
 	column[2] = GLKVector4Add(column[2], GLKVector4MultiplyScalar(column[0], -shearRatios.y));
 	shearRatios.z = GLKVector4DotProduct(column[1], column[2]);
 	column[2] = GLKVector4Add(column[2], GLKVector4MultiplyScalar(column[1], -shearRatios.z));

 	/* Next, get Z scale and normalize 3rd column. */
 	scale.z = GLKVector4Length(column[2]);
 	column[2] = GLKVector4DivideScalar(column[2], scale.z);
 	shearRatios.y /= scale.z;
 	shearRatios.z /= scale.z;
    
 	/* Next take care of translation. */
    translation = (GLKVector3){column[3].x, column[3].y, column[3].z};
    column[3].x = column[3].y = column[3].z =0;
        
 	/* At this point, the matrix (in columns[]) is orthonormal.
 	 * Check for a coordinate system flip.  If the determinant
 	 * is -1, then negate the matrix and the scaling factors.
 	 */
 	if ( GLKVector4DotProduct( column[0], GLKVector4CrossProduct( column[1], column[2]) ) < 0 )
 		for ( i = 0; i < 3; i++ ) {
 			scale.v[i] *= -1;
 			column[i].x *= -1;
 			column[i].y *= -1;
 			column[i].z *= -1;
 		}
 
    //
    // Get rotation out as a quaternion
    // What the unmatrix algorithim called rows GLKit calls columns.
    //
    GLKMatrix3 rotationMatrix = GLKMatrix4GetMatrix3(GLKMatrix4MakeWithColumns(column[0], column[1], column[2], column[3]));
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
    if (scaleOut != nil) *scaleOut = scale;
    if (shearRatiosOut != nil) *shearRatiosOut = shearRatios;
    if (rotationOut != nil) *rotationOut = rotation;
    if (quaternionOut != nil) *quaternionOut = quaternion;
    if (translationOut != nil) *translationOut = translation;
    if (perspectiveOut != nil) *perspectiveOut = perspective;
 	return true;
}

@end
