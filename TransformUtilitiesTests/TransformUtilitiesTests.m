//
//  TransformUtilitiesTests.m
//  TransformUtilitiesTests
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

#import "TransformUtilitiesTests.h"
#import <GLKit/GLKit.h>
#import "TUTransformUtilities.h"

@implementation TransformUtilitiesTests {
    // Declare inputs.
    GLKVector3 scaleIn, shearRatiosIn, rotationAxisIn, translationIn;
    GLKVector3 eyeIn, centerIn, upIn;
    GLKVector4 frustrumBoundsIn;
    float rotationAngleIn, fieldOfViewIn, aspectRatioIn, nearClipIn, farClipIn;
    GLKMatrix4 transformIn; // Computed from other input componenets as input for decomposition
    GLKMatrix4 projectionIn, projectionOut;
    // Declare vectors to receive results.
    GLKVector3 scaleOut, shearRatiosOut, rotationOut, translationOut;
    GLKQuaternion quaternionIn, quaternionOut;
    GLKVector4 perspectiveOut;
    GLKMatrix4 transformOut;    // Reconstructed from outputs for comparison to input.
    bool bypassScaleCheck, bypassTranslationCheck, bypassShearCheck;
    float expectedAccuracy, maxError;
}

static const GLKVector3 invalidGLKVector3 = {-1.0,-1.0,-1.0};
static const GLKVector3 zeroGLKVector3 = {0.0,0.0,0.0};
static const GLKVector3 unitGLKVector3 = {1.0,1.0,1.0};
static const GLKVector4 invalidGLKVector4 = {-1.0,-1.0,-1.0,-1.0};
static const GLKQuaternion invalidGLKQuaternion = {0.0, 0.0, 0.0, 0.0};
static const float coverage = 0.05; // Fraction of possible tests to run in long running cases.
static const float selectionThreshold = coverage * RAND_MAX;

- (void)setUp
{
    [super setUp];
    
    expectedAccuracy = 1.0e-4;  // I expect this is sufficient for most practical uses of this class
    srand(42);  // Seed a consistent value to the random number generator.
    maxError = 0.0;
    //
    // Initialize test inputs.
    //
    scaleIn = shearRatiosIn = rotationAxisIn = translationIn = eyeIn = centerIn = upIn =invalidGLKVector3;
    frustrumBoundsIn = invalidGLKVector4;
    rotationAngleIn = fieldOfViewIn = aspectRatioIn = nearClipIn = farClipIn = 0.0;
    quaternionIn = invalidGLKQuaternion;
    //
    // Set up test transform.
    //
    transformIn = GLKMatrix4Identity;
    //
    // Initialize vectors to receive results.
    //
    scaleOut = shearRatiosOut = rotationOut = translationOut = invalidGLKVector3;
    quaternionOut = invalidGLKQuaternion;
    perspectiveOut = invalidGLKVector4;
    transformOut = GLKMatrix4Identity;
    //
    // Reset bypass flags
    //
    bypassScaleCheck = false;
    bypassTranslationCheck = false;
    bypassShearCheck = false;
    //
    // Add a space to the log between tests
    //
    NSLog(@" - ");
}

- (void)tearDown
{
    NSLog(@"Maximum element error magnitude was %f.", maxError);
    
    // Tear-down code here (none needed).
    
    [super tearDown];
}

#pragma mark - Core test support methods

- (bool)coreDecompositionTest
{
    bool result = false;
    [self composeInputTranform];
    float errorSize;
    result = [TUTransformUtilities decomposeTransform:transformIn intoScale:&scaleOut  shearRatios:&shearRatiosOut rotation:&rotationOut quaternion:&quaternionOut translation:&translationOut perspective:&perspectiveOut];
    if (result) {
        GLKVector3 delta;
        if (!bypassTranslationCheck) {
            if (!GLKVector3AllEqualToVector3(translationIn, invalidGLKVector3)) {
                delta = GLKVector3Subtract(translationIn, translationOut);
            }
            else {
                delta = translationOut;
            }
            errorSize = fabsf(GLKVector3Length(delta));
            if (errorSize > expectedAccuracy) {
                result = false;
                STFail(@"Difference in translations of %f exceeds expectations.", errorSize);
            }
        }
        if (!bypassScaleCheck) {
            if (!GLKVector3AllEqualToVector3(scaleIn, invalidGLKVector3)) {
                delta = GLKVector3Subtract(scaleIn, scaleOut);
            }
            else {
                delta = GLKVector3Subtract(unitGLKVector3,scaleOut);
            }
            errorSize = fabsf(GLKVector3Length(delta));
            if (errorSize > expectedAccuracy) {
                result = false;
                STFail(@"Difference in scale of %f exceeds expectations.", errorSize);
            }
        }
        if (!bypassShearCheck) {
            if (!GLKVector3AllEqualToVector3(shearRatiosIn, invalidGLKVector3)) {
                delta = GLKVector3Subtract(shearRatiosIn, shearRatiosOut);
            }
            else {
                delta = GLKVector3Subtract(zeroGLKVector3,shearRatiosOut);
            }
            errorSize = fabsf(GLKVector3Length(delta));
            if (errorSize > expectedAccuracy) {
                result = false;
                STFail(@"Difference in shear ratios of %f exceeds expectations.", errorSize);
            }
        }

        [self composeOutputMatrix];
        result = [self compareOutputTransformWithInput] && result;
        if (!result) {
            [self logTestConditions];
        }
    }
    else {
        STFail(@"Input transform seen as singular; false returned from decomposeTransform.");
    }
    return result;
}

- (void)composeInputTranform
{
    transformIn = GLKMatrix4Identity;
    if (!GLKVector3AllEqualToVector3(translationIn, invalidGLKVector3)) {
        transformIn = GLKMatrix4TranslateWithVector3(transformIn, translationIn);
    }
    if (!GLKVector3AllEqualToVector3(rotationAxisIn, invalidGLKVector3)) {
        transformIn = GLKMatrix4RotateWithVector3(transformIn, rotationAngleIn, rotationAxisIn);
        //
        // Quaternion constructor requires normalized axis.
        //
        quaternionIn = GLKQuaternionMakeWithAngleAndVector3Axis(rotationAngleIn, GLKVector3Normalize(rotationAxisIn));
    }
    if (!GLKVector3AllEqualToVector3(shearRatiosIn, invalidGLKVector3)) {
        GLKMatrix4 shearMatrix = GLKMatrix4Identity;
        shearMatrix.m21 = shearRatiosIn.z;  //YZ shear
        shearMatrix.m20 = shearRatiosIn.y;  //XZ shear
        shearMatrix.m10 = shearRatiosIn.x;  //XY shear
        transformIn = GLKMatrix4Multiply(transformIn, shearMatrix);
    }
    if (!GLKVector3AllEqualToVector3(scaleIn, invalidGLKVector3)) {
        transformIn = GLKMatrix4ScaleWithVector3(transformIn, scaleIn);
    }
//    NSLog(@"Input trans -proj: %@", NSStringFromGLKMatrix4(transformIn));
    if ([self composeInputProjection]) {
        transformIn = GLKMatrix4Multiply(projectionIn, transformIn);
    }
}

-(bool)composeInputProjection
{
    bool result = true;
    if (fieldOfViewIn != 0.0) {
        projectionIn = GLKMatrix4MakePerspective(fieldOfViewIn, aspectRatioIn, nearClipIn, farClipIn);
    }
    else if (!GLKVector4AllEqualToVector4(frustrumBoundsIn, invalidGLKVector4)) {
        projectionIn = GLKMatrix4MakeFrustum(frustrumBoundsIn.x, frustrumBoundsIn.y, frustrumBoundsIn.z, frustrumBoundsIn.w, nearClipIn, farClipIn);
    }
    else {
        projectionIn = GLKMatrix4Identity;
        result = false;
    }
    return result;
}

- (void)composeOutputMatrix
{
    transformOut = GLKMatrix4Identity;
    transformOut = GLKMatrix4TranslateWithVector3(transformOut, translationOut);
    GLKMatrix4 rotationMatrix = GLKMatrix4MakeWithQuaternion(quaternionOut);
    transformOut = GLKMatrix4Multiply(transformOut, rotationMatrix);
    //
    // The following code can be used to verify the extraction of Euler angles to represent rotation.  However,
    // the Euler angle approach is considerably less numerically accurate than using quaternions, and quaternion
    // or axis angle representations are likely to be more useful in most cases, so this code is disabled but
    // left for reference.
    //
    /*
    GLKMatrix4 rotationMatrix2 = GLKMatrix4Identity;
    rotationMatrix2 = GLKMatrix4RotateZ(rotationMatrix2, rotationOut.z);
    rotationMatrix2 = GLKMatrix4RotateY(rotationMatrix2, rotationOut.y);
    rotationMatrix2 = GLKMatrix4RotateX(rotationMatrix2, rotationOut.x);
    GLKMatrix4 differences = GLKMatrix4Subtract(rotationMatrix2, rotationMatrix);
    for (int i=1; i<16; i++) {
        STAssertEqualsWithAccuracy(differences.m[i], 0.0f, expectedAccuracy, @"Difference in transform element %i exceeds expectations", i);
        if (fabsf(differences.m[i]) > 1.0e-3) {
            NSLog(@"Difference exceeded");
        }
    }
    //
    // When using Euler rotations, this is one canonical order of composition.
    transformOut = GLKMatrix4RotateZ(transformOut, rotationOut.z);
    transformOut = GLKMatrix4RotateY(transformOut, rotationOut.y);
    transformOut = GLKMatrix4RotateX(transformOut, rotationOut.x);
    */
    GLKMatrix4 shearMatrix = GLKMatrix4Identity;
    shearMatrix.m21 = shearRatiosOut.z;  //YZ shear
    shearMatrix.m20 = shearRatiosOut.y;  //XZ shear
    shearMatrix.m10 = shearRatiosOut.x;  //XY shear
    transformOut = GLKMatrix4Multiply(transformOut, shearMatrix);
    transformOut = GLKMatrix4ScaleWithVector3(transformOut, scaleOut);
    // Compose the projection matrix representing the perspective output.
//    NSLog(@"perspective = %f, %f, %f, %f", perspectiveOut.x, perspectiveOut.y, perspectiveOut.z,perspectiveOut.w);
    projectionOut = GLKMatrix4Identity;
    projectionOut.m03 = perspectiveOut.x;
    projectionOut.m13 = perspectiveOut.y;
    // Adjust the z scaling to match output of GLKMatrix4MakeFrustrum and GLKMatrix4MakePerspective
    if (perspectiveOut.z != 0.0) {
        projectionOut.m22 = -1.0/perspectiveOut.z;
        projectionOut.m23 = -1.0;
        projectionOut.m32 = -1.0 * perspectiveOut.w / perspectiveOut.z;
        projectionOut.m33 = 0.0;
    }
    transformOut = GLKMatrix4Multiply(projectionOut, transformOut);
    
}

- (bool) compareOutputTransformWithInput
{
    bool result = true;
    float elementError = 0.0;
    GLKMatrix4 differences = GLKMatrix4Subtract(transformIn, transformOut);
    for (int i=1; i<16; i++) {
        STAssertEqualsWithAccuracy(differences.m[i], 0.0f, expectedAccuracy, @"Difference in transform element %i exceeds expectations", i);
        elementError = fabsf(differences.m[i]);
        if (elementError > maxError) {
            maxError = elementError;
        }
        if ( elementError > expectedAccuracy) {
            result = false;
        }
    }
    return result;
}

- (void) logTestConditions {
    NSString *logString = @"\n\nInputs:";
    logString = [logString stringByAppendingFormat:@"\nPerspective inputs - FOV: %f, apsect ratio: %f, near clip: %f, far clip: %f.", fieldOfViewIn, aspectRatioIn, nearClipIn, farClipIn];
    logString = [logString stringByAppendingFormat:@"\nQuaternion for input rotation = %@", NSStringFromGLKQuaternion(quaternionIn)];
    logString = [logString stringByAppendingFormat:@"\nscale = %f, %f, %f", scaleIn.x, scaleIn.y, scaleIn.z];
    logString = [logString stringByAppendingFormat:@"\ntranslation = %f, %f, %f", translationIn.x, translationIn.y, translationIn.z];
    logString = [logString stringByAppendingFormat:@"\n\nOutputs:"];
    logString = [logString stringByAppendingFormat:@"\nperspective = %f, %f, %f, %f", perspectiveOut.x, perspectiveOut.y, perspectiveOut.z,perspectiveOut.w];
    logString = [logString stringByAppendingFormat:@"\nshearRatios = %f, %f, %f", shearRatiosOut.x, shearRatiosOut.y, shearRatiosOut.z];
    logString = [logString stringByAppendingFormat:@"\nscale = %f, %f, %f", scaleOut.x, scaleOut.y, scaleOut.z];
    logString = [logString stringByAppendingFormat:@"\ntranslation = %f, %f, %f", translationOut.x, translationOut.y, translationOut.z];
    logString = [logString stringByAppendingFormat:@"\nquaternion = %@", NSStringFromGLKQuaternion(quaternionOut)];
    logString = [logString stringByAppendingFormat:@"\n\nMatrices:"];
    logString = [logString stringByAppendingFormat:@"\nInput projection:  %@", NSStringFromGLKMatrix4(projectionIn)];
    logString = [logString stringByAppendingFormat:@"\nInput transform:   %@", NSStringFromGLKMatrix4(transformIn)];
    logString = [logString stringByAppendingFormat:@"\nOutput projection: %@", NSStringFromGLKMatrix4(projectionOut)];
    logString = [logString stringByAppendingFormat:@"\nOutput transform:  %@", NSStringFromGLKMatrix4(transformOut)];
    NSLog(@"%@\n\n", logString);
}

#pragma mark - Level 1 tests
//
// Tests of resolution of independent transform types

- (void)testTranslationDecomposition
{
    for (int i=-100; i<=100; i+=5) {
        translationIn.x = (float)i;
        for (int j=-100; j<=100; j+=5) {
            translationIn.y = (float)j;
            for (int k=-100; k<=100; k+=5) {
                translationIn.z = (float)k;
                bool result = [self coreDecompositionTest];
                if (!result) return;
//                STAssertTrue(result, @"false returned from coreDecompositionTest, indicating decomposeTransform saw input transform as singular while testing translation %@", NSStringFromGLKVector3(translationIn));
            }
        }
    }
}

// Scale must be non-zero, and negative scales may be confounded with rotations
- (void)testScaleDecomposition
{
    for (int i=1; i<=100; i+=5) {
        scaleIn.x = (float)i/10.0;
        for (int j=1; j<=100; j+=5) {
            scaleIn.y = (float)j/10.0;
            for (int k=1; k<=100; k+=5) {
                scaleIn.z = (float)k/10.0;
                bool result = [self coreDecompositionTest];
                if (!result) return;
//                STAssertTrue(result, @"false returned from coreDecompositionTest, indicating decomposeTransform saw input transform as singular while testing scale %@", NSStringFromGLKVector3(scaleIn));
            }
        }
    }
}

- (void)testShearDecomposition
{
    for (int i=-100; i<=100; i+=5) {
        shearRatiosIn.x = (float)i/10.0;
        for (int j=-100; j<=100; j+=5) {
            shearRatiosIn.y = (float)j/10.0;
            for (int k=-100; k<=100; k+=5) {
                shearRatiosIn.z = (float)k/10.0;
                bool result = [self coreDecompositionTest];
                if (!result) return;
//                STAssertTrue(result, @"false returned from coreDecompositionTest, indicating decomposeTransform saw input transform as singular while testing hear ratios %@", NSStringFromGLKVector3(shearRatiosIn));
            }
        }
    }
}

- (void)testRotationDecomposition
{
    expectedAccuracy = 1.0e-5;
    for (rotationAngleIn = -M_PI; rotationAngleIn <= M_PI; rotationAngleIn += M_PI/360.0) { // Half degree increments
        for (int i=-20; i<=20; i+=5) {
            rotationAxisIn.x = (float)i;
            for (int j=-20; j<=20; j+=5) {
                rotationAxisIn.y = (float)j;
                for (int k=-20; k<=20; k+=5) {
                    rotationAxisIn.z = (float)k;
                    if ((i!=0)||(j!=0)||(j!=0)) {
                        bool result = [self coreDecompositionTest];
                        if (!result) return;
//                        STAssertTrue(result, @"false returned from coreDecompositionTest, indicating decomposeTransform saw input transform as singular while testing rotation around axis %@ with angle %f", NSStringFromGLKVector3(rotationAxisIn), rotationAngleIn);
                    }
                }
            }
        }
    }
}

- (void)testProjectionDecomposition
{
    translationIn = (GLKVector3){0.0, 0.0, 100.0};  // A non-zero translation is required for perspective.
    bypassScaleCheck = true;    // Projection and scale changes cannot be disambiguated.
    for (fieldOfViewIn = M_PI/18.0; fieldOfViewIn <= (17.0/18.0) * M_PI; fieldOfViewIn += M_PI/18.0) {
        for (aspectRatioIn = 0.25; aspectRatioIn <= 4.0; aspectRatioIn *= 2.0) {
            for (nearClipIn = 5; nearClipIn <= 50; nearClipIn += 5) {
                for (farClipIn = 100; farClipIn <= 100000; farClipIn *= 10) {
                    bool result = [self coreDecompositionTest];
                    if (!result) return;
//                    STAssertTrue(result, @"false returned from coreDecompositionTest, indicating decomposeTransform saw input transform as singular while testing projection with field of view: %f, aspect ratio: %f, near clip: %f, and far clip %f.",fieldOfViewIn, aspectRatioIn, nearClipIn, farClipIn );
                }
            }
        }
    }
    bypassScaleCheck = false;
}

#pragma mark - Level 2 tests
//
// Tests of resolution of pairs of transform types

- (void)testTranslationAndScaleDecomposition
{
    for (int i=-100; i<=100; i+=20) {
        translationIn.x = (float)i;
        for (int j=-100; j<=100; j+=20) {
            translationIn.y = (float)j;
            for (int k=-100; k<=100; k+=20) {
                translationIn.z = (float)k;
                for (int l=1; l<=100; l+=20) {
                    scaleIn.x = (float)l/10.0;
                    for (int m=1; m<=100; m+=20) {
                        scaleIn.y = (float)m/10.0;
                        for (int n=1; n<=100; n+=20) {
                            scaleIn.z = (float)n/10.0;
                            if (rand() < selectionThreshold) {
                                bool result = [self coreDecompositionTest];
                                if (!result) return;
//                                STAssertTrue(result, @"false returned from coreDecompositionTest, indicating decomposeTransform saw input transform as singular while testing translation %@ with scale %@", NSStringFromGLKVector3(translationIn), NSStringFromGLKVector3(scaleIn));
                            }
                        }
                    }
                }
            }
        }
    }
}

- (void)testTranslationAndShearDecomposition
{
    expectedAccuracy = 1.0e-3;
    for (int i=-100; i<=100; i+=20) {
        translationIn.x = (float)i;
        for (int j=-100; j<=100; j+=20) {
            translationIn.y = (float)j;
            for (int k=-100; k<=100; k+=20) {
                translationIn.z = (float)k;
                for (int l=1; l<=100; l+=20) {
                    shearRatiosIn.x = (float)l/10.0;
                    for (int m=1; m<=100; m+=20) {
                        shearRatiosIn.y = (float)m/10.0;
                        for (int n=1; n<=100; n+=20) {
                            shearRatiosIn.z = (float)n/10.0;
                            if (rand() < selectionThreshold) {
                                bool result = [self coreDecompositionTest];
                                if (!result) return;
//                                STAssertTrue(result, @"false returned from coreDecompositionTest, indicating decomposeTransform saw input transform as singular while testing translation %@ with shear %@", NSStringFromGLKVector3(translationIn), NSStringFromGLKVector3(shearRatiosIn));
                            }
                        }
                    }
                }
            }
        }
    }
}

- (void)testTranslationAndRotationDecomposition
{
    expectedAccuracy = 5.0e-4;
    for (int i=-100; i<=100; i+=25) {
        translationIn.x = (float)i;
        for (int j=-100; j<=100; j+=25) {
            translationIn.y = (float)j;
            for (int k=-100; k<=100; k+=25) {
                translationIn.z = (float)k;
                for (rotationAngleIn = -M_PI; rotationAngleIn <= M_PI; rotationAngleIn += M_PI/18.0) { // Ten degree increments
                    for (int l=-20; l<=20; l+=5) {
                        rotationAxisIn.x = (float)l;
                        for (int m=-20; m<=20; m+=5) {
                            rotationAxisIn.y = (float)m;
                            for (int n=-20; n<=20; n+=5) {
                                rotationAxisIn.z = (float)n;
                                if (rand() < selectionThreshold) {
                                    if ((l!=0)||(m!=0)||(n!=0)) {
                                        bool result = [self coreDecompositionTest];
                                        if (!result) return;
//                                        STAssertTrue(result, @"false returned from coreDecompositionTest, indicating decomposeTransform saw input transform as singular while testing translation %@ with rotation around axis %@ with angle %f", NSStringFromGLKVector3(translationIn), NSStringFromGLKVector3(rotationAxisIn), rotationAngleIn);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

- (void)testTranslationAndProjectionDecomposition
{
    bypassScaleCheck = true;    // Projection and scale changes cannot be disambiguated?
    bypassTranslationCheck = true;
    bypassShearCheck = true;
    expectedAccuracy = 5.0e-4;
    for (int i=-100; i<=100; i+=25) {
        translationIn.x = (float)i;
        for (int j=-100; j<=100; j+=25) {
            translationIn.y = (float)j;
            for (int k=100; k<=200; k+=25) { // Keep y translation positive for comatibility with projection
                translationIn.z = (float)k;
                for (fieldOfViewIn = M_PI/18.0; fieldOfViewIn <= (17.0/18.0) * M_PI; fieldOfViewIn += M_PI/18.0) {
                    for (aspectRatioIn = 0.25; aspectRatioIn <= 4.0; aspectRatioIn *= 2.0) {
                        for (nearClipIn = 5; nearClipIn <= 50; nearClipIn += 5) {
                            for (farClipIn = 100; farClipIn <= 100000; farClipIn *= 10) {
                                if (rand() < selectionThreshold) {
                                    bool result = [self coreDecompositionTest];
                                    if (!result) return;
//                                    STAssertTrue(result, @"false returned from coreDecompositionTest, indicating decomposeTransform saw input transform as singular while testing translation %@ with projection with field of view: %f, aspect ratio: %f, near clip: %f, and far clip %f.",NSStringFromGLKVector3(translationIn), fieldOfViewIn, aspectRatioIn, nearClipIn, farClipIn );
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    bypassScaleCheck = false;
    bypassTranslationCheck = false;
    bypassShearCheck = false;
}

// Scale must be non-zero, and negative scales may be confounded with rotations
- (void)testScaleAndShearDecomposition
{
    for (int i=1; i<=100; i+=25) {
        scaleIn.x = (float)i/10.0;
        for (int j=1; j<=100; j+=25) {
            scaleIn.y = (float)j/10.0;
            for (int k=1; k<=100; k+=25) {
                scaleIn.z = (float)k/10.0;
                for (int l=-100; l<=100; l+=25) {
                    shearRatiosIn.x = (float)l/10.0;
                    for (int m=-100; m<=100; m+=25) {
                        shearRatiosIn.y = (float)m/10.0;
                        for (int n=-100; n<=100; n+=25) {
                            shearRatiosIn.z = (float)n/10.0;
                            if (rand() < selectionThreshold) {
                                bool result = [self coreDecompositionTest];
                                if (!result) return;
//                                STAssertTrue(result, @"false returned from coreDecompositionTest, indicating decomposeTransform saw input transform as singular while testing scale %@ and shear ratios %@", NSStringFromGLKVector3(scaleIn), NSStringFromGLKVector3(shearRatiosIn));
                            }
                        }
                    }
                }
            }
        }
    }
}

- (void)testScaleAndRotationDecomposition
{
    expectedAccuracy = 5.0e-4;
    for (int i=1; i<=100; i+=25) {
        scaleIn.x = (float)i/10.0;
        for (int j=1; j<=100; j+=25) {
            scaleIn.y = (float)j/10.0;
            for (int k=1; k<=100; k+=25) {
                scaleIn.z = (float)k/10.0;
                for (rotationAngleIn = -M_PI; rotationAngleIn <= M_PI; rotationAngleIn += M_PI/18.0) { // Ten degree increments
                    for (int l=-20; l<=20; l+=5) {
                        rotationAxisIn.x = (float)l;
                        for (int m=-20; m<=20; m+=5) {
                            rotationAxisIn.y = (float)m;
                            for (int n=-20; n<=20; n+=5) {
                                rotationAxisIn.z = (float)n;
                                if (rand() < selectionThreshold) {
                                    if ((l!=0)||(m!=0)||(n!=0)) {
                                        bool result = [self coreDecompositionTest];
                                        if (!result) return;
//                                        STAssertTrue(result, @"false returned from coreDecompositionTest, indicating decomposeTransform saw input transform as singular while testing scale %@ with rotation around axis %@ with angle %f", NSStringFromGLKVector3(scaleIn), NSStringFromGLKVector3(rotationAxisIn), rotationAngleIn);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

- (void)testScaleAndProjectionDecomposition
{
    translationIn = (GLKVector3){0.0, 0.0, 100.0};  // A non-zero translation is required for perspective.
    bypassScaleCheck = true;    // Projection and scale changes cannot be disambiguated.
    for (int i=1; i<=100; i+=25) {
        scaleIn.x = (float)i/10.0;
        for (int j=1; j<=100; j+=25) {
            scaleIn.y = (float)j/10.0;
            for (int k=1; k<=100; k+=25) {
                scaleIn.z = (float)k/10.0;
                for (fieldOfViewIn = M_PI/18.0; fieldOfViewIn <= (17.0/18.0) * M_PI; fieldOfViewIn += M_PI/18.0) {
                    for (aspectRatioIn = 0.25; aspectRatioIn <= 4.0; aspectRatioIn *= 2.0) {
                        for (nearClipIn = 5; nearClipIn <= 50; nearClipIn += 5) {
                            for (farClipIn = 100; farClipIn <= 100000; farClipIn *= 10) {
                                if (rand() < selectionThreshold) {
                                    bool result = [self coreDecompositionTest];
                                    if (!result) return;
//                                    STAssertTrue(result, @"false returned from coreDecompositionTest, indicating decomposeTransform saw input transform as singular while testing scale %@ with projection with field of view: %f, aspect ratio: %f, near clip: %f, and far clip %f.",NSStringFromGLKVector3(scaleIn), fieldOfViewIn, aspectRatioIn, nearClipIn, farClipIn );
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    bypassScaleCheck = false;
}

- (void)testShearAndRotationDecomposition
{
    expectedAccuracy = 5.0e-4;
    for (int i=-100; i<=100; i+=25) {
        shearRatiosIn.x = (float)i/10.0;
        for (int j=-100; j<=100; j+=25) {
            shearRatiosIn.y = (float)j/10.0;
            for (int k=-100; k<=100; k+=25) {
                shearRatiosIn.z = (float)k/10.0;
                for (rotationAngleIn = -M_PI; rotationAngleIn <= M_PI; rotationAngleIn += M_PI/18.0) { // Ten degree increments
                    for (int l=-20; l<=20; l+=5) {
                        rotationAxisIn.x = (float)l;
                        for (int m=-20; m<=20; m+=5) {
                            rotationAxisIn.y = (float)m;
                            for (int n=-20; n<=20; n+=5) {
                                rotationAxisIn.z = (float)n;
                                if (rand() < selectionThreshold) {
                                    if ((l!=0)||(m!=0)||(n!=0)) {
                                        bool result = [self coreDecompositionTest];
                                        if (!result) return;
//                                        STAssertTrue(result, @"false returned from coreDecompositionTest, indicating decomposeTransform saw input transform as singular while testing shear ratios %@ with rotation around axis %@ with angle %f", NSStringFromGLKVector3(shearRatiosIn), NSStringFromGLKVector3(rotationAxisIn), rotationAngleIn);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

- (void)testShearAndProjectionDecomposition
{
    expectedAccuracy = 0.1;  //   Much lower numeric stability on device!
    translationIn = (GLKVector3){0.0, 0.0, 100.0};  // A non-zero translation is required for perspective.
    bypassScaleCheck = true;    // Projection and scale changes cannot be disambiguated.
    bypassShearCheck = true;
    for (int i=-100; i<=100; i+=10) {
        shearRatiosIn.x = (float)i/10.0;
        for (int j=-100; j<=100; j+=10) {
            shearRatiosIn.y = (float)j/10.0;
            for (int k=-100; k<=100; k+=10) {
                shearRatiosIn.z = (float)k/10.0;
                for (fieldOfViewIn = M_PI/18.0; fieldOfViewIn <= (17.0/18.0) * M_PI; fieldOfViewIn += M_PI/18.0) {
                    for (aspectRatioIn = 0.25; aspectRatioIn <= 4.0; aspectRatioIn *= 2.0) {
                        for (nearClipIn = 5; nearClipIn <= 50; nearClipIn += 5) {
                            for (farClipIn = 100; farClipIn <= 100000; farClipIn *= 10) {
                                if (rand() < selectionThreshold) {
                                    bool result = [self coreDecompositionTest];
                                    if (!result) return;
//                                    STAssertTrue(result, @"false returned from coreDecompositionTest, indicating decomposeTransform saw input transform as singular while testing shear ratios %@ with projection with field of view: %f, aspect ratio: %f, near clip: %f, and far clip %f.",NSStringFromGLKVector3(shearRatiosIn), fieldOfViewIn, aspectRatioIn, nearClipIn, farClipIn );
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    bypassScaleCheck = false;
}

- (void)testRotationAndProjectionDecomposition
{
    expectedAccuracy = 5.0e-4;
    bypassShearCheck = true;
    translationIn = (GLKVector3){0.0, 0.0, 100.0};  // A non-zero translation is required for perspective.
    bypassScaleCheck = true;    // Projection and scale changes cannot be disambiguated.
    for (rotationAngleIn = -M_PI; rotationAngleIn <= M_PI; rotationAngleIn += M_PI/18.0) { // Ten degree increments
        for (int i=-20; i<=20; i+=5) {
            rotationAxisIn.x = (float)i;
            for (int j=-20; j<=20; j+=5) {
                rotationAxisIn.y = (float)j;
                for (int k=-20; k<=20; k+=5) {
                    rotationAxisIn.z = (float)k;
                    if ((i!=0)||(j!=0)||(j!=0)) {
                        for (fieldOfViewIn = M_PI/18.0; fieldOfViewIn <= (17.0/18.0) * M_PI; fieldOfViewIn += M_PI/18.0) {
                            for (aspectRatioIn = 0.25; aspectRatioIn <= 4.0; aspectRatioIn *= 2.0) {
                                for (nearClipIn = 5; nearClipIn <= 50; nearClipIn +=10) {
                                    for (farClipIn = 150; farClipIn <= 100000; farClipIn *= 20) {
                                        if (rand() < selectionThreshold) {
                                            bool result = [self coreDecompositionTest];
                                            if (!result) return;
//                                            STAssertTrue(result, @"false returned from coreDecompositionTest, indicating decomposeTransform saw input transform as singular while testing rotation around axis %@ with angle %f with projection with field of view: %f, aspect ratio: %f, near clip: %f, and far clip %f.",NSStringFromGLKVector3(rotationAxisIn), rotationAngleIn, fieldOfViewIn, aspectRatioIn, nearClipIn, farClipIn );
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    bypassScaleCheck = false;
}

#pragma mark - Debug tests

- (void)testMatrixDecomposition
{
    translationIn = (GLKVector3){0.0, 0.0, 100.0};
//    shearRatiosIn = (GLKVector3){0.5, 0.0, 0.0};
//    rotationAngleIn = M_PI_2/2.0;
//    rotationAxisIn = (GLKVector3){1.0, 1.0, 0.0};
//    scaleIn = (GLKVector3){1.0, 2.0, 1.0};
//    fieldOfVewIn = GLKMathDegreesToRadians(35.0f);
//    aspectRatioIn = 1.0;
//    nearClipIn = 50.0f;
//    farClipIn = 1000.0f;   //(M_PI_2, 1.0, 10.0, 10000.0);
//
    frustrumBoundsIn = (GLKVector4){-4.0, 16.0, -5.0, 25.0};
//    fieldOfViewIn = GLKMathDegreesToRadians(90.0f);
//    aspectRatioIn = 1.0;
    nearClipIn = 10.0f;
    farClipIn = 100.0f;
//    eyeIn = (GLKVector3){0.0, 0.0, 1.0};
//    centerIn = (GLKVector3){0.0, 0.0, 0.0};
//    upIn = (GLKVector3){0.0, 1.0, 0.0};
//    bool result =
    [self coreDecompositionTest];
//    STAssertTrue(result, @"false returned from coreDecompositionTest.");
}

@end
