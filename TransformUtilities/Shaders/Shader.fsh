//
//  Shader.fsh
//  TransformUtilities
//
//  Created by Steven McGrath on 5/8/13.
//  Copyright (c) 2013 Steven McGrath. All rights reserved.
//

varying lowp vec4 colorVarying;

void main()
{
    gl_FragColor = colorVarying;
}
