Copyright (C) 2016 Minsik Lee

    Division of Electrical Engineering
    Hanyang University, Korea.

This code is an implementation of the methods described in:

    Minsik Lee, Jungchan Cho, and Songhwai Oh,
    "Consensus of Non-Rigid Reconstructions,"
    CVPR 2016, Las Vegas, Nevada, June 26-July 1, 2016.

This software is distributed WITHOUT ANY WARRANTY. Use of this software is 
granted for research conducted at research institutions only. Commercial use
of this software is not allowed. Corporations interested in the use of this
software should contact the authors. If you use this code for a scientific
publication, please cite the above paper.

USAGE:

Please see the demo files "demo.m" for usage information. These scripts were
tested with MATLAB versions R2015b.

FEEDBACK:

Your feedback is greatly welcome. Please send bug reports, suggestions, and/or
new results to:

    mlee.paper@gmail.com

In the future, updated versions of this software will be available at:

    http://vml.hanyang.ac.kr/

CONTENTS:

    README.txt:                         This file.
    gpl.txt:                            License information.
    demo.m:                             Demo program.
    NRSFM_Consensus.m:                  Implementation of the method.
    select_idx.m:                       Sample trajectory groups.
    reconstruct.m:                      Weak reconstructor.
    part_reflection.m:                  Resolve reflection ambiguities between groups.
    combine.m:                          Obtain strong reconstruction.
    pout_trans.m:                       Eliminate translation component.
    mse.m:                              Calculate mean squared error.
