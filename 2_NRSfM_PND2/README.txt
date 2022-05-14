Copyright (C) 2016 Minsik Lee

    Division of Electrical Engineering
    Hanyang University, Korea.

This code is an implementation of the methods described in:

    Minsik Lee, Jungchan Cho, and Songhwai Oh,
    "Procrustean Normal Distribution for Non-Rigid Structure from Motion,"
    IEEE Trans. Pattern Analysis and Machine Intelligence, to appear.

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

    http://vml.hanyang.ac.kr/research/procrustean/

CONTENTS:

    README.txt:                         This file.
    gpl.txt:                            License information.
    download_and_rearrange_data.m       Download and rearrange test data sets.
    demo.m:                             Demo program.
    NRSFM_PND2.m:                       Implementation of the method.
    initialize.m:                       Calculate initial rotations.
    GPTA.m:                             Modified generalized Procrustes Analysis.
    EM_PND2.m:                          Fit PND to data using EM algorithm.
    update_parameters:                  Update parameters for M-step.
    plot_NRSfM.m:                       Plot reconstructed results.
    mse.m:                              Calculate mean squared error.
