% Download and rearrange NRSfM data sets.
%
% Data sets:
%
% shark, face, walking sequences from http://www.cs.dartmouth.edu/~lorenzo/nrsfm.html
% Ref: L. Torresani, A. Hertzmann, and C. Bregler,
% "Nonrigid Structure-from-Motion: Estimating Shape and Motion with Hierarchical Priors,"
% IEEE Trans. Pattern Analysis and Machine Intelligence, vol. 30, no. 5, pp. 878–892, May 2008.
%
% drink, pickup, yoga, stretch, dance sequences from http://cvlab.lums.edu.pk/nrsfm/
% Ref: I. Akhter, Y. Seikh, S. Khan, and T. Kanade,
% "Trajectory Space: A Dual Representation for Nonrigid Structure from Motion,"
% IEEE Trans. Pattern Analysis and Machine Intelligence, vol. 33, no. 7, pp. 1442–1456, July 2011.
%
% FRGC manual annotation from http://vml.hanyang.ac.kr/research/procrustean/
% Ref: Minsik Lee, Jungchan Cho, Chong-Ho Choi, and Songhwai Oh,
% "Procrustean Normal Distribution for Non-Rigid Structure from Motion,"
% CVPR 2013, Portland, Oregon, June 23-28, 2013.
%
% Implemented by Minsik Lee (mlee.paper@gmail.com)
% Last update: 2016-09-07


clear; close all; clc;

disp('Downloading... (This may take some time. Please wait.)');

mkdir('Data');


% Data from Torresani et al.
seq = {'jaws', 'face', 'walking'};
url = 'http://www.cs.dartmouth.edu/~lorenzo/Data/';

for k=1:numel(seq)
    if strcmp(seq{k}, 'jaws')
        unzip([url 'shark.zip']);
    else
        unzip([url seq{k} '.zip']);
    end
    load(seq{k});
    delete([seq{k} '.mat']);
    [kn, p] = size(P3_gt);
    X = reshape(reshape(P3_gt, [], 3*p)', 3, p, []);
    save(['Data/' seq{k} '_rearranged.mat'], 'X');
    disp([seq{k} ' done.']);
end
delete('runme.m');


% Data from Akhter et al.
seq = {'drink', 'pickup', 'yoga', 'stretch', 'dance'};
url = 'https://cvlab.lums.edu.pk/nrsfm_/';

for k=1:numel(seq)
    while true
        try
            websave([seq{k} '.mat'], [url seq{k} '.mat'], weboptions('KeyName', 'Authorization', 'KeyValue', 'mykey', 'TimeOut', 30));
            load(seq{k});
            break;
        catch
            pause(10);
        end
    end
    delete([seq{k} '.mat']);
    [kn, p] = size(S);
    n = kn/3;
    
    if strcmp(seq{k}, 'dance')
        Rs = repmat([0 1 0; 0 0 1; 1 0 0], [1 1 n]);
    else
        temp = permute(reshape(Rs, 2, n, 3), [1 3 2]);
        Rs = zeros(3, 3, n);
        Rs(1:2, :, :) = temp;
        
        for i=1:n
            Rs(3, :, i) = cross(Rs(1, :, i), Rs(2, :, i));
        end
    end
    
    X = permute(reshape(S, 3, n, p), [1 3 2]);
    for i=1:n
        X(:, :, i) = Rs(:, :, i)*X(:, :, i);
        if mse(W(2*i-[1 0], :) - X(1:2, :, i)) ~= 0
            disp('err');
        end
    end
    save(['Data/' seq{k} '_rearranged.mat'], 'X');
    disp([seq{k} ' done.']);
end


% Manual annotation of FRGC 2.0 database
unzip('http://vml.hanyang.ac.kr/wp-content/uploads/2016/09/FRGC_annotation.zip');
movefile('FRGC_annotation.mat','Data/FRGC_rearranged.mat');
disp('FRGC done.');

clear;
disp('Completed.');

