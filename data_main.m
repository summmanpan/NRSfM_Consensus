
seq = 'dinosaur_real'; % face_mocap; face_real; dinosaur_real; ogre_synthetic

% load(fprintf('./data_set/falta_rearenged/%s.mat',seq))
load(['./data_set/data_cv/falta_rearenged/' seq]) 

X1 = GT.shape3D;
X = adapt_rearrange(X1,seq);

%% LOAD DATA SET

clear all
close all
clc
load data_set\without_rot\yoga.mat
dataname = 'yoga';
%%
ang = 60;

S_60 = addRotation(ang,S,dataname,list,0); % para plotear, añadir el flag de 1
seq = ['rot_60/' dataname '_60'];
X = adapt_rearrange(S_60,seq,1);

%%

ang = 90;
S_90 = addRotation(ang,S,dataname,list,0);
seq = ['rot_90/' dataname '_90'];
X = adapt_rearrange(S_90,seq,1);

%%

ang = 120;
S_120 = addRotation(ang,S,dataname,list,0);
seq = ['rot_120/' dataname '_120'];
X = adapt_rearrange(S_120,seq,1);


