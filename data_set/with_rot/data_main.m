
% create init dataset
clear all
close all
clc
load data_set\without_rot\pickup.mat
%%
ang = 60;
dataname = 'pickup';
S_yoga_60_2 = addRotation(ang,S,dataname,list,1); % para plotear, añadir el flag de 1
save('./data_set/with_rot/S_yoga_60_2',"S_yoga_60_2")
%%

ang = 90;
dataname = 'pickup';
S_pickup_90_2 = addRotation(ang,S,dataname,list);
save('./data_set/with_rot/S_pickup_90_2',"S_pickup_90_2")
%%

ang = 120;
dataname = 'pickup';
S_yoga_120_2 = addRotation(ang,S,dataname,list);
save('./data_set/with_rot/S_yoga_120_2',"S_yoga_120_2")


%%

% adapt data form
clear all
close all
clc
load data_set\with_rot\S_yoga_60_2.mat
X_yoga_60_2 = adapt_rearrange(S_yoga_60_2);
save('./data_set/with_rot/adapt_datas/X_yoga_60_2',"X_yoga_60_2")
%%
load data_set\with_rot\S_pickup_90_2.mat
X_pickup_90_2 = adapt_rearrange(S_pickup_90_2);
save('./data_set/with_rot/adapt_datas/X_pickup_90_2',"X_pickup_90_2")
%%
load data_set\with_rot\S_yoga_120_2.mat
X_yoga_120_2 = adapt_rearrange(S_yoga_120_2);
save('./data_set/with_rot/adapt_datas/X_yoga_120_2',"X_yoga_120_2")


