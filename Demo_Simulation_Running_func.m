%% simulation: generate the two fourth-order tensors of space, 
%% topo frequency, time and subject
clear all
clc
close all
startup; % import EEGLab & tensor toolbox
load 'SyntheticData.mat'
%% Mode1-Topo: Common 2 + Indiv-SynHC 1 + Indiv-SynMDD 2
%Utopo = bsxfun(@rdivide,Utopo,sqrt(sum(Utopo.^2)));
U0{1}{1} = Utopo(:,[1 2 3]);
U0{2}{1} = Utopo(:,[1 2 4 5]);
%% Mode2-Spectrum: Common 2 + Indiv-SynHC 1 + Indiv-SynMDD 2
%Uspec = bsxfun(@rdivide,Uspec,sqrt(sum(Uspec.^2)));
U0{1}{2} = Uspec(:,[1 2 3]);
U0{2}{2} = Uspec(:,[1 2 4 5]);
%% Mode3-Time: Common 2 + Indiv-SynHC 1 + Indiv-SynMDD 2
%Utime = bsxfun(@rdivide,Utime,sqrt(sum(Utime.^2)));
U0{1}{3} = Utime(:,[1 2 3]);
U0{2}{3} = Utime(:,[1 2 4 5]);
%% Mode4-Subject: SynHC-19 SynMDD-20
U0{1}{4} = Usub(1:19,[1 2 3]);
U0{2}{4} = Usub(20:end,[1 2 3 4]);
%% generate the synthetic fourth-order tensors
noiseSynHC  = tensor(rand(64,130,500,19));
noiseSynMDD = tensor(rand(64,130,500,20));
signalSynHC  = full(ktensor(U0{1}));
signalSynMDD = full(ktensor(U0{2}));

Z.snr = 10;
Z.sig = 1;
Z.noi = Z.sig/(10^(Z.snr/10));
Z.object{1} = Z.sig*signalSynHC + Z.noi*noiseSynHC/norm(noiseSynHC)*norm(signalSynHC); % SynHC-group
Z.object{2} = Z.sig*signalSynMDD + Z.noi*noiseSynMDD/norm(noiseSynMDD)*norm(signalSynMDD); % SynMDD-group
Z.R = [3 4];
Z.C = [2 2 2 0];
Z.tol  = 1e-6;
Z.maxIter = 1000;
Z.U0 = U0;
Z.eps = 1e-16;
%% perform the coupled factorization
Mont = 1;
for mont = 1:Mont
    outt{mont} = f_NCTF_ADMM(Z);
end
filename = 'results.mat';
save(filename,'outt','U0')