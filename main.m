%main function to calculate the functional connectivity matrix  by using
%cal_assocmatrix

% Input:
% Dn: N*1 cell, data matrix of normalized fMRI series of N subjects
% lamlist: set the values of lambda to be tested

% Output:  
% M: Association matrix
% Para: obj&reconstruction error, for results control
clc;
clear;


lamlist=[0.1,0.05,0.01,0.12,0.25]; %regularizer parameter lambda



load('ROISignals/Dn.mat'); %Dn:Nsub*Ntime*Nnode
Nsub=length(Dn);
disp(['Total number of subjects is ',num2str(Nsub)]);

NewFolder='Results_M'; %Save the results
if ~exist(NewFolder,'dir')
mkdir(NewFolder);
end



for lam=lamlist
    disp(['Start lambda=',num2str(lam),' for ',num2str(Nsub),' subjects!'])
    [M,Para]=cal_assocmatrix(Dn,lam);
    save(['Results_M/M_l_',num2str(lam),'.mat'],'M');
    save(['Results_M/Para_l_',num2str(lam),'.mat'],'Para');
    
end