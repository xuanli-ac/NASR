function [M,Para]=cal_assocmatrix(Dn,lam)
% calculate association matrices by NASR

%input: D_n normalized fMRI data, a N*1 cell containing the time series of N subjects.
% D_n{i} is t*d, where t is the number of time points and d it the number of
% ROIs/channels/nodes/voxels
% 
% lam:value of lambda (parameter)
% Output:M NASR association matrix

Nsub=length(Dn);
M=cell(1,Nsub);
Para=cell(1,Nsub);


for i=1:Nsub
    [Coeff,paravec] = solve_nntl( Dn{i},lam,0.03, 0);  % l=0.1,mu=2.5/90
    W=0.5*(Coeff+Coeff');
    M{i}=W;
    Para{i}=paravec;
    disp(['Subject', num2str(i),'with lambda=',num2str(lam),'done!']);
end


