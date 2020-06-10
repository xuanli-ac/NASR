function [upbd_all,upbd_min] = get_init_lambda ( X )

%This function calculates the upper bound of lambda.
%The value of lambda used in NASR should be smaller than the upperbound.

% Input: X is a t*d time series, t time points and d channels/ROIs.
% Output: lambda_upperbound, upper bound for each sample in X

[dim,num] = size(X) ;
upbd_all=zeros(num,1);

for i = 1 : num
    y = X(:,i) ;
    if i == 1
        ind = 2:num ;
    elseif i == num
        ind = 1:num-1 ;
    else
        ind = [1:i-1 i+1:num] ;
    end
    riX = X(:,ind) ; 

    
% calculate the upper bound for each sample in X  
    [U1,sigma1,V1] = svd(riX);
    op_n = max(diag(sigma1));
    inf_n=norm(riX'*y,Inf);
    upbd_all(i)=op_n*inf_n;
end
upbd_min=min(upbd_all);   
