function [Coeff,paravec] = solve_nntl( X , lambda ,mu0, display  )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This routine solves the following Non-negative trace lasso optimization problem
% for each data point (probe sample) y=x_i (i=1,...,n) 
% based on X=[x_1,...,x_{i-1},x_{i+1},...x_n] as the dictionary

% Formulation: min 0.5 * ||y-Xw||_2^2 + lambda * ||X*Diag(w)||_*

% The optimization problem can be solve by using Alternating Direction Method (ADM)
% We add a non-negative constraint on w, thus make a little adjustment on the optimization procedure
% 4 input augments

%%%%%%%%% Input: X: m timepoints by n nodes, each column is normalized 
%%%%%%%%%        lambda & mu0 : two parameters involved in the optimization
%%%%%%%%%        display: value 1 or 0, 1: show the results     

%%%%%%%%% Output: 
%%%%%%%%% Coeff: Coefficient matrix (n by n)
%%%%%%%%% paravec: save the values of reconstruction error, regularizer, and
%%%%%%%%% the objective function for optimization (.err, .norm, .obj)

%%%%%%%%% By Xuan Li 2015/12/09

if nargin<4
    display = false ;
end
if nargin < 2
    lambda = 1 ;
end

% if nargin < 3
%     mu = 1 ;
% end

[dim,num] = size(X) ;
Coeff = zeros( num , num ) ;
tol = 1e-8 ;
maxIter = 1000 ; % # of max iteration
tol2 = 1e-4 ;
rho0 = 1.1 ;
max_mu = 1e10 ;

for i = 1 : num
    y = X(:,i) ; % timeseries of the i th node as the probe sample
    if i == 1
        ind = 2:num ;
    elseif i == num
        ind = 1:num-1 ;
    else
        ind = [1:i-1 i+1:num] ;
    end
    riX = X(:,ind) ; % Dictionary
    XtX = riX'*riX ;
    diagXtX = diag(diag(XtX)) ;
    Xty = riX'*y ;
    
    
    %% Initializing optimization variables
    mu = mu0 ; %2.5/min(num,dim) ;
    w = zeros(num-1,1) ; % coefficient vector
%     w = inv(riX'*riX+lambda*eye(num-1,num-1)) * riX' * y ;
    Z = zeros(dim,num-1) ;
    Y = zeros(dim,num-1) ;
    iter = 0 ;   
    
    
   %% Updating variables  
    while iter<maxIter
        iter = iter + 1; 
        w_old = w ;
        Z_old = Z ;       
        
        %update Z
        temp = riX*diag(w) - Y/mu ;
        [U,sigma,V] = svd(temp,'econ');
        sigma = diag(sigma);
        svp = length( find( sigma>lambda/mu ) ) ;
        if svp>=1
            sigma = sigma(1:svp)-lambda/mu ;
        else
            svp = 1 ;
            sigma = 0 ;
        end
        Z = U(:,1:svp)*diag(sigma)*V(:,1:svp)' ;

        %udpate w        
        A = XtX + mu*diagXtX ;
        b = Xty + diag((Y+mu*Z)'*riX) ;
        w = A\b ;                                     % \ means inverse(A)*b
       
        %%%%% non-negative constraint: w=max(A\b,0), according to the KKT condition
        w(w<0)=0;                                
       
  
        ymxw = y - riX*w ;                            % y-Xw
        leq = Z - riX*diag(w) ;                       % J-XDiag(w)      
        stopC = max(max(abs(leq))) ;
        err = norm(ymxw) ;
        reg = nuclearnorm(riX*diag(w)) ;
        obj(iter) = 0.5*err^2 + lambda*reg + trace(Y'*leq)+mu*norm(leq,'fro')/2 ;   %%%%% calculate the value of objective function
        
        %% display
        if display && (iter==1 || mod(iter,1)==0 || stopC<tol)
            err = norm(ymxw) ;                         % reconstruction error
            reg = nuclearnorm(riX*diag(w)) ;
            obj(iter) = 0.5*err^2 + lambda*reg + trace(Y'*leq)+mu*norm(leq,'fro')/2 ;
%             obj(iter) = 0.5*err^2 + lambda*reg ;
            disp(['iter ' num2str(iter) ',mu=' num2str(mu,'%2.1e') ...
                ',rank=' num2str(rank(riX*diag(w))) ',stopALM=' num2str(stopC,'%2.3e') ...
                ',err=' num2str(err) ',norm=' num2str(reg) ',obj=' num2str(obj(iter)) ]);
        end
        
        if stopC<tol
            break;
        else
            Y = Y + mu*leq;
            if max( max(abs(w-w_old)) , max(max(abs(Z-Z_old))) ) > tol2
                rho = rho0 ;
            else
                rho = 1 ;
            end
            mu = min(max_mu,mu*rho);
%             mu = min(max_mu,mu*rho(iter));
        end 
         
    end   
    
    %% saving results 

    paravec.err=err;
    paravec.norm=reg;
    paravec.obj=obj(iter);

    Coeff(ind,i) = w ;
    
    % show progress
    disp([ 'solve' num2str(i) 'sample'] );
    
end


% calculate the nuclear norm, |X|_*
function nnorm = nuclearnorm( X )
s = svd(X) ;
nnorm = sum( s ) ;
