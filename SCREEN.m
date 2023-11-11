function [ ParaOut ] = SCREEN( ParaIn )
% SCREEN is a general sparse sampling NMR spectroscopy reconstruction
% method based on the Compressed Sensing theory
% Author: Enping Lin   (646332908@qq.com)
%-------------------------------------------------------------------------
K = ParaIn.K;
KMatrixYes = logical(ParaIn.KMatrixYes);% 1: matrix ; 0£ºOperator
KSize = ParaIn.KSize;
if logical(ParaIn.KComplexYes)    
    if ~isfield(ParaIn,'beta')
      beta = 0.5;
    else
      beta = ParaIn.beta ;
    end  
    if KMatrixYes
       K = [real(K),imag(K);-imag(K),real(K)];    
    end
    sn = [real(ParaIn.s);imag(ParaIn.s)];
else
   sn = ParaIn.s; 
end 

if ~isfield(ParaIn,'lambda')
    lambda = 0.1;
else
    lambda = ParaIn.lambda ;
end   

if ~isfield(ParaIn,'rel_tol')
    rel_tol = 1e-5;
else
    rel_tol = ParaIn.rel_tol;
end 
if ~isfield(ParaIn,'nonneg')
    nonneg = 0;
else
    nonneg = ParaIn.nonneg;
end 

if logical(nonneg)  % Check whether non-negative constraint on x is required
  if KMatrixYes  
    [x,xdiff]=l1_ls_nonneg(K,sn,lambda,rel_tol);
  else
    [x,xdiff]=l1_ls_nonneg(K,K',KSize(1),KSize(2),sn,lambda,rel_tol);  
  end
else
  if KMatrixYes  
    [x,xdiff]=l1_ls(K,sn,lambda,rel_tol);
  else
    [x,xdiff]=l1_ls(K,K',KSize(1),KSize(2),sn,lambda,rel_tol);  
  end
end

if logical(ParaIn.KComplexYes)  
   d = length(x)/2;
   ParaOut.x = beta*x(1:d)+ 1/(1-beta)*1j*x(d+1:end);
else
   ParaOut.x = x;  
end
ParaOut.xdiff = xdiff;
end

