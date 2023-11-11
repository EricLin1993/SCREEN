function [ KC,yc  ] = IPSVDCompress( K,y, n )
% IPSVDCOMPRESS employs SVD (Singular Value Decomposition) to compress
% the Inverse Problem  K*x = y as Kc*x = yc with n principal components for the purpose of saving
% time and memory.
% Author: Enping Lin, (646332908@qq.com)  
    [u,sv,v]=svd(K,'econ');
%   figure,plot(diag(sv),'ro');
    u = u(:,1:n);
    v = v(:,1:n);
    yc = u'*y;
    KC = u'*K; 

end

