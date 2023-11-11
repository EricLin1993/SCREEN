function [A_final,X,diffe]=fista(A, integration_corrected, lambda, matAl,iterk)
t=1/max(2*eig(matAl*matAl'));
%Calculating the 1/L (Lipschitz function)?(-1)
y=A';
x1=y;
x=A';
s=1;
iii = 1;

xlast = x;
for k=1:iterk
c=y-2*t*matAl'*(matAl*y-integration_corrected);
for l=1:size(c,1)
soft_thresholding(l)=(abs(c(l))-lambda*t)*sign(c(l));
end
not_zero=size(find(x));
x1=soft_thresholding';
s1=(1+sqrt(1+4*s^2))/2;

for i=1:size(x1,1)
if x1(i)<0
x(i)=0;
x1(i)=0;
end
end
y=x1+(s-1)/(s1)*(x1-x);

s=s1;

xlast = x;
x=(x1);
diffe(k) = norm(x-xlast)/norm(x);
if diffe(k)<1e-5
    break; 
end    
    
if mod(k,100) == 1
 
     X(:,iii) = x; % history of x in iteration
     iii = iii +1;
   
end
end
sprintf('End of ITAMeD procedure, non zero elements: %i',not_zero(1))
A_final=x;
end

