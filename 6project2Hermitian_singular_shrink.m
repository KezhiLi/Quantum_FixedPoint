function [ rho2 ] = project2Hermitian_singular_shrink( rho1 ,tao)
%{
    this function is used to project the matrix rho1 to a low rank and
    Conjugate symmetric matrix rho2.
%}
[m,n]=size(rho1);
d=sqrt(m);
rho1=reshape(rho1,d,d);
 X=(rho1+rho1')/2;
 [V,D]=eig(X);
 diagD=diag(D);
 diagDD=max(diagD-tao,0);
 diagDT=diagDD+min(diagD+tao,0);
 rho2=V*diag(diagDT)*V';
 
end

