function [ rho,S,result_3 ] =Robust_Quantum_fixedpoint ( y,A,maxite,tol_1 ,X_true)
%{
   this is the funtion of reconstructe the density matrix by FP_ADMM algorithm. X_hat is the 
  reconstructed matrix, S is the reconstructed outlier. If the simulation without outlier, set
  it to zero.
%}
[y_row,y_col]  = size(y);
[A_row,A_col]  = size(A); 
d = round(sqrt(A_col)); %  the nearest integers 
lambda=1/sqrt(d); %can be tuned
rho2=zeros(d,d); %rho2的初始值
S1=zeros(d*d,1);
 
u=0.5/norm(y,'fro'); % u的初始值 can be tuned
sita=1.5; %can be tuned

Y = zeros(y_row, 1);
u_bar=u*1e7;
Temp=eye(A_col)-A'*A; %d^2*d^2
converged=0;
ite=0;

while ~converged
    ite=ite+1;
    
    if mod(ite,10)==0
        disp(['iterations:' num2str(ite)]);
    end
    
    S2=S1;
    temp=Temp*reshape(S1,d*d,1)+A'*(y-Y/u-A*reshape(rho2,d*d,1));
    temp=abs(temp);
    S1=max(temp-lambda/u,0);
    S1=S1+min(temp+lambda/u,0); % d^2*1
  
    
    rho1=Temp*reshape(rho2,d*d,1)+A'*(y-Y/u-A*reshape(S1,d*d,1)); % d^2*1
    rho2=project2Hermitian_singular_shrink(rho1,1/u); % d*d
    rho3=reshape(rho2,d*d,1);
   
    resid=A*(rho3+S1)-y;
    Y=Y+u*resid;

%      u=min(sita*u, u_bar);% 有干扰时建议u用这种更新方式，u不变结果不收敛
%      if u * norm(S1-S2, 'fro') / norm(y, 'fro') < 1e-5
%          u = 1.05* u;
%      else
%          u = 0.75* u;
%      end
     
    stop=norm(resid,'fro') / norm(y, 'fro');
    
    if stop<tol_1
        converged=1;
    end
     
    if ~converged && ite>maxite
        disp('maximum iteration reached');
        converged=1;
    end
    
end

S=reshape(S1,d,d);
rho=rho2;
result_3=norm(X_true-rho,'fro')^2/norm(X_true,'fro')^2;

end

