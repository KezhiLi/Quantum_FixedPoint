function [rho, S,result_3]=Robust_Quantum_ADMM3(y, A, flagS, R,tao2,lambda3,abs_tol,X_true)

[y_row,y_col]  = size(y);  
[A_row,A_col]  = size(A); 
d = round(sqrt(A_col)); %  the nearest integers 密度矩阵的维度d*d

S2=zeros(d,d);   
sq_error = inf;        
ti=0;
u = sparse(zeros(y_row,1));   
rho2 = A'*y;            
sq_before=0;
while(true)
    ti=ti+1;
    sq_error_old = sq_error;
    
    if ti>30  %%%  循环30次，若ti大于30则，不收敛，误差不满足约束    
        disp('[iteration failed]: doesnot meet the error bound')
        break
    end
    
    rho1 = reshape(inv(A'*A+0.001*eye(A_col))*A'*(y-u-A*reshape(S2,d*d,1)),d,d);  %迭代求rho  eye(A_col)：返回A_col*A_col的单位矩阵 
    rho2 = func_project_matrix(rho1, 6, R);
    
    Vec_S1 =  inv(A'*A+0.001*eye(A_col))*A'*(y-u-A*reshape(rho2,d*d,1));%迭代求S：A'*A+0.00001*eye(A_col)，矩阵加一个微小的矩阵
     S2 = sparse(real(reshape(sign(real(Vec_S1)).*max(abs(real(Vec_S1))-tao2, 0),d,d))); %reshape：重新把矩阵确定为32*32  sign：判断符号  tao2=0.01固定参数！！！
        
    u = sparse(u+ lambda3*(A*reshape(rho2+S2,d*d,1)-y));       %迭代求u  lambda3=0.3;
    
        %Compute residual
    %-----------------------------
    resid    = y - A*reshape( rho2+S2, d* d,1 );
    sq_error = resid'*resid;                %误差！！
     
     if  sq_error > sq_before
       lambda3=lambda3-0.01;
    else
        lambda3=lambda3+0.01;
     end
  
    if (sq_error < 0.1*abs_tol )  %|| (sq_error_old/sq_error < rel_tol)     %判断误差是否小于性能指标abs_tol
        sq_error = sq_error_old;   %判断误差是否小于性能指标abs_tol
        break; %terminate
    end
    
 %   
end
rho=rho2;
S=S2;
result_3=norm(X_true-rho,'fro')^2/norm(X_true,'fro')^2 ; %% F范数计算：sqrt(sum(diag(A'*A))).      
end

