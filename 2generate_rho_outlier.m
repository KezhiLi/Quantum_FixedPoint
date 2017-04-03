function [outlier,X_true]=generate_rho_outlier(N,P,R,kk,tt)
%{
this fuction is used to generate the numerical density matrix of rank R 
and the outlier 
%}
LL=randn(N,R)+1i*randn(N,R);
RR=LL';
X_true=LL*RR;
X_true=X_true/trace(X_true);

KK=ceil(kk*N*P);%取整
loc=randperm(N*P);%最大值为NP的随机序列
loc=sort(loc(1:KK));%前KK个值排序
outlier=sparse(zeros(N*P,1));
outlier(loc) = sparse( tt * norm(X_true) * randn(KK,1) );
end