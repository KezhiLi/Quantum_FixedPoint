function [ A,y] = generate_A_y( AA,eta,X_true,outlier,N,P )
%{
    this function is used to generate A and y from AA according the measurment
    rate eta. It's a simulatoin of the process of measurement.
%}
M=ceil(eta*N*P);%就近取整
rows=randperm(N*P);
select_row=sort(rows(1:M));
A=AA(select_row,:);
y = A * ( reshape(X_true,N*P,1) + outlier);
end

