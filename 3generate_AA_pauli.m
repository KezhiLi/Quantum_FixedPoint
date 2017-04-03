function [ AA ] = generate_AA_pauli( N,P ,qubit)
%{
    this fuction is used to generate the complete measuement set. 
    the measurement matrix is a subset of AA
%}
 A1=[1,0;0,1];
    A2=[0,1;1,0];
    A3=[0,-1i;1i,0];
    A4=[1,0;0,-1];
    A0=[A1;A2;A3;A4];
  
    AA=[];
    for qq=0:(N*P)-1;
    
    Qun_ran=zeros(qubit,1);
    str1=dec2base(qq, 4);
    for ii=1:length(str1)
    Qun_ran(qubit-ii+1)=str2num(str1(length(str1)-ii+1));
    end

    A_block = [A0(Qun_ran(1)*2+1,:);
           A0(Qun_ran(1)*2+2,:)];
    for ii=2:qubit;
        A_block= kron(A_block,[A0(Qun_ran(ii)*2+1,:); A0(Qun_ran(ii)*2+2,:)]);
    end
    AA=[AA; reshape(A_block,1,N*P)];
    end
AA=AA/N^0.5;
end

