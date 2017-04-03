function [ AA ] = generate_AA_pauli( N,P ,qubit)
% m1=[sqrt(3)/3,sqrt(3)/3,sqrt(3)/3];
% m2=[-sqrt(3)/3,sqrt(3)/3,sqrt(3)/3];
% m3=[-sqrt(3)/3,-sqrt(3)/3,sqrt(3)/3];
% m4=[sqrt(3)/3,-sqrt(3)/3,sqrt(3)/3];
m1=[-sqrt(6)/3,-sqrt(2)/3,-1/3];
m2=[sqrt(6)/3,-sqrt(2)/3,-1/3];
m3=[0,sqrt(8)/3,-1/3];
m4=[0,0,1];
sita1=[1,0;0,1];
sita2=[0,1;1,0];
sita3=[0,-1i;1i,0];
sita4=[1,0;0,-1];
A0=(m1(1)*sita2+m1(2)*sita3+m1(3)*sita4 +sita1)/2;
A1=(m2(1)*sita2+m2(2)*sita3+m2(3)*sita4 +sita1)/2;
A2=(m3(1)*sita2+m3(2)*sita3+m3(3)*sita4 +sita1)/2;
A3=(m4(1)*sita2+m4(2)*sita3+m4(3)*sita4 +sita1)/2;
A4=[A0;A1;A2;A3];

AA=rand(N*P,N*P);
AA=complex(AA);

for i=0:N*P-1  %AAµÄµÚ0 1 2 ....N*P-1ÐÐ
    str=dec2bin(i,2*qubit);
    str_4=zeros(1,qubit);
    
    Kronecker=1;
    for j=1:qubit
        str_4(j)=str2num(str(2*j))+2*str2num(str(2*j-1));
        Kronecker=kron(Kronecker,A4(2*str_4(j)+1:2*str_4(j)+2,:));
    end
    temp=reshape(Kronecker,1,N*P);
    AA(i+1,:)=temp;
end
AA=AA / sqrt(N );
end

