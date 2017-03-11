function [ output_args, eigv ,numiter,gpcrrnt] = fw( input_args, k,xini, datm )
% My implementation of the Frank-Wolfe alg, JJ Z
testiter=200;
gpcrrnt=zeros(testiter,1);
Q=input_args;
if nargin ==4
    A=datm;
end
[m n]=size(Q);
if nargin >2
    x=xini;
else
    x=zeros(n,1);
    x(1)=1;
j=1;
switch nargin
    case {2,3}
        qx=Q * x;
    case 4
        qx=(-A') * (A * x);
end
xqx=-1;
for j=1:testiter
    if j>2000
        break
    end
    x=proj(-qx,k);
    switch nargin
        case {2,3}
            qx=Q * x;
            
        case 4
            qx=(-A') * (A * x);
    end
    xqx=x'*qx;
    gpcrrnt(j,1)=-xqx;
end
output_args=x;
eigv=-xqx;
numiter=j;
end
