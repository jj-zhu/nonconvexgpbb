function [ output_args, eigv ,numiter,gpcrrnt] = gpbbls( input_args, k,xini, datm )
% Non-convex gradient projection with BB approx., w/ non-monotone line search
% as in the paper. JJ Z
maxiter=1000;
testiter=200;
gpcrrnt=zeros(testiter,1);
Q=input_args;
ai=1;
if nargin == 4
    A=datm;
end
[m n]=size(Q);
if nargin >2
    x=xini;
else
    x=zeros(n,1);
    x(1)=1;
end
oldx=zeros(n,1);
switch nargin
    case {2,3}
        qx=Q * x;
        gradfx=2* qx;
        gradfxold=gradfx;
    case {4,5}
        qx=(-A') * (A * x);
        gradfx=2* qx;
        gradfxold=gradfx;
end
j=1;
xqx=1;
for j=1:testiter
    if j>maxiter
        break
    end
    if j>2
        ai=((x-oldx)'*(gradfx-gradfxold))/((x-oldx)'*(x-oldx));
    end
    oldx=x;
    gradfxold=gradfx;
    if gradfx==0
        break
    end
    i_search=0;
    while i_search==0 || (j>2 && (-xqx < gpcrrnt(j-1,1)) )
        i_search=i_search+1;
        if i_search >30
            break
        end
        x=proj(ai*x - gradfx,k);
        switch nargin
            case {2,3}
                qx=Q*x;
            case {4,5}
                qx=(-A') * (A * x);
        end
        xqx=x'*qx;
        ai=.25*ai;
    end
    gradfx=2* qx;
    gpcrrnt(j,1)=-xqx;
end
output_args=x;
eigv=-xqx;
numiter=j;
end
