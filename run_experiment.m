% Code to generate the plot in the paper https://arxiv.org/abs/1404.4132
% JJ Z
clear all
load('seed1.mat');
rng(s);
zz=0;
sumds=0;
sumgp=0;
sumfw=0;
rsumds=0;
rsumgp=0;
rsumfw=0;
accu_gpbb=zeros(200,1);
accu_gpbbls=zeros(200,1);
accu_fw=zeros(200,1);
for irun=1:10 ;
    % prep data
    m=250; n=500;
    A=normrnd(0,1,m,n);
    % centering
    A=A-repmat((mean(A,1)),m,1);
    Q=-(A')*A;
    % initialization
    [pc, D]=eig(Q);
    id=eye(size(Q));
    card = 500 ;
    ratio1=zeros(n,1);
    ratio2=ratio1;
    % algorithm 1, GPBB
    V=zeros(n,n);
    vari=zeros(n,1);
    t=cputime;
    for k=1:1
        Qk=Q;
        for i=1:1
            [V(:, i), vari(i),~,gpbb_current]=gpbb(Qk,card); 
        end
        ratio1(card)=vari(1);
    end
    t1=cputime-t;
    % algorithm 2, Frank-Wolfe
    V=zeros(n,n);
    vari=zeros(n,1);  
    t=cputime;
    for k=1:1
        Qk=Q;
        for i=1:1
            [V(:, i), vari(i),~,fw_current]=fw(Qk,card);   % Use Frank-Wolfe alg
        end
        ratio2(card)=vari(1);
    end
    t2=cputime-t;
    % algorithm 3, GPBB w/ non-monotone line search
    V=zeros(n,n);
    vari=zeros(n,1);
    t=cputime;
    for k=1:1
        Qk=Q;
        for i=1:1
            [V(:, i), vari(i),~,gpbb_ls_current]=gpbbls(Qk,card);
        end
    end
    t3=cputime-t;
    % Take averages
    accu_gpbb=accu_gpbb+gpbb_current/(-D(1));
    accu_fw=accu_fw+fw_current/(-D(1));
    accu_gpbbls=accu_gpbbls+gpbb_ls_current/(-D(1));
end
accu_gpbb = abs(1 - accu_gpbb/irun) ;
I = find (accu_gpbb == 0) ;
accu_gpbb(I) = 1.e-16 ;
plot1=plot(log10(accu_gpbb),'-o');
hold on
plot2=plot(log10(1-accu_fw/irun),'-+r');
accu_gpbbls = abs(1 - accu_gpbbls/irun) ;
I = find (accu_gpbbls == 0) ;
accu_gpbbls(I) = 1.e-16 ;
plot3=plot(log10(accu_gpbbls),'-sk');
hold off
title('Random data, m=250, n=500, cardinality=500.','FontSize',16)
xlabel('Iteration Number','FontSize',16)
ylabel('LOG10 (Error)','FontSize',16)
leg1=legend([plot1;plot2;plot3],'GPBB','Frank-Wolfe','GPBB with line search');
set(leg1,'Location','SouthWest')
set(gca,'Fontsize',16)