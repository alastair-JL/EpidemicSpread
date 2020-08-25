F=0;

%
%path(path,'/home/babblefish/UtilityPrograms/MatlabScripts/2019/matlab_bgl/');
%path(path,'/home/babblefish/UtilityPrograms/MatlabScripts/2019/Scale free/');
%
laticeSize=30;
N=laticeSize*laticeSize;

while(any(sum(F))==0)
    Vec=zeros(1,laticeSize)
    Vec(2)=1;
    Vec(end)=1;
    C = gallery('circul',Vec);
    F=kron(C,eye(laticeSize))+kron(eye(laticeSize),C);
end
F=F./sum(F);


gam=0.001;
alpha= 0.5+0.2*rand(N,1)+0.2*rand(N,1)-0.2*rand(N,1)-0.2*rand(N,1); 
recoveryRate=0.1+0.05*rand(N,1)+0.05*rand(N,1)-0.05*rand(N,1)-0.05*rand(N,1);

alpha= 0.5*ones(N,1); 
recoveryRate=0.1*ones(N,1);


Q = gam*F' + diag(alpha);
Q=[Q;recoveryRate'];

if(any(Q)<0)
    error('negative transition rates detected');
end


T= linspace(0,450,2000);
clear bob

tic()
AnalyticSurvivalCurve=zeros(N,length(T));

 
    iii=1
y0=ones(N,1);
y0(iii)=0;
deriv= @(t,y) SurvivalDeriv(t,y,gam*F,alpha,iii,recoveryRate);
[~,ode45y]=ode45(deriv, T, y0,odeset('NonNegative',1));


odeTime=toc();

figure()
subplot(1,2,2)
plot(T,ode45y(:,:));
ylabel('S(t)');
xlabel('time');

odeTime
