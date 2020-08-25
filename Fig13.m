N=3000;
F=0;

%
%path(path,'/home/babblefish/UtilityPrograms/MatlabScripts/2019/matlab_bgl/');
%path(path,'/home/babblefish/UtilityPrograms/MatlabScripts/2019/Scale free/');
%

while(any(sum(F))==0)
    seed =[0 1 0 0 1;1 0 0 1 0;0 0 0 1 0;0 1 1 0 0;1 0 0 0 0];
Net = double(SFNG(N, 5, seed));
Nk=(sum(Net).*10000.*exp(0.1*randn(1,N)))'; %Okay, so this is a kind of derpy interpretation of size...
%Nk=(1000.*exp(0.1*randn(1,NumCities)))'; %Okay, so this is a kind of derpy interpretation of size...

Net=Net.*Nk; %% This increases our rate of flying TO large airports.
Net=Net./sum(Net).*exp(0.1*randn(1,N)); %% and here we make things normalISH.
%PL_Equation = PLplot(round(Net))
F = sparse(double(Net));
end
F=F./sum(F);


gam=0.001;
alpha= 0.5+0.2*rand(N,1)+0.2*rand(N,1)-0.2*rand(N,1)-0.2*rand(N,1); 
recoveryRate=0.1+0.05*rand(N,1)+0.05*rand(N,1)-0.05*rand(N,1)-0.05*rand(N,1);

alpha= 0.5*ones(N,1); 
recoveryRate=0.1*ones(N,1);
gam=10^-5;

Q = gam*F' + diag(alpha);
Q=[Q;recoveryRate'];

if(any(Q)<0)
    error('negative transition rates detected');
end


T= linspace(0,150,1000);
clear bob

tic()
AnalyticSurvivalCurve=zeros(N,length(T));

for(iii=2:N)    
    iii
y0=ones(N,1);
y0(iii)=0;
deriv= @(t,y) SurvivalDeriv(t,y,gam*F,alpha,iii,recoveryRate);
[~,ode45y]=ode45(deriv, T, y0,odeset('NonNegative',1));
AnalyticSurvivalCurve(iii,:)=ode45y(:,1)';
end

odeTime=toc();

figure()
plot(T,AnalyticSurvivalCurve);
ylabel('Analytic Survival function');
xlabel('time');

odeTime
