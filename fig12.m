N=30;
F=0;

%
%path(path,'/home/babblefish/UtilityPrograms/MatlabScripts/2019/matlab_bgl/');
path(path,'/home/babblefish/UtilityPrograms/MatlabScripts/2019/Scale free/');
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


gam=0.1;
alpha= 0.5+0.2*rand(N,1)+0.2*rand(N,1)-0.2*rand(N,1)-0.2*rand(N,1); 
recoveryRate=0.1+0.05*rand(N,1)+0.05*rand(N,1)-0.05*rand(N,1)-0.05*rand(N,1);

alpha= 0.5*ones(N,1); 
recoveryRate=0.1*ones(N,1);


Q = gam*F' + diag(alpha);
Q=[Q;recoveryRate'];

if(any(Q)<0)
    error('negative transition rates detected');
end


SimulationRuns=1500
observedTimes=-ones(N,SimulationRuns);
observedTimes(1,:)=0;

observedExpTimes=observedTimes;

startingNumber=4;

tic()
for(qqq=1:size(observedTimes,2))
    qqq
    P=startingNumber*eye(N,1);
    t=0;
    expt=0;

while(any(observedTimes(:,qqq)<0) & any(P>0) )
    eventRate= sum(sum( Q.*P'));
    t=t+ exprnd(1/eventRate);
    expt=expt+1./eventRate;
    
    select= rand()*eventRate;
    select= sum((cumsum(reshape(Q.*P',[(N*(N+1)),1]))<select))+1;
    [i,j]= deal(mod(select-1,(N+1))+1, floor((select-1)/(N+1))+1);
    
    if(i<=N)
        P(i)=P(i)+1;
    end
    
    if(i~=j)
       P(j)=P(j)-1; 
    end
    
    if(any(P<0))
        error('wut')
    end
        
    
    SetTimes = (observedTimes(:,qqq)<0 & P>0);
    observedTimes(SetTimes,qqq)=t;
    observedExpTimes(SetTimes,qqq)=expt;
end

    SetTimes=(observedTimes(:,qqq)<0);
    observedTimes(SetTimes,qqq)=inf;
    observedExpTimes(SetTimes,qqq)=inf;

end

simTime=toc()

bob=observedTimes(:);
bob=bob(bob<inf);
T= linspace(0,ceil(max(bob)*0.75),300);
clear bob

survivalCurve=zeros(N,length(T));

ExpectedsurvivalCurve=zeros(N,length(T));


for(iii=1:N)
    survivalCurve(iii,:)= mean((T<observedTimes(iii,:)'));
    ExpectedsurvivalCurve(iii,:)= mean((T<observedExpTimes(iii,:)'));
end

tic()
AnalyticSurvivalCurve=zeros(N,length(T));

for(iii=2:N)
    
y0=ones(N,1);
y0(iii)=0;
deriv= @(t,y) SurvivalDeriv(t,y,gam*F,alpha,iii,recoveryRate)
[~,ode45y]=ode45(deriv, T, y0,odeset('NonNegative',1));
AnalyticSurvivalCurve(iii,:)=ode45y(:,1)';
end

odeTime=toc();


figure()
subplot(1,2,1)
plot(T,survivalCurve);
ylabel('Observed Survival function');
xlabel('time');

subplot(1,2,2)
plot(T,AnalyticSurvivalCurve.^startingNumber);
ylabel('Analytic Survival function');
xlabel('time');

figure()
plot(T,AnalyticSurvivalCurve.^startingNumber-survivalCurve);
xlabel('time');
ylabel('analytic-simulation Survival');


odeTime
simTime