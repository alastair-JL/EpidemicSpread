%%This code encapsulates the weird graph used in Figure 4, runs our
%%simulations, and prints it all out to some figures.


gam=0.1;


alpha=[0,0.02,4,0.01,0,0,0]';
recoveryRate=[0,0.01,0,0,0,0,0]';
N=length(alpha);

F= [0 0 0 0 0 0 0;
   10 0 0 0 0 0 0;
   10 0 0 0 0 0 0;
   10 0 0 0 0 0 0;
    0 0.02 10^-0 10^-4 0 0 0;
    0 0.014 10^-4 10^-2 0 0 0;
    0 0.001 10^-2 10^-0 0 0 0]'

Q = gam*F' + diag(alpha);
Q=[Q;recoveryRate'];


startingNumber=1;


tic()
T=[0,logspace(-1,3,5000)];

AnalyticSurvivalCurve=zeros(N,length(T));


for(iii=2:N)
    
y0=ones(N,1);
y0(iii)=0;
deriv= @(t,y) SurvivalDeriv(t,y,gam*F,alpha,iii,recoveryRate)
[~,ode45y]=ode45(deriv, T, y0,odeset('NonNegative',1,'RelTol',10^-8));
AnalyticSurvivalCurve(iii,:)=ode45y(:,1)';
end

odeTime=toc();


figure()
subplot(1,2,2)
plot(T,AnalyticSurvivalCurve.^startingNumber);
ylabel('Analytic Survival function');
xlabel('time');
toc()


SimulationRuns=25000;
observedTimes=-ones(N,SimulationRuns);
observedTimes(1,:)=0;

observedExpTimes=observedTimes;

tic()
for(qqq=1:size(observedTimes,2))
    qqq
    P=startingNumber*eye(N,1);
    t=0;
    expt=0;

while(any(observedTimes(5:7,qqq)<0) & any(P(1:4)>0) )
    
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
    
    SetTimes = (observedTimes(:,qqq)<0 & P>0);
    observedTimes(SetTimes,qqq)=t;
    observedExpTimes(SetTimes,qqq)=expt;
        
end

    SetTimes=(observedTimes(:,qqq)<0);
    observedTimes(SetTimes,qqq)=inf;
    observedExpTimes(SetTimes,qqq)=inf;
    save('work.mat')
end

simTime=toc()

bob=observedTimes(:);
bob=bob(bob<inf);
%T= linspace(0,ceil(max(bob)*0.8),1300);
clear bob

survivalCurve=zeros(N,length(T));

for(iii=1:N)
    survivalCurve(iii,:)= mean((T<observedTimes(iii,:)'));
end
% 
% tic()
% AnalyticSurvivalCurve=zeros(N,length(T));
% 
% for(iii=2:N)
%     
% y0=ones(N,1);
% y0(iii)=0;
% deriv= @(t,y) SurvivalDeriv(t,y,gam*F,alpha',iii,recoveryRate)
% [~,ode45y]=ode45(deriv, T, y0,odeset('NonNegative',1));
% AnalyticSurvivalCurve(iii,:)=ode45y(:,1)';
% end
% 
% odeTime=toc();


% figure()
subplot(2,2,1)
plot(T,survivalCurve);
ylabel('Observed Survival function');
xlabel('time');
set(gca, 'XScale', 'log')

subplot(2,2,2)
plot(T,AnalyticSurvivalCurve.^startingNumber);
ylabel('Analytic Survival function');
xlabel('time');
set(gca, 'XScale', 'log')

subplot(2,2,3)
plot(T,AnalyticSurvivalCurve.^startingNumber-survivalCurve);
xlabel('time');
ylabel('Error');
set(gca, 'XScale', 'log')



subplot(2,2,4)


numSamples= [100,300,1000,3000,10000];

%numSamples= [100,300,1000];


errorNum=zeros(40,length(numSamples));

dts=  ([T(1:end-1),T(end-1)]-[0,T(1:end-1)])/2;


for(sss=1:length(numSamples))
    
    for(qqq=1:size(errorNum,1))
        sample = randsample(size(observedTimes,2),numSamples(sss));
        for(iii=1:N)
            survivalCurve(iii,:)= mean((T<observedTimes(iii,sample)'));
        end
    
    errorNum(qqq,sss)=sum(sum((AnalyticSurvivalCurve.^startingNumber-survivalCurve).^2.*dts));
    
    end
end
boxplot(errorNum);
