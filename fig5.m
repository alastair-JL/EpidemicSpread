N=155;
F=0;

while(any(sum(F))==0)
    F= rand(N).*(rand(N)<4/N);
    F=F-diag(diag(F));
    F=F+F';
end

gam=10^-3;
alpha= 0.5+0.2*rand(N,1)+0.2*rand(N,1)-0.2*rand(N,1)-0.2*rand(N,1); 
alpha=0*alpha+0.5
recoveryRate=0*alpha;
            %Okay, lets test with variable alpha
%alpha=0;%alpha/15;

Q = gam*F' + diag(alpha);

SimulationRuns=5000
observedTimes=-ones(N,SimulationRuns);
observedTimes(1,:)=0;

observedExpTimes=observedTimes;


tic()
for(qqq=1:size(observedTimes,2))
    qqq
    P=eye(N,1);
    t=0;
    expt=0;
    
    
   StepMax=max(alpha)/200; %Set stepmax so that one in 200 things duplicates. Negligible chance of double duplication.
dt=0;

randset=rand(2,10000); 
randset(1,:)=-log(randset(1,:));
randIndx=0;

maxing=0;
stepping=0;

while(any(observedTimes(:,qqq)<0))
    EventBox=Q.*P';
    CommonEvents = find(EventBox>4);
    SpecialEvents=EventBox(CommonEvents);
    EventBox(CommonEvents)=0;
    
    eventRate= sum(sum(EventBox));
    randIndx=randIndx+1;
    if(randIndx>length(randset))
        randset=rand(2,10000); 
        randset(1,:)=-log(randset(1,:));
        randIndx=1;
    end
    
    DesiredStep=randset(1,randIndx)/eventRate;
    
    if(DesiredStep<StepMax)
        stepping=stepping+1;
        dt=DesiredStep;
        t=t+ dt;
    
        select= randset(2,randIndx)*eventRate;
        select= sum((cumsum(reshape(EventBox,[N^2,1]))<select))+1;
        [i,j]= deal(mod(select-1,N)+1, floor((select-1)/N)+1);
    
        P(i)=P(i)+1;
        if(i~=j)
            P(j)=P(j)-1; 
        end
    
    else
        maxing=maxing+1;
        dt=StepMax;
        t=t+ dt;
    end

    numEvents=poissrnd(dt*SpecialEvents);
            
    for(www=1:length(CommonEvents) )
        select=CommonEvents(www);
        [i,j]= deal(mod(select-1,N)+1, floor((select-1)/N)+1);

        P(i)=P(i)+numEvents(www);
        if(i~=j)
            P(j)=P(j)-numEvents(www); 
        end
    end
    
    
    SetTimes = (observedTimes(:,qqq)<0 & P>0);
    observedTimes(SetTimes,qqq)=t;
    observedExpTimes(SetTimes,qqq)=expt;
end


end

simTime=toc()

if(SimulationRuns>0)

T= linspace(0,ceil(max(max(observedTimes)))*1.1,1300);

survivalCurve=zeros(N,length(T));

ExpectedsurvivalCurve=zeros(N,length(T));


for(iii=1:N)
    survivalCurve(iii,:)= mean((T<observedTimes(iii,:)'));
   
end

else
    T= linspace(0,50,1300);
    survivalCurve=0;
   
end

tic()
AnalyticSurvivalCurve=zeros(N,length(T));

for(iii=2:N) 
y0=ones(N,1);
y0(iii)=0;
deriv= @(t,y) SurvivalDeriv(t,y,gam*F,alpha,iii)
[~,ode45y]=ode113(deriv, T, y0,odeset('NonNegative',1));
AnalyticSurvivalCurve(iii,:)=ode45y(:,1)';
end

odeTime=toc();

M= gam*F +diag(alpha) - gam*diag(sum(F,2));
M=M';

LogExpSurvivalCurve=zeros(N,length(T));

for(ttt=1:length(T) )
    LogExpSurvivalCurve(:,ttt)=1./(1+(eye(1,N)*expm(M*T(ttt)))');
end


muVals=zeros(N,1);
muVals2=zeros(N,1);

for(nnn=2:N)
muse= muCalcsFunc(gam*F,alpha(1),nnn);
muVals(nnn)=muse(1);

muse= muCalcsMethod2(gam*F,alpha(1),nnn);
muVals2(nnn)=muse(1);
end

MuseSurvivalCurve2= 1./(1+exp(alpha(1)*(T-muVals2)));
MuseSurvivalCurve2(1,:)=0;


MuseSurvivalCurve= 1./(1+exp(alpha(1)*(T-muVals)));

MuseSurvivalCurve(1,:)=0;
LogExpSurvivalCurve(1,:)=0;


subplot(3,3,1)
plot(T,survivalCurve);
ylabel('S(t)');
xlabel('time');
title('Gillespie Simulations')
xlim([T(1),T(end)])

subplot(3,3,2)
plot(T,AnalyticSurvivalCurve);
ylabel('S(t)');
xlabel('time');
title('Numerical solution ODE')
xlim([T(1),T(end)])

subplot(3,3,3)
plot(T,AnalyticSurvivalCurve-survivalCurve);
ylabel('error');
xlabel('time');
title('Numerical solution ODE')
xlim([T(1),T(end)])

subplot(3,3,5)
plot(T,MuseSurvivalCurve);
ylabel('S(t)');
xlabel('time');
title('Analytic \mu calculation')
xlim([T(1),T(end)])


subplot(3,3,6)
plot(T,MuseSurvivalCurve-survivalCurve);
ylabel('error');
xlabel('time');
title('Analytic \mu calculation')
xlim([T(1),T(end)])

subplot(3,3,8)
plot(T,LogExpSurvivalCurve);
ylabel('S(t)');
xlabel('time');
title('Matrix exponent approximation')
xlim([T(1),T(end)])

subplot(3,3,9)
plot(T,LogExpSurvivalCurve-survivalCurve);
ylabel('error');
xlabel('time');
title('Matrix exponent approximation')
xlim([T(1),T(end)])

subplot(3,3,[4,7])
plot(T,T ,'--k')
hold on
scatter(mean(observedTimes,2),sum(AnalyticSurvivalCurve,2)*(T(2)-T(1)),'o')
scatter(mean(observedTimes,2),sum(MuseSurvivalCurve,2)*(T(2)-T(1)),'*')
scatter(mean(observedTimes,2),sum(LogExpSurvivalCurve,2)*(T(2)-T(1)),'v')
ylabel('predicted');
xlabel('time');
xlim([T(1),T(end)])
ylim([T(1),T(end)])
