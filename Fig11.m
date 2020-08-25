N=20;
F=0;

while(any(sum(F))==0)
    F= rand(N).*(rand(N)<5/N);
    F=F-diag(diag(F));
    F=F+F';
end
F=F./sum(F);

gam=10^-1.2;
layers=5;

alphaBase= 0.5+0.2*rand(N,1)+0.2*rand(N,1)-0.2*rand(N,1)-0.2*rand(N,1);

alphaMult=[1,0.2,1.5,0.2,1.5,0.2,1.5,0.2,1.5,0.2,1.5,0.2,1.5,0.2,1.5,0.2,1.5];
alpha=alphaBase*alphaMult(1);

            %Okay, lets test with variable alpha
%alpha=0;%alpha/15;



            
Q = gam*F' + diag(alpha);

SimulationRuns=120
observedTimes=-ones(N,SimulationRuns);
observedTimes(1,:)=0;

observedExpTimes=observedTimes;


TimeThresholds= [0,5,10,15,20,25,30,35,40,45,inf];
nextThresholdIndex=2;


tic()
for(qqq=1:size(observedTimes,2))
    qqq
    P=eye(N,1);
    t=0;
    expt=0;
    nextThresholdIndex=1;

while(any(observedTimes(:,qqq)<0))
    eventRate= sum(sum( Q.*P'));
    t=t+ exprnd(1/eventRate);
    
    if(t>TimeThresholds(nextThresholdIndex))
        t=TimeThresholds(nextThresholdIndex);
        alpha=alphaBase*alphaMult(nextThresholdIndex);
        Q = gam*F' + diag(alpha);
        nextThresholdIndex=nextThresholdIndex+1;        
        continue;
    end
    
    
    expt=expt+1./eventRate;
    
    select= rand()*eventRate;
    select= sum((cumsum(reshape(Q.*P',[N^2,1]))<select))+1;
    [i,j]= deal(mod(select-1,N)+1, floor((select-1)/N)+1);
    
    P(i)=P(i)+1;
    if(i~=j)
       P(j)=P(j)-1; 
    end
    
    SetTimes = (observedTimes(:,qqq)<0 & P>0);
    observedTimes(SetTimes,qqq)=t;
    observedExpTimes(SetTimes,qqq)=expt;
end


end

simTime=toc()

T= linspace(0,ceil(max(max(observedTimes)))*0.6,300);

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
        AnalyticSurvivalCurve(iii,1)=y0(1);
        for(www=2:length(T))
            [www,iii]
            deriv= @(t,y) SurvivalDerivPiecewise(T(www)-t,y,gam*F,alphaBase,iii,0*alpha,alphaMult,TimeThresholds);
            [~,ode45y]=ode45(deriv, [0,T(www)], y0,odeset('NonNegative',1,'Reltol',1e-13,'AbsTol',1e-14));
            AnalyticSurvivalCurve(iii,www)=ode45y(end,1)';
        
        end
    end


odeTime=toc();




figure()
subplot(1,2,1)
plot(T,survivalCurve);
hold on
for(rrr=1:length(TimeThresholds))
plot([TimeThresholds(rrr),TimeThresholds(rrr)],[0,1],':','color',[0,0,0]+0.5);
end
ylabel('Observed Survival function');
xlabel('time');



subplot(1,2,2)
plot(T,AnalyticSurvivalCurve);
ylabel('Analytic Survival function');
xlabel('time');
hold on
for(rrr=1:length(TimeThresholds))
plot([TimeThresholds(rrr),TimeThresholds(rrr)],[0,1],':','color',[0,0,0]+0.5);
end


figure()
plot(T,AnalyticSurvivalCurve-survivalCurve);
xlabel('time');
ylabel('analytic-simulation Survival');



tic()
BadAnalyticSurvivalCurve=zeros(N,length(T));

    for(iii=2:N)   
        y0=ones(N,1);
        y0(iii)=0;
         deriv= @(t,y) SurvivalDerivPiecewise(t,y,gam*F,alphaBase,iii,0*alpha,alphaMult,TimeThresholds);
         [~,ode45y]=ode45(deriv, T, y0,odeset('NonNegative',1,'Reltol',1e-13,'AbsTol',1e-14));
         BadAnalyticSurvivalCurve(iii,:)=ode45y(:,1)';
    end

figure()
subplot(2,2,1)
plot(T,survivalCurve);
xlim([0,25]);
hold on
for(rrr=1:length(TimeThresholds))
plot([TimeThresholds(rrr),TimeThresholds(rrr)],[0,1],':','color',[0,0,0]+0.5);
end
ylabel('Observed Survival function');
xlabel('time');



subplot(2,2,3)
plot(T,AnalyticSurvivalCurve);
ylabel('Analytic Survival function (\alpha(\tau-t))');
xlabel('time');
xlim([0,25]);
hold on
for(rrr=1:length(TimeThresholds))
plot([TimeThresholds(rrr),TimeThresholds(rrr)],[0,1],':','color',[0,0,0]+0.5);
end
    
    
subplot(2,2,2)
plot(T,BadAnalyticSurvivalCurve);
xlim([0,25]);
ylabel('Analytic Survival function (\alpha(t))');
xlabel('time');
hold on
for(rrr=1:length(TimeThresholds))
plot([TimeThresholds(rrr),TimeThresholds(rrr)],[0,1],':','color',[0,0,0]+0.5);
end
    
subplot(2,2,4)
for(rrr=1:length(TimeThresholds))
plot([0,25],[TimeThresholds(rrr),TimeThresholds(rrr)],':','color',[0,0,0]+0.5);
end
hold on
patch([0,25,0],[0,25,25],[0.75,0.75,0.75])
plot([0,25],[0,25],'k','linewidth',3);

ylim([0,25]);
xlim([0,25]);
ylabel('\tau');
xlabel('time');


odeTime
simTime