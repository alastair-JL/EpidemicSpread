%%This is the code used to make Figure 2 of the paper. Here we consider a
%%Simple chain graph.


gam=0.1;

N=8;
alpha=0.5*ones(N,1);
recoveryRate=zeros(N,1);
N=length(alpha);

F= diag(ones(1,N-1),1);

Q = gam*F' + diag(alpha);
Q=[Q;recoveryRate'];


startingNumber=1;


tic()
T=linspace(0,40);

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
plot(T,AnalyticSurvivalCurve.^startingNumber);
ylabel('Analytic Survival function');
xlabel('time');
toc()


SimulationRuns=5000;
observedTimes=-ones(N,SimulationRuns);
observedTimes(1,:)=0;

observedExpTimes=observedTimes;

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
    
    SetTimes = (observedTimes(:,qqq)<0 & P>0);
    observedTimes(SetTimes,qqq)=t;
    observedExpTimes(SetTimes,qqq)=expt;
end

    SetTimes=(observedTimes(:,qqq)<0);
    observedTimes(SetTimes,qqq)=inf;
    observedExpTimes(SetTimes,qqq)=inf;

end

simTime=toc();

subplot(1,2,2)

for(iii=1:N)
scatter(sum(AnalyticSurvivalCurve(iii,:))*T(2) ,mean(observedTimes(iii,:)), 'LineWidth',2);
hold on    
end

ylabel('Observed Mean AT');
xlabel('Predicted mean AT');