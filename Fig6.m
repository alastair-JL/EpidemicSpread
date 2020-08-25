N=2;
F=0;
T=linspace(0,15)

F=ones(N,N)
F=F-diag(diag(F));
T=linspace(0,15)


gam=0.000003;
            %Okay, lets test with variable alpha
alpha=[1;1];
recoveryRate=0.0;
            
Q = gam*F' + diag(alpha);
Q=[Q;recoveryRate*ones(1,N)];


startingNumber=1;


observedTimes=-ones(N,SimulationRuns);
observedTimes(1,:)=0;

observedExpTimes=observedTimes;


Tsmooth=linspace(0,15,1000)

tic()
AnalyticSurvivalCurve=zeros(N,length(Tsmooth));

for(iii=2:N)    
y0=ones(N,1);
y0(iii)=0;
deriv= @(t,y) SurvivalDeriv(t,y,gam*F,alpha,iii,recoveryRate)
[~,ode45y]=ode45(deriv, Tsmooth, y0,odeset('NonNegative',1));
AnalyticSurvivalCurve(iii,:)=ode45y(:,1)';
end

odeTime=toc();


for(iii=2:N)    
y0=ones(N,1);
y0(iii)=0;
deriv= @(t,y) SurvivalDeriv(t,y,gam*F,alpha,iii,recoveryRate)
[~,ode45y]=ode45(deriv, Tsmooth, y0,odeset('NonNegative',1));
AnalyticSurvivalCurve(iii,:)=ode45y(:,1)';
end

odeTime=toc();


AnalyticSurvivalCurve=AnalyticSurvivalCurve.^50;

gumbelApproximation= exp(-exp(((Tsmooth+4)*(alpha(1)-recoveryRate(1))+log(gam))))

figure()
subplot(2,1,1)
ylabel('S(t)');
xlabel('time');
hold on
plot(Tsmooth,gumbelApproximation,':','LineWidth',2);
plot(Tsmooth,AnalyticSurvivalCurve(2,:),'LineWidth',2);
legend({'Gumbel','Logistic (p=50)'})

subplot(2,1,2)
ylabel('P.D.F');
xlabel('time');
hold on
plot(Tsmooth(2:end)-Tsmooth(2)/2,diff(1-gumbelApproximation)./(Tsmooth(2)-Tsmooth(1)),':','LineWidth',2);
plot(Tsmooth(2:end)-Tsmooth(2)/2,diff(1-AnalyticSurvivalCurve(2,:))./(Tsmooth(2)-Tsmooth(1)) ,'LineWidth',2);
