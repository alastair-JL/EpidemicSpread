N=135;
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

d_jk= (F+eye(N,N)>0)*1.00;
d_jk=sparse(d_jk);
D_dumb= -ones(N,1);
sss=0;
while(any(D_dumb<0));
    Arrived=((d_jk^sss*eye(N,1))>0);
    D_dumb(Arrived&D_dumb<0)=sss;
    sss=sss+1
    
    if(sss>N)
        error('Disconnected Graph');
    end
end

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
    CommonEvents = find(EventBox>10*StepMax);
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

ObservedMeanArrivalTime=mean(observedTimes,2);

T= linspace(0,ceil(max(max(observedTimes)))*2.1,4300);


alphaPerturb=[0.2,1,5,0.2,1,5,0.2,1,5];
gammaPerturb=[0.2,0.2,0.2,1,1,1,5,5,5];


alphaPerturb=[0.1,1,10,0.1,1,10,0.1,1,10];
gammaPerturb=[0.1,0.1,0.1,1,1,1,10,10,10];


for(jjj=1:9)

T= linspace(0,ceil(max(max(observedTimes)))*2.1/alphaPerturb(jjj),4300)
AnalyticSurvivalCurve=zeros(N,length(T));

for(iii=2:N) 
y0=ones(N,1);
y0(iii)=0;
deriv= @(t,y) SurvivalDeriv(t,y,gammaPerturb(jjj)*gam*F,alphaPerturb(jjj)*alpha,iii)
[~,ode45y]=ode113(deriv, T, y0,odeset('NonNegative',1));
AnalyticSurvivalCurve(iii,:)=ode45y(:,1)';
end

AnalyticArrivalTime= sum(AnalyticSurvivalCurve,2)*(T(2)-T(1));

figure(1)
subplot( 3 , 3 , jjj ) 
md_brock = fitlm(AnalyticArrivalTime,ObservedMeanArrivalTime)
txt = ['R^2 = ', num2str(md_brock.Rsquared.Adjusted)];
oneone=[0,min(max(ExpFakeBreakoutTimes(keepOnlyValid)),max(TimeBreakOut(keepOnlyValid)))];
plot(oneone,oneone,'--','color',[0.7,0.7,0.7]);
hold on
FitLine=[0,max(ExpFakeBreakoutTimes(keepOnlyValid))];
plot(FitLine,md_brock.Coefficients.Estimate(1)+md_brock.Coefficients.Estimate(2)*FitLine,':k');

scatter(AnalyticArrivalTime,ObservedMeanArrivalTime,15,D_dumb,'LineWidth',2);

text(0.2,0.15,txt,'Units','normalized')
txt = ['\gamma x', num2str(gammaPerturb(jjj))];
text(0.15,0.85,txt,'Units','normalized')
txt = ['\alpha x', num2str(alphaPerturb(jjj))];
text(0.15,0.70,txt,'Units','normalized')
txt = ['\tau', num2str(corr(AnalyticArrivalTime,ObservedMeanArrivalTime,'Type','Kendall'))];
text(0.15,0.55,txt,'Units','normalized')


xlabel('Predicted AT')
ylabel('Observed AT')


end



xiVal=[0,0.1,1];

for(jjj=1:3)

T= linspace(0,ceil(max(max(observedTimes)))*2.1,4300)
AnalyticSurvivalCurve=zeros(N,length(T));

Fpert=F.*exp(sqrt(xiVal(jjj))*randn(N,N));
Fpert=(Fpert+Fpert')/2;

for(iii=2:N) 
y0=ones(N,1);
y0(iii)=0;
deriv= @(t,y) SurvivalDeriv(t,y,gam*Fpert,alpha,iii)
[~,ode45y]=ode113(deriv, T, y0,odeset('NonNegative',1));
AnalyticSurvivalCurve(iii,:)=ode45y(:,1)';
end

AnalyticArrivalTime= sum(AnalyticSurvivalCurve,2)*(T(2)-T(1));

figure(2)
subplot( 1 , 3 , jjj ) 
md_brock = fitlm(AnalyticArrivalTime,ObservedMeanArrivalTime)
txt = ['R^2 = ', num2str(md_brock.Rsquared.Adjusted)];
scatter(AnalyticArrivalTime,ObservedMeanArrivalTime,15,D_dumb,'LineWidth',2);
text(0.15,0.85,txt,'Units','normalized')
txt = ['\xi =', num2str(xiVal(jjj))];
text(0.1,0.75,txt,'Units','normalized')
txt = ['\tau', num2str(corr(AnalyticArrivalTime,ObservedMeanArrivalTime,'Type','Kendall'))];
text(0.15,0.55,txt,'Units','normalized')

xlabel('Predicted Mean AT')
ylabel('Mean Observed AT')

figure(3)
subplot( 1 , 3 , jjj ) 
md_brock = fitlm(AnalyticArrivalTime,observedTimes(:,1))
txt = ['R^2 = ', num2str(md_brock.Rsquared.Adjusted)];
scatter(AnalyticArrivalTime,observedTimes(:,1),15,D_dumb,'LineWidth',2);
text(0.10,0.85,txt,'Units','normalized')
txt = ['\xi =', num2str(xiVal(jjj))];
text(0.1,0.75,txt,'Units','normalized')
xlabel('Predicted Mean AT')
ylabel('Mean Observed AT')
end
