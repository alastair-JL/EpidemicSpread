N=130;
F=0;


%path(path,'/home/babblefish/UtilityPrograms/MatlabScripts/2019/matlab_bgl/');
%path(path,'/home/babblefish/UtilityPrograms/MatlabScripts/2019/Scale free/');


while(any(sum(F))==0)
    F= rand(N).*(rand(N)<4/N);
    F=F-diag(diag(F));
    F=F+F';
end

gam=10^-2;
alpha= 0.5+0.2*rand(N,1)+0.2*rand(N,1)-0.2*rand(N,1)-0.2*rand(N,1); 
alpha=0*alpha+0.5
recoveryRate=0*alpha;
            %Okay, lets test with variable alpha
%alpha=0;%alpha/15;

maxPops= floor(5+exp(1.5+2*rand(N,1)))

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
diagIndices= (1:(N+1):N^2);

while(any(observedTimes(:,qqq)<0))
    EventBox=Q.*P';
    EventBox=EventBox.*((maxPops-P))./maxPops;
    CommonEvents = find(EventBox>inf);
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
    
    if(any(isnan(P)))
       error('Nans encountered.') 
    end
    
    SetTimes = (observedTimes(:,qqq)<0 & P>0);
    observedTimes(SetTimes,qqq)=t;
    observedExpTimes(SetTimes,qqq)=expt;
end


end



simTime=toc()

ObservedMeanArrivalTime=mean(observedTimes,2);

T= linspace(0,ceil(max(max(observedTimes)))*2.1,4300);

AnalyticSurvivalCurve=zeros(N,length(T));

for(iii=2:N) 
y0=ones(N,1);
y0(iii)=0;
deriv= @(t,y) SurvivalDeriv(t,y,gam*F,alpha,iii)
[~,ode45y]=ode113(deriv, T, y0,odeset('NonNegative',1));
AnalyticSurvivalCurve(iii,:)=ode45y(:,1)';
end

AnalyticArrivalTime= sum(AnalyticSurvivalCurve,2)*(T(2)-T(1));

md_brock = fitlm(AnalyticArrivalTime,ObservedMeanArrivalTime)
txt = ['R^2 = ', num2str(md_brock.Rsquared.Adjusted)];
scatter(AnalyticArrivalTime,ObservedMeanArrivalTime,15,D_dumb,'LineWidth',2);
hold on
plot([0,500],[0,500],'k:')
plot([0,500],[0,500]*md_brock.Coefficients.Estimate(2)+md_brock.Coefficients.Estimate(1),'r:')

text(0.15,0.15,txt,'Units','normalized')
xlabel('Predicted mean AT')
ylabel('Observed mean AT')

