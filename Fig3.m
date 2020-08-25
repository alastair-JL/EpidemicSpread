%%This is the code used to make Figure 3 of the paper. Here we consider a
%%"spinning top" graph in which one node is accesible via a single path,
%%and another node is accesible via many "narrower" paths.



% This code makes use of the "matlab_bgl" package by David Gleich for
% shortest distance calculation. In order to use this code, please download the package from 
% https://dgleich.github.io/matlab-bgl/ and replace the  following line with the correct path.
%%%
path(path,'/home/babblefish/UtilityPrograms/MatlabScripts/2019/matlab_bgl/');
%%%


tic();

k=80;

NumCities=4+k;
WidthWindow=1200;
Threshold=1;
repetition=150;
dt=0.01;


Nk= ones(NumCities,1)*80000000;

Network = speye(NumCities,NumCities);
Network(1,2)=1;
Network(2,3)=1;
Network(4:end-1,3)=1/k;
Network(end,4:end)=1;
%% This makes A network. Not nessesarily the most amazing one, but hey, its a network, it does its thing.

%Network= exp(rand(NumCities,NumCities)*4); %Lets claim that this is the number of people transfered each day...
%maxs= maxk(Network,5);
%Network= (Network.* (Network> maxs(end,:)));


[i,j,links] = find(Network)


%Network=(Network+Network');%%Make Symetic, yay.
Network=Network-diag(diag(Network));
%Network=Network.*Nk';
%Network=Network./sum(Network).*exp(rand(1,NumCities)*3);
setUpTime=toc();
'Set Up complete'

tic();
Ptrans=Network./sum(Network);
d_jk= 1 - log(Ptrans);
d_jk(d_jk==inf)=0;
d_jk=sparse(d_jk);
[D_brock,zob] = all_shortest_paths(d_jk,struct('algname','floyd_warshall'));
D_brock=D_brock';

d_jk= (Ptrans>0)*1.00;
d_jk=sparse(d_jk);
[D_dumb,zob] = all_shortest_paths(d_jk,struct('algname','floyd_warshall'));
D_dumb=D_dumb'

basicDistanceTime=toc();
 'have made short paths'

alpha=0.5;
gam=0.01;

AdjustW_jk=sparse(gam*(Network-diag(sum(Network))));

startNode=3;

HistoryInfectionLevel= zeros(NumCities,WidthWindow);
HistoryInfectionLevel(startNode,1)=1;

PowerFaker=HistoryInfectionLevel;

RecoveredLevel=zeros(NumCities,1);

recoveryRate=0;


tic()
T=linspace(0,40);

AnalyticSurvivalCurve=zeros(k+4,length(T));

FunctionFilter=[1:k+4]
for(iii=FunctionFilter)    
y0=ones(size(Network,1),1);
y0(iii)=0;
deriv= @(t,y) SurvivalDeriv(t,y,gam*Network',alpha,iii,recoveryRate)
[~,ode45y]=ode45(deriv, T, y0,odeset('NonNegative',1));
AnalyticSurvivalCurve(iii,:)=ode45y(:,startNode)';
end
AnalyticSurvivalCurve=AnalyticSurvivalCurve.^2;

AnalyticSurvivalCurveConditioned=(AnalyticSurvivalCurve-AnalyticSurvivalCurve(:,end))./(AnalyticSurvivalCurve(:,1)-AnalyticSurvivalCurve(:,end))
ExpFakeBreakoutTimes=sum(AnalyticSurvivalCurveConditioned')*T(2)


tic();

TimeBreakOuts=zeros(NumCities,repetition);

for (rrr=1:repetition)
    rrr

RecoveredLevel=zeros(NumCities,1);
InfectedLevel=zeros(NumCities,1);
InfectedLevel(startNode)=2;

TimeBreakOut=-ones(NumCities,1);
PopBreakOut=-ones(NumCities,1);
time=0;


tic();
ttt=0;
while(true)
    
    ttt=ttt+1;
    time=time+dt;
    
    nowRecov= min(poissrnd(recoveryRate*dt*InfectedLevel),InfectedLevel);
    nowSick= poissrnd(alpha.*dt.*InfectedLevel.* (Nk-InfectedLevel-RecoveredLevel)./Nk);
    nowMoveForwardSeed= poissrnd(gam*links.*dt.*InfectedLevel(i) );
    nowMoveReverse= poissrnd(gam*links.*dt.*InfectedLevel(j) );
    motionMatrix= sparse(i,j,nowMoveForwardSeed)+sparse(j,i,nowMoveReverse);
    moveGain=sum(motionMatrix,1)';
    moveLoss=sum(motionMatrix,2);
    InfectedLevel= max(InfectedLevel+ moveGain-moveLoss+nowSick-nowRecov,0);
    RecoveredLevel=RecoveredLevel+nowRecov;
    
    NowOutbreak= (TimeBreakOut==-1) & (InfectedLevel>=Threshold);
    TimeBreakOut(NowOutbreak)=time-dt;
    
    
    if(all( TimeBreakOut>=0) )
        break;
    end
    
    
    if(all(InfectedLevel==0) ) %%Simulation failure. Try again.
        InfectedLevel=zeros(NumCities,1);
        InfectedLevel(startNode)=1;
        TimeBreakOut=-ones(NumCities,1);
        time=0;
        ttt=0;
        'fail'
    end
    
    if(any( InfectedLevel<0) ) %%Simulation failure. Try again.
        InfectedLevel=zeros(NumCities,1);
        InfectedLevel(startNode)=1;
        TimeBreakOut=-ones(NumCities,1);
        time=0;
        ttt=0;
        'fail'
    end
    %PowerFaker(:,ttt)=  expm(dt*(ttt-1)*(AdjustW_jk + (alpha-recoveryRate)*eye(NumCities) ))*PowerFaker(:,1);
end


TimeBreakOuts(:,rrr) = TimeBreakOut;
end

TimeBreakOut=mean(TimeBreakOuts,2);

TimeBreakOutSingle=TimeBreakOuts(:,end);
simulationTime=toc();

    
c = [0.5,0.8,1.0; 1,0,0;0.0,0.0,0; ones(k,1)*[1,0.8,0]; 0.5,0.8,1.0];

ExpFakeBreakoutTimes(3)=0;
tic()

keepOnlyValid=(~isnan(ExpFakeBreakoutTimes))

%close all
figure(1)
subplot( 2 , 2 , 2 ) 
scatter(ExpFakeBreakoutTimes(1),TimeBreakOut(1),40,c(1,:),'lineWidth',1,'marker','^');
hold on
scatter(ExpFakeBreakoutTimes(2),TimeBreakOut(2),40,c(2,:),'lineWidth',1,'marker','d');
scatter(ExpFakeBreakoutTimes(3),TimeBreakOut(3),40,c(3,:),'lineWidth',1,'marker','o','MarkerFaceColor','k');
scatter(ExpFakeBreakoutTimes(4:(end-1)),TimeBreakOut(4:(end-1)),40,c(4:(end-1),:),'lineWidth',1,'marker','+');
scatter(ExpFakeBreakoutTimes(end),TimeBreakOut(end),40,c(end,:),'lineWidth',1,'marker','v');

xlabel('BP AT')
ylabel('Mean AT')
metric=ExpFakeBreakoutTimes(keepOnlyValid);
p=polyfit(metric',TimeBreakOut(keepOnlyValid),1)
yfit = polyval(p,[0,ceil(max(metric)*1.05)]);
hold on;
plot([0,ceil(max(metric)*1.05)],yfit,'k:')

txt = ['\gamma = 10^{', num2str(log10(gam))  '}'];
text(0.5,0.20,txt,'Units','normalized')

keepOnlyValid=true(k+4,1)

subplot( 2 , 2 , 1 ) 
scatter(D_brock(startNode,1),TimeBreakOut(1),40,c(1,:),'lineWidth',1,'marker','^');
hold on
scatter(D_brock(startNode,2),TimeBreakOut(2),40,c(2,:),'lineWidth',1,'marker','d');
scatter(D_brock(startNode,3),TimeBreakOut(3),40,c(3,:),'lineWidth',1,'marker','o','MarkerFaceColor','k');
scatter(D_brock(startNode,4:(end-1)),TimeBreakOut(4:(end-1)),40,c(4:(end-1),:),'lineWidth',1,'marker','+');
scatter(D_brock(startNode,end),TimeBreakOut(end),40,c(end,:),'lineWidth',1,'marker','v');

xlabel('ED AT')
ylabel('Mean AT')
metric=D_brock(startNode,(keepOnlyValid))';
p=polyfit(metric,TimeBreakOut(keepOnlyValid),1)
yfit = polyval(p,[0,ceil(max(metric)*1.05)]);
hold on;
plot([0,ceil(max(metric)*1.05)],yfit,'k:')

txt = ['\gamma = 10^{', num2str(log10(gam))  '}'];
text(0.5,0.20,txt,'Units','normalized')

figuretime=toc();





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



tic();

keepOnlyValid=true(k+4,1)

NumCities=4+k;
WidthWindow=1200;
Threshold=1;
repetition=150;

Nk= ones(NumCities,1)*80000000;

Network = speye(NumCities,NumCities);
Network(1,2)=1;
Network(2,3)=1;
Network(4:end-1,3)=1/k;
Network(end,4:end)=1;
%% This makes A network. Not nessesarily the most amazing one, but hey, its a network, it does its thing.

%Network= exp(rand(NumCities,NumCities)*4); %Lets claim that this is the number of people transfered each day...
%maxs= maxk(Network,5);
%Network= (Network.* (Network> maxs(end,:)));


[i,j,links] = find(Network)


%Network=(Network+Network');%%Make Symetic, yay.
Network=Network-diag(diag(Network));
%Network=Network.*Nk';
%Network=Network./sum(Network).*exp(rand(1,NumCities)*3);
setUpTime=toc();
'Set Up complete'

tic();
Ptrans=Network./sum(Network);
d_jk= 1 - log(Ptrans);
d_jk(d_jk==inf)=0;
d_jk=sparse(d_jk);
[D_brock,zob] = all_shortest_paths(d_jk,struct('algname','floyd_warshall'));
D_brock=D_brock';

d_jk= (Ptrans>0)*1.00;
d_jk=sparse(d_jk);
[D_dumb,zob] = all_shortest_paths(d_jk,struct('algname','floyd_warshall'));
D_dumb=D_dumb'

basicDistanceTime=toc();
 'have made short paths'

alpha=0.5;
gam=0.00001;

AdjustW_jk=sparse(gam*(Network-diag(sum(Network))));

startNode=3;

HistoryInfectionLevel= zeros(NumCities,WidthWindow);
HistoryInfectionLevel(startNode,1)=1;

PowerFaker=HistoryInfectionLevel;

RecoveredLevel=zeros(NumCities,1);

recoveryRate=0;


tic()
T=linspace(0,40);

AnalyticSurvivalCurve=zeros(k+4,length(T));

FunctionFilter=[1:k+4]
for(iii=FunctionFilter)    
y0=ones(size(Network,1),1);
y0(iii)=0;
deriv= @(t,y) SurvivalDeriv(t,y,gam*Network',alpha,iii,recoveryRate)
[~,ode45y]=ode45(deriv, T, y0,odeset('NonNegative',1));
AnalyticSurvivalCurve(iii,:)=ode45y(:,startNode)';
end
AnalyticSurvivalCurve=AnalyticSurvivalCurve.^2;

AnalyticSurvivalCurveConditioned=(AnalyticSurvivalCurve-AnalyticSurvivalCurve(:,end))./(AnalyticSurvivalCurve(:,1)-AnalyticSurvivalCurve(:,end))
ExpFakeBreakoutTimes=sum(AnalyticSurvivalCurveConditioned')*T(2)


tic();

TimeBreakOuts=zeros(NumCities,repetition);

for (rrr=1:repetition)
    rrr

RecoveredLevel=zeros(NumCities,1);
InfectedLevel=zeros(NumCities,1);
InfectedLevel(startNode)=2;

TimeBreakOut=-ones(NumCities,1);
PopBreakOut=-ones(NumCities,1);
time=0;


tic();
ttt=0;
while(true)
    
    ttt=ttt+1;
    time=time+dt;
    
    nowRecov= min(poissrnd(recoveryRate*dt*InfectedLevel),InfectedLevel);
    nowSick= poissrnd(alpha.*dt.*InfectedLevel.* (Nk-InfectedLevel-RecoveredLevel)./Nk);
    nowMoveForwardSeed= poissrnd(gam*links.*dt.*InfectedLevel(i) );
    nowMoveReverse= poissrnd(gam*links.*dt.*InfectedLevel(j) );
    motionMatrix= sparse(i,j,nowMoveForwardSeed)+sparse(j,i,nowMoveReverse);
    moveGain=sum(motionMatrix,1)';
    moveLoss=sum(motionMatrix,2);
    InfectedLevel= max(InfectedLevel+ moveGain-moveLoss+nowSick-nowRecov,0);
    RecoveredLevel=RecoveredLevel+nowRecov;
    
    NowOutbreak= (TimeBreakOut==-1) & (InfectedLevel>=Threshold);
    TimeBreakOut(NowOutbreak)=time-dt;
    
    
    if(all( TimeBreakOut>=0) )
        break;
    end
    
    
    if(all(InfectedLevel==0) ) %%Simulation failure. Try again.
        InfectedLevel=zeros(NumCities,1);
        InfectedLevel(startNode)=1;
        TimeBreakOut=-ones(NumCities,1);
        time=0;
        ttt=0;
        'fail'
    end
    
    if(any( InfectedLevel<0) ) %%Simulation failure. Try again.
        InfectedLevel=zeros(NumCities,1);
        InfectedLevel(startNode)=1;
        TimeBreakOut=-ones(NumCities,1);
        time=0;
        ttt=0;
        'fail'
    end
    %PowerFaker(:,ttt)=  expm(dt*(ttt-1)*(AdjustW_jk + (alpha-recoveryRate)*eye(NumCities) ))*PowerFaker(:,1);
end


TimeBreakOuts(:,rrr) = TimeBreakOut;
end

TimeBreakOut=mean(TimeBreakOuts,2);

TimeBreakOutSingle=TimeBreakOuts(:,end);
simulationTime=toc();

    
c = [0.5,0.8,1.0; 1,0,0;0.0,0.0,0; ones(k,1)*[1,0.8,0]; 0.5,0.8,1.0];

ExpFakeBreakoutTimes(3)=0;
tic()

%close all
figure(1)
subplot( 2 , 2 , 4 ) 
scatter(ExpFakeBreakoutTimes(1),TimeBreakOut(1),40,c(1,:),'lineWidth',1,'marker','^');
hold on
scatter(ExpFakeBreakoutTimes(2),TimeBreakOut(2),40,c(2,:),'lineWidth',1,'marker','d');
scatter(ExpFakeBreakoutTimes(3),TimeBreakOut(3),40,c(3,:),'lineWidth',1,'marker','o','MarkerFaceColor','k');
scatter(ExpFakeBreakoutTimes(4:(end-1)),TimeBreakOut(4:(end-1)),40,c(4:(end-1),:),'lineWidth',1,'marker','+');
scatter(ExpFakeBreakoutTimes(end),TimeBreakOut(end),40,c(end,:),'lineWidth',1,'marker','v');


xlabel('BP AT')
ylabel('Mean AT')
metric=ExpFakeBreakoutTimes(keepOnlyValid);
p=polyfit(metric',TimeBreakOut(keepOnlyValid),1)
yfit = polyval(p,[0,ceil(max(metric)*1.05)]);
hold on;
plot([0,ceil(max(metric)*1.05)],yfit,'k:')

txt = ['\gamma = 10^{', num2str(log10(gam))  '}'];
text(0.5,0.20,txt,'Units','normalized')

keepOnlyValid=true(k+4,1)

subplot( 2 , 2 , 3 ) 
scatter(D_brock(startNode,1),TimeBreakOut(1),40,c(1,:),'lineWidth',1,'marker','^');
hold on
scatter(D_brock(startNode,2),TimeBreakOut(2),40,c(2,:),'lineWidth',1,'marker','d');
scatter(D_brock(startNode,3),TimeBreakOut(3),40,c(3,:),'lineWidth',1,'marker','o','MarkerFaceColor','k');
scatter(D_brock(startNode,4:(end-1)),TimeBreakOut(4:(end-1)),40,c(4:(end-1),:),'lineWidth',1,'marker','+');
scatter(D_brock(startNode,end),TimeBreakOut(end),40,c(end,:),'lineWidth',1,'marker','v');

xlabel('ED AT')
ylabel('Mean AT')
metric=D_brock(startNode,(keepOnlyValid))';
p=polyfit(metric,TimeBreakOut(keepOnlyValid),1)
yfit = polyval(p,[0,ceil(max(metric)*1.05)]);
hold on;
plot([0,ceil(max(metric)*1.05)],yfit,'k:')

txt = ['\gamma = 10^{', num2str(log10(gam))  '}'];
text(0.5,0.20,txt,'Units','normalized')

figuretime=toc();
