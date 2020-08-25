
NumPoints=30;
finalT=10;
x=[0,sort(finalT*rand(1,NumPoints)),finalT];
v=([rand(1,NumPoints+1),0]).^5;

extinctionProbability=rand()*0.2;
xq=linspace(0,finalT,8000);
vq = interp1(x,v,xq,'nearest');
%vq = interp1(x,v,xq);

vq=-vq/(sum(vq)*xq(2))*(1-extinctionProbability)

yq=1+cumsum(vq)*xq(2);



MatchingFineness=500;
MatchingPoints= linspace(1,extinctionProbability,MatchingFineness+1);
MatchingPoints=MatchingPoints+(1-extinctionProbability)/(MatchingFineness+2);
MatchingPoints=MatchingPoints(2:end);

IndexPerBucket= sum(MatchingPoints'-yq<0,2);
IndexPerBucket(1)=1;

slopes=vq(IndexPerBucket);


means=linspace(0,finalT*1,MatchingFineness)
alphas= 0*means+MatchingFineness./finalT;

T=linspace(0,finalT*1.1,1000)'

gams= exp(-means.*alphas);
for(qqq=1:100)
gams= gams- (gams-(2*gams+alphas).*exp(-means.*(alphas+gams)))./(1+(2*gams+alphas).*means.*exp(-means.*(alphas+gams)));
end
Q= (gams+alphas)./(alphas + (2*gams+alphas).*exp((alphas+gams).*(T-means) )  );

WeightPerBucket= interp1(x,v,means,'nearest');

WeightPerBucket=WeightPerBucket./sum(WeightPerBucket).*(1-extinctionProbability)

figure()
plot(T,sum(WeightPerBucket.*Q,2)+extinctionProbability,'g','lineWidth',2)
hold on
plot(xq,yq,'k')
legend({'Approximation','f(t)'})