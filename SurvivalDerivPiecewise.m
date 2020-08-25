function [deriv] = SurvivalDerivPiecewise(t,y,F,alpha,pit,recoveryRate,alphaMult,TimeThresholds)

     if(~exist('recoveryRate','var'))
        recoveryRate=0; 
     end
     
     index=sum(TimeThresholds<t);
     index=max(index,1);
     index=min(index,length(alphaMult));
    alpha=alpha*alphaMult(index);
     
    deriv=   F*y+ alpha.*y.*y +recoveryRate -(sum(F,2).*y+alpha.*y+recoveryRate.*y);
    deriv(pit)=0;
    
end