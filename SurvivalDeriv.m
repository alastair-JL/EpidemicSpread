function [deriv] = SurvivalDeriv(t,y,F,alpha,pit,recoveryRate)

     if(~exist('recoveryRate','var'))
        recoveryRate=0; 
     end

    deriv=   F*y+ alpha.*y.*y +recoveryRate -(sum(F,2).*y+alpha.*y+recoveryRate.*y);
    deriv(pit)=0;
    
end