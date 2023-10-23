function [taunew,entropybeta,logVbeta] = lasso1(km,lambdabar,betanew,kbetanew)

taunew = zeros(km,1);
% entropybeta = zeros(km,1);
logVbeta = zeros(km,1);

taunew = (sqrt(lambdabar./((betanew.^2) + kbetanew)));
% logVbeta = log(1./taunew);

entropybeta = iGentropy(1./taunew,1./(sqrt(lambdabar./((betanew.^2) + kbetanew))) + 1./lambdabar,lambdabar);
 logVbeta = ElogV(1./taunew,1./(sqrt(lambdabar./((betanew.^2) + kbetanew))),lambdabar);
%  parfor tv = 1:km 
% %  entropybeta(tv,1) = iGentropy(1./taunew(tv,1),1./(sqrt(lambdabar(tv)/((betanew(tv)^2) + kbetanew(tv)))) + 1/lambdabar(tv),lambdabar(tv));
%   logVbeta(tv,1) = ElogV(1./taunew(tv,:),1./(sqrt(lambdabar(tv)/((betanew(tv)^2) + kbetanew(tv)))),lambdabar(tv));
%  end



end