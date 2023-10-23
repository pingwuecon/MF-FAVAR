function q = ElogV(y,mu,lambda)
%% m = 1 Lasso, Adaptive Lasso
% fun = @(y,mu,lambda) pdf('InverseGaussian',y,mu,lambda).*log(y);
fun = @(y,mu,lambda) (((lambda./(2*pi.*y.^3)).^0.5).*exp((-lambda.*(y-mu).^2)./(2*(mu.^2).*y))).*log(y);
q = integral(@(y)fun(y,mu,lambda),0,Inf,'ArrayValued',true);  
end