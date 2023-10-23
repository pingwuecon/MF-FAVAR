function q = iGentropy(y,mu,lambda)
fun = @(y,mu,lambda) -((((lambda./(2*pi.*y.^3)).^0.5).*exp((-lambda.*(y-mu).^2)./(2*(mu.^2).*y))).*(0.5*log(lambda./(2*pi*y.^3)) ...
- (lambda.*(y - mu).^2)./(2*(mu.^2).*y))) ;
% q = integral(@(y)fun(y,mu,lambda),0,Inf);
q = integral(@(y)fun(y,mu,lambda),0,Inf,'ArrayValued',true);
end



