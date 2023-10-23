function q = ElogV1(y,a1,b1,x,a0,b0)
%% Horseshoe
%  fun = @(x,y) invgammpdf(y,a1,b1).*invgammpdf(x,a0,b0).*log(x.*y);    
% q = integral2(fun,0,Inf,0,Inf);    
% 
% fun1 = @(y) invgammpdf(y,a1,b1).*log(y);
% q1 = integral(fun1,0,10^5); 
% 
% fun2 = @(x) invgammpdf(x,a0,b0).*log(x);
% q2 = integral(fun2,0,10^5); 

% fun3 = @(y) invgammpdf(y,a1,b1);
% q3 = integral(fun3,0,Inf); 
% 
% fun4 = @(x) invgammpdf(x,a0,b0);
% q4 = integral(fun4,0,Inf); 

% q = inversegamcdfgam(y,a1,b1)*q2 + q1*inversegamcdfgam(x,a0,b0);


fun1 = @(y) exp(linvgammpdf(y,a1,b1)).*log(y);
q1 = integral(fun1,0,Inf,'ArrayValued',true); 

fun2 = @(x) exp(linvgammpdf(x,a0,b0)).*log(x);
q2 = integral(fun2,0,Inf,'ArrayValued',true); 

fun3 = @(y) exp(linvgammpdf(y,a1,b1));
q3 = integral(fun3,0,Inf,'ArrayValued',true); 

fun4 = @(x) exp(linvgammpdf(x,a0,b0));
q4 = integral(fun4,0,Inf,'ArrayValued',true); 

q = q3.*q2 + q1.*q4;
end