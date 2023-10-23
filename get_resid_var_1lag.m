function sig2 = get_resid_var_1lag(Y0,Y)
[T,n] = size(Y);
sig2 = zeros(n,1);
tmpY = [Y0; Y];
t = size(tmpY,1);
for i=1:n
    Z = [ones(t-1,1) tmpY(1:end-1,i)];
    tmpb = (Z'*Z)\(Z'*tmpY(2:end,i));
    sig2(i) = mean((tmpY(2:end,i)-Z*tmpb).^2);
end