function lden = lgammpdf(y,a,b)
lden = a.*log(b) - gammaln(a) - (a-1) .* log(y) - b.*y;
end
