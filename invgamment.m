function lden = invgamment(a,b)
lden = a + log(b) + gammaln(a) - (a+1).*psi(a);
end
