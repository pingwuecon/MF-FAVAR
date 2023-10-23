function lden = gamment(a,b)
lden = a - log(b) + gammaln(a) + (1-a).*psi(a);
end
