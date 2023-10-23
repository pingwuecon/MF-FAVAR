function [F,X] = TP3(X,r1,r3)
[T,N] = size(X);
w = ~isnan(X);
Tallidtmp = sum(w,1);
Tallid = find(Tallidtmp==T);
XTall = X(:,Tallid);
[TTall,NTall] = size(XTall);
if r1>NTall
r1=NTall;
end
[FTall,LTall] = extract2(XTall,r1);
for i=1:N
    Xi = X(:,i);
    Li = find(~isnan(Xi))';
    Mi = find(isnan(Xi))';
    X_ob = Xi(Li,1);
    Ftmp = FTall(Li,:);
    Lambda =  (Ftmp'*Ftmp)\(Ftmp'*X_ob);
    CMiss = FTall*Lambda;
    X(Mi,i) = CMiss(Mi,1);
end

[F,L] = extract2(X,r3);

end