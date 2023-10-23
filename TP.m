function [F,Z] = TP(X,r1,r3)
[T,N] = size(X);
w = ~isnan(X);
Tallidtmp = sum(w,1);
Tallid = find(Tallidtmp==T);
XTall = X(:,Tallid);
[TTall,NTall] = size(XTall);

if NTall==0
X1 = X(1:end-2,:);
X1(:, all(isnan(X1),1)) = [];
[F,Z1] = TP3(X1,r1,r3);
[F,Z] = TP3(Z1,r1,r3);
else
[F,Z] = TP3(X,r1,r3);
end

end