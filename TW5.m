function [F,X] = TW5(X,r1,r2,r3)
[T,N] = size(X);
w = ~isnan(X);
Tallidtmp = sum(w,1);
Tallid = find(Tallidtmp==T);
XTall = X(:,Tallid);
[TTall,NTall] = size(XTall);
Wideidtmp = sum(w,2);
Wideid = find(Wideidtmp==N);
XWide = X(Wideid,:);
[TWide,NWide] = size(XWide);
if r1>NTall
r1=NTall;
end
[FTall,LTall] = extract2(XTall,r1);
[FWide,LWide] = extract2(XWide,r2);

LWidetmp = LWide(Tallid,:);
% HMiss = (LWidetmp'*LWidetmp)\(LWidetmp'*LTall);
HMiss = (LTall'*LTall)\(LTall'*LWidetmp);
CMiss = FTall*HMiss*LWide';
X(~w) = CMiss(~w);

[F,L] = extract2(X,r3);
end