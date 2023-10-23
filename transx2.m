function [y]=transx2(X,Tcode)
%    Transform x
%    Return Series with same dimension and corresponding dates
%    Missing values where not calculated
%    -- Tcodes:
%             1 Level
%             5 Log-First-Difference

%  Translated from the Gauss procs of Stock&Watson(2005),'Implications of
%  dynamic factor models for VAR analysis'
%  Dimitris Korobilis, June 2007

[n,N]=size(X);
y=zeros(n,N);        %storage space for y

for i=1:N
    x = X(:,i);
    tcode = Tcode(i);
if tcode == 1
    y(:,i)=x;
elseif tcode == 5
    y(1,i)=nan;
    x(x==0,:)=1.0e-04;
    for t=2:n
        x1 = x(t);
        x0 = x(t-1);
        if isnan(x1) 
            y(t,i)=nan;
        elseif isnan(x0)
            y(t,i)=nan;
        else
            y(t,i)=log(x1)-log(x0);
        end
    end
end
end