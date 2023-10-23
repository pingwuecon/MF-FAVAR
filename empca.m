function [F,Z,iter] = empca(X,nf1,nf)
%PCA coupled with EM, by Stock and Watson
[T,N] = size(X);
Z = X;
A = cell(N,1);
X_ob = cell(N,1);
I = eye(T);

for i=1:N
    Li = find(isnan(X(:,i)))';
    Ai = I;
    Ai(Li,:)=[];
    A{i} = Ai;
    Xi = X(:,i);
    Xi(Li,:)=[];
    X_ob{i} = Xi;
%     Z(Li,i) = randn(length(Li),1);
    Z(Li,i) = Xi(1);
end

% Initial values
[F,L] = extract2(Z,nf1);Zold=Z;
flag=0;i=1;iter=1;
while i<2000 && flag==0
    % Expectation    
    PredG=F*L';  
    for j=1:N
        Xtmp = X_ob{j};
        ntmp = length(Xtmp);
        if ntmp<T
            Atmp = A{j};
            Z(:,j) = PredG(:,j) + (Atmp'/(Atmp*Atmp'))*(Xtmp-Atmp*PredG(:,j));
        end
    end
    perchange = (Z-Zold)./Zold;
    conver=max(perchange);
    if abs(conver)<.01
        flag=1;iter=i;
    else
        Zold=Z;i=i+1;flag=0;
    end
    % Maximization
    [F,L] = extract2(Z,nf1);
end
[F,L] = extract2(Z,nf);