function [vecY,k,k2,M_o,M_u,Y_confinal,M_a, A,idn]=create_matrixVBtestwithRSTIQ1test(ind,Y_2,lag,m,RSTI,ITL1shares,Y_con)
RSTI = RSTI(182:end,:);
[n,T]=size(Y_2);
R = 12;

% Y_con(46:end,1:12) = NaN;

vecY=Y_2(ind==1);
k=length(vecY);
k2=T*n-k;
M_o=sparse(find(ind==1),1:k,ones(k,1));
M_o=[M_o;sparse(T*n-size(M_o,1),k)];
M_u=sparse(find(ind==0),1:(k2),ones(k2,1));
M_u=[M_u;sparse(T*n-size(M_u,1),k2)];
A=cell(lag+1);
A{1}=[sparse(T-lag,lag),speye(T-lag)];
for j=1:lag
A{j+1}=[sparse(T-lag,lag-j),speye(T-lag),sparse(T-lag,j)];
end
Tannual = size(Y_con,1);
k3=(Tannual-1)*(n-m);
% Y_con = Y_2(m+1:end,:)';
% Y_con( ~any(Y_con,2), : ) = [];


%% ITL1 cross-sectional restrictions
ITL1 = [ITL1shares zeros(T,(n-m) - R)]; % Cross-sectional ILT2 
Ma_ITL1cs = SURform3(ITL1,1);
Y_ITL1cs = Y_2(1,:)';

%% RSTI
A1 = [ones(T,1) zeros(T,11) zeros(T,(n-m) - R)]; % NE
A2 = [zeros(T,1) ones(T,1) zeros(T,10) zeros(T,(n-m) - R)]; % NW
A3 = [zeros(T,2) ones(T,1) zeros(T,9) zeros(T,(n-m) - R)]; % YH
A4 = [zeros(T,3) ones(T,1) zeros(T,8) zeros(T,(n-m) - R)]; % EM
A5 = [zeros(T,4) ones(T,1) zeros(T,7) zeros(T,(n-m) - R)]; % WM
A6 = [zeros(T,5) ones(T,1) zeros(T,6) zeros(T,(n-m) - R)]; % NE
A7 = [zeros(T,6) ones(T,1) zeros(T,5) zeros(T,(n-m) - R)]; % NE
A8 = [zeros(T,7) ones(T,1) zeros(T,4) zeros(T,(n-m) - R)]; % NE
A9 = [zeros(T,8) ones(T,1) zeros(T,3) zeros(T,(n-m) - R)]; % NE
A10 = [zeros(T,9) ones(T,1) zeros(T,2) zeros(T,(n-m) - R)]; % NE
A11 = [zeros(T,10) ones(T,1) zeros(T,1) zeros(T,(n-m) - R)]; % NE
A12 = [zeros(T,11) ones(T,1) zeros(T,(n-m) - R)]; % NE

A1new = SURform3(A1,1);
A2new = SURform3(A2,1);
A3new = SURform3(A3,1);
A4new = SURform3(A4,1);
A5new = SURform3(A5,1);
A6new = SURform3(A6,1);
A7new = SURform3(A7,1);
A8new = SURform3(A8,1);
A9new = SURform3(A9,1);
A10new = SURform3(A10,1);
A11new = SURform3(A11,1);
A12new = SURform3(A12,1);

Ma_RSTI = [A1new(182:end-1,:);A2new(182:end-1,:);A3new(182:end-1,:);A4new(182:end-1,:);A5new(182:end-1,:);A6new(182:end-1,:);A7new(182:end-1,:);A8new(182:end-1,:);...
    A9new(182:end-1,:);A10new(182:end-1,:);A11new(182:end,:);A12new(182:end,:)];
Y_RSTI =[reshape(RSTI(1:end-2,1:10),size(RSTI(1:end-2,1:10),1)*10,1);reshape(RSTI(1:end-1,11:12),size(RSTI(1:end-1,11:12),1)*2,1)];

% Y_RSTI = reshape(RSTI,size(RSTI,1)*R,1);


%% Temporal restrictions
Temp_con=[kron(speye(Tannual-1),kron([0,1/4, 1/2,3/4],speye(n-m))),sparse(k3,4*(n-m))]...
    +[sparse(k3,4*(n-m)),kron(speye(Tannual-1),kron([1,3/4,1/2,1/4],speye(n-m)))]; % Temporal restrictions

Y_temphat=reshape(Y_con(2:end,:)',k3,1);
idn1 = find(isnan(Y_temphat)==0);
Y_temp = Y_temphat(idn1,:);
Ma_temp = Temp_con(idn1,:);
% Y_temp = [reshape(Y_con(2:32,1:12)',31*12,1);reshape(Y_con(33:Tannual,:)',size(Y_con(33:Tannual,:),1)*(n-m),1)]; 
Ma_temp = [Ma_temp zeros(size(Ma_temp,1),T*(n-m) - length(Ma_temp))];
Mafinal = [Ma_ITL1cs;Ma_RSTI;Ma_temp];
Y_confinal = [Y_ITL1cs;Y_RSTI;Y_temp];
% [nr nl] = size(Temp_con);
% 
% Temp_con = [Temp_con zeros(nr,(T*R) - nl) ];
% Y_temp = reshape(Y_con(2:Tannual-1+1,:)',(Tannual-1)*12,1); 
% 
% Mafinal = [Ma_ITL1cs;Ma_RSTI;Temp_con];
% Y_confinal = [Y_ITL1cs;Y_RSTI;Y_temp];

idn = [ones(length(Y_ITL1cs),1); reshape(repmat(2:12,T,1),T*(R-1),1); zeros(length(Y_RSTI),1);zeros(length(Y_temp),1)];

M_a=sparse(Mafinal)';