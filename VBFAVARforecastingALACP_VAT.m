%% Variational bayes for the Adaptive Lasso with palleralisation (with missing data in the beginning) with cross-sectional on ITL1 regions
clear all; clc;
warning off
addpath('./Data')

load('Annual2.mat')
load('Annualdat.mat')
load('ITL1weights.mat')
load('Quarterly.mat')
load('RSTIdata.mat')
load('fcast.mat')
load('ITL1fcast2.mat')
load('Quarterlydat.mat')
load('VAT.mat')
load('ITL1id.mat')

icount = 1;
Q1count = 1;
Q2count = 1;
Q3count = 1;

Q1countITL1 = 1;
Q2countITL1 = 1;
Q3countITL1 = 1;

store_nowcast_error = []; 
store_backcast_error = [];
store_forecast_error = [];

for R = 1:32

YQ0 = Quarterly{R,1};
RSTInew = RSTIdata{R+1,1};

RSTIann = [];

for tv = 7:size(RSTInew,1)
   RSTIann = [ RSTIann ;.25*RSTInew(tv,:) + .5*RSTInew(tv-1,:) + .75*RSTInew(tv-2,:) + RSTInew(tv-3,:) + .75*RSTInew(tv-4,:) + .5*RSTInew(tv-5,:) + .25*RSTInew(tv-6,:)];
end
RSTIann = [zeros(6,12) ; RSTIann(:,2:end)];


for sR = 1:3

    store_nowcastITL1 = [];
    store_backcastITL1 = [];
    store_forecastITL1 = [];


if Annualdat(icount) == 2012
      ITL1shares =  ITL1weights{1,1};
      YA0 = Annual2{1,1};
      Tint = find(YQ0(:,1)==2013);

elseif Annualdat(icount) == 2013
      ITL1shares =  ITL1weights{2,1};
      YA0 = Annual2{2,1};
      Tint = find(YQ0(:,1)==2013);

elseif Annualdat(icount) == 2014
       ITL1shares =  ITL1weights{3,1};
      YA0 = Annual2{3,1};
      Tint = find(YQ0(:,1)==2013);

elseif Annualdat(icount) == 2015
       ITL1shares =  ITL1weights{4,1};
      YA0 = Annual2{4,1};
      Tint = find(YQ0(:,1)==2013);

 elseif Annualdat(icount) == 2016
       ITL1shares =  ITL1weights{5,1};
      YA0 = Annual2{5,1};
      Tint = find(YQ0(:,1)==2013);

   elseif Annualdat(icount) == 2017
       ITL1shares =  ITL1weights{6,1};
      YA0 = Annual2{6,1};
      Tint = find(YQ0(:,1)==2013);
      
    elseif Annualdat(icount) == 2018
       ITL1shares =  ITL1weights{7,1};
      YA0 = Annual2{7,1};
       Tint = find(YQ0(:,1)==2013);
    else
      ITL1shares =  ITL1weights{8,1};
      YA0 = Annual2{8,1};
      Tint = find(YQ0(:,1)==2013);
    end


if sR == 1
YQ = YQ0(1:end-1,2:end);
enddate = YQ0(end-1,1);
else
YQ = YQ0(1:end,2:end);
enddate = YQ0(end,1);
end
T = length(YQ);
ITL1shares = kron(ITL1shares,ones(4,1));

na = 12; % no. of annual variables

p = 1; % no .of lags (endogenous variables)
YA = kron(YA0(:,2:end),[0 0 0 1]');
Tnew = length(YQ) - length(YA);
if Tnew > 0
   YA = [YA;ones(Tnew,1)*YA(end,:)];
   ITL1shares = [ITL1shares;ones(Tnew,1)*ITL1shares(end,:)];
end
%% eigen
FACTOR = 1;
nf=3;
nq = 5+12*nf; % no. of quarterly variables

regionalpredic;
nf1 = 2;
nf2 = 3;
vatlag = 5;
Xpredictors = cellfun( @(a) a(1:189+R,:), regionalpred,'UniformOutput',false);
if R>31
    VATpredictor = VAT{31};
else
    VATpredictor = VAT{R};
end
[Ficount corrx FX] = Factor_est_vat(Xpredictors,idtranlag,nf1,nf2,3,sR,ITL1id,VATpredictor,vatlag);
Fmatrix=[];
F_sel=1;
if F_sel==1 % eig
    for iF=1:12
        Ftmp=Ficount{iF};
        Fmatrix=[Fmatrix Ftmp(:,1)];
    end
    if nf>1
        for iF=1:12
            Ftmp=Ficount{iF};
            Fmatrix=[Fmatrix Ftmp(:,2)];
        end
    end
    if nf>2
        for iF=1:12
            Ftmp=Ficount{iF};
            Fmatrix=[Fmatrix Ftmp(:,3)];
        end
    end
    if nf>3
        for iF=1:12
            Ftmp=Ficount{iF};
            Fmatrix=[Fmatrix Ftmp(:,4)];
        end
    end
    if nf>4
        for iF=1:12
            Ftmp=Ficount{iF};
            Fmatrix=[Fmatrix Ftmp(:,5)];
        end
    end
elseif F_sel==2 % corr
    for iF=1:12
        Ftmp=Ficount{iF};
        corr = zeros(nf2,1);
        for idf=1:nf2
            corrtmp = corrcoef(RSTInew(182:T,iF),Ftmp(182:T,idf));
            corr(idf) = corrtmp(2,1);
        end
        [~,corrorder] = sort(abs(corr),'descend');
        Fmatrix=[Fmatrix Ftmp(:,corrorder(1))];
    end
    if nf>1
        for iF=1:12
            Ftmp=Ficount{iF};
            corr = zeros(nf2,1);
            for idf=1:nf2
                corrtmp = corrcoef(RSTInew(182:T,iF),Ftmp(182:T,idf));
                corr(idf) = corrtmp(2,1);
            end
            [~,corrorder] = sort(abs(corr),'descend');
            Fmatrix=[Fmatrix Ftmp(:,corrorder(2))];
        end
    end
    if nf>2
        for iF=1:12
            Ftmp=Ficount{iF};
            corr = zeros(nf2,1);
            for idf=1:nf2
                corrtmp = corrcoef(RSTInew(182:T,iF),Ftmp(182:T,idf));
                corr(idf) = corrtmp(2,1);
            end
            [~,corrorder] = sort(abs(corr),'descend');
            Fmatrix=[Fmatrix Ftmp(:,corrorder(3))];
        end
    end
    if nf>3
        for iF=1:12
            Ftmp=Ficount{iF};
            corr = zeros(nf2,1);
            for idf=1:nf2
                corrtmp = corrcoef(RSTInew(182:T,iF),Ftmp(182:T,idf));
                corr(idf) = corrtmp(2,1);
            end
            [~,corrorder] = sort(abs(corr),'descend');
            Fmatrix=[Fmatrix Ftmp(:,corrorder(4))];
        end
    end
end
FMatrix=Fmatrix;
Fmatrix=FMatrix(1:T,:);
Ydata = [YQ Fmatrix YA];

ILT1w = ITL1shares;
T = length(Ydata);
ind = [ones(T,nq) zeros(T,na)]';
if sR == 1
[vecY,kdim,k2,M_o,M_u,Y_con,M_a,A,idn]=create_matrixVBtestwithRSTIQ1test(ind,Ydata',p,nq,RSTInew(1:end-1,2:end),ITL1shares,YA0(:,2:end));
elseif sR ==2
[vecY,kdim,k2,M_o,M_u,Y_con,M_a,A,idn]=create_matrixVBtestwithRSTIQ2test(ind,Ydata',p,nq,RSTInew(1:end-1,2:end),ITL1shares,YA0(:,2:end));
else
[vecY,kdim,k2,M_o,M_u,Y_con,M_a,A,idn]=create_matrixVBtestwithRSTIQ3test(ind,Ydata',p,nq,RSTInew(1:end-1,2:end),ITL1shares,YA0(:,2:end));
end

y0 = Ydata(1:p,:);
y = Ydata(p+1:end,:);
[T,n] = size(y);
longy = reshape(y',T*n,1);
Ridn=cell(n,1);
createRidn;
% 
k = n+p*n^2;
np = n*p;
m = n*(n - 1)/2;
Lid = find(tril(ones(n,n),-1))';
L = eye(n);
warning off
%% Construct X
X = zeros(T,np); 
tempY = [y0; y];
for i=1:p
    X(:,(i-1)*n+1:i*n) = tempY(p-i+1:end-i,:);
end
X = [ones(T,1) X];

Xqlag = zeros(T,nq*p); 
for i=1:p
    Xqlag(:,(i-1)*nq+1:i*nq) = tempY(p-i+1:end-i,1:nq);
end

%% Priors
cserr  = .00001;

Vbeta = 1;alphabeta = 0*ones(k,1);
Va = 1;
nu0 = 5 ; S0 = .01*(nu0-1)*ones(n,1);
lambdaa = 1; taua = 10; 
a0 = .1 ; b0 = .001;
S0cs = .01; nucs = 1000;

s2 = get_resid_var_1lag(y0,y);
%% Variational bayes
%  Sbar = S0;
Sbar_store = ones(n,1)*10;
invsig2_store = ones(n,1)*10;
nubar = T+n;
check  = 10;
pyhat =[];

pyhat0 = zeros(T,na);
pyhat1= 0;
count = 0;
Kbetabar_store = cell(n,1);
betabar_store = cell(n,1);
invVbeta_store = cell(n,1);
taubar_store = cell(n,1);
lambdabar_store = cell(n,1);
llike_store = zeros(n,1);
store_diff = zeros(n,1);

for iv = 1:n
    indic =  Ridn{iv,1} + 1;
    km = length(indic);
lambdabar_store{iv,1} = lambdaa;
taubar_store{iv,1} = taua*ones(km,1);
if iv ==1
invVbeta_store{1,1} = diag([1./Vbeta*ones(km,1)]);
else
invVbeta_store{iv,1} = diag([1./Vbeta*ones(length(indic),1);0.04./s2(1:iv-1,1)]);
end
end

tic()

count1 = 0;
chck1 = -1;
while chck1<0
parfor ii = 1:n
warning off    
check  = 10;
pyhat =[];     
if ii ==1

while check>1    
    
indic =  Ridn{ii,1} + 1;
zt =  X(:,indic);  
km = size(zt,2);     
% zt = [X];   
% km = k/n;
Sbar = Sbar_store(ii,1);

invVbeta = invVbeta_store{ii,1};
lambdabar = lambdabar_store{ii,1};
taubar = taubar_store{ii,1}; 

Kbetabar = (invVbeta + (nu0 + T/2)/Sbar*zt'*zt)\speye(km);
betabar =(nu0 + T/2)/Sbar*(Kbetabar*zt'*y(:,ii));
Sbar = S0(ii) + 0.5*norm(y(:,ii) - zt*betabar)^2 + 0.5*trace(zt'*zt*Kbetabar);

% Lasso priors
betanew = betabar;
lambdabar = (a0+1)./(b0 + (taubar)/2);
kbetanew = diag(Kbetabar );
[taunew,entropybeta,logVbeta] = lasso1(km,lambdabar,betanew,kbetanew);
taubar = 1./taunew(:,:);
invVbeta = diag(taunew);

Kbetabar = (Kbetabar + Kbetabar')/2;
gg = eig(full(Kbetabar))>0;
if sum(gg) == km 
CholKbeta = real(logdet(Kbetabar));
else
   CholKbeta = 0;
end

lliketau = sum(log(lambdabar/2) - lambdabar/2.*taubar) - sum(entropybeta);
llikelam = sum(lgammpdf(lambdabar,a0,b0) - gamment((a0+1),(b0 + (taubar)/2)));
    
llike1 = -T/2*log(2*pi) - T/2*log(Sbar/(nu0 + T/2)) - ((y(:,ii) - zt*betabar)'*(y(:,ii) - zt*betabar) + trace(zt*Kbetabar*zt'))/(2*Sbar/(nu0 + T/2)); % llikehood
llike2 = -km/2*log(2*pi)  - .5*sum(logVbeta) - 0.5*(betabar'*invVbeta*betabar + trace(invVbeta*Kbetabar)) ; % prior on theta

llike3 = linvgammpdf((Sbar/(nu0 + T/2)),nu0,S0(ii)); % prior on the sigma
llike4 = km/2 + km/2*log(2*pi) + 0.5*CholKbeta; % entropy on theta
llike5 = invgamment(nu0 + T/2,Sbar); % entropy on sigma
llike = llike1 + llike2 + llike3 - ( llike4 + llike5) + lliketau +  llikelam ;  
llike = llike/(T*n);

Kbetabar_store{ii,1} = Kbetabar;
betabar_store{ii,1} = betabar;
Sbar_store(ii,1) = Sbar;    
llike_store(ii,1) = llike;     
invVbeta_store{ii,1} = invVbeta ;
taubar_store{ii,1} = taubar ;
lambdabar_store{ii,1} =  lambdabar;
    
diff = abs(llike - pyhat);

if diff < 10^-4
    check = -10;
    store_diff(ii,:) = diff;
else
    pyhat = llike;
end
    
end
    

else
    

 while check>1       
    
indic =  Ridn{ii,1} + 1;  
km1 = length(indic);
zt = [X(:,indic) -y(:,1:ii-1)];  
km = size(zt,2);    
% km = k/n + ii-1;

% zt = [X -y(:,1:ii-1)];      
Sbar = Sbar_store(ii,1);
invVbeta = invVbeta_store{ii,1};
lambdabar = lambdabar_store{ii,1};
taubar = taubar_store{ii,1}; 

Kbetabar = (invVbeta + (nu0 + T/2)/Sbar*zt'*zt)\speye(km);
betabar =(nu0 + T/2)/Sbar*(Kbetabar*zt'*y(:,ii));
Sbar = S0(ii) + 0.5*norm(y(:,ii) - zt*betabar)^2 + 0.5*trace(zt'*zt*Kbetabar);
err=y(:,ii) - zt*betabar;
% Lasso priors
betanew = betabar(1:km1);
lambdabar = (a0+1)./(b0 + (taubar)/2);
kbetanew = diag(Kbetabar(1:km1,1:km1));
[taunew,entropybeta,logVbeta] = lasso1(km1,lambdabar,betanew,kbetanew);
taubar = 1./taunew(:,:);
invVbetanew = diag(taunew);
invVbeta(1:km1,1:km1) = invVbetanew;
invVbeta2 = invVbeta(km1+1:end,km1+1:end);
beta2 = betabar(km1+1:end);

Kbetabar = (Kbetabar + Kbetabar')/2;
Kbetabar1 = Kbetabar(1:km1,1:km1);
Kbetabar2 = Kbetabar(km1+1:end,km1+1:end);
gg = eig(full(Kbetabar))>0;
if sum(gg) == km 
CholKbeta = real(logdet(Kbetabar));
else
   CholKbeta = 0;
end
lliketau = sum(log(lambdabar/2) - lambdabar/2.*taubar) - sum(entropybeta);
llikelam = sum(lgammpdf(lambdabar,a0,b0) - gamment((a0+1),(b0 + (taubar)/2)));


vi = 1 + ii/2;
Si = s2(ii)/2;
invsig2bar = (vi + T/2)/Sbar;
llike1 = -T/2*log(2*pi) - T/2*log(Sbar/(nu0 + T/2)) - ((y(:,ii) - zt*betabar)'*(y(:,ii) - zt*betabar) + trace(zt*Kbetabar*zt'))/(2*Sbar/(nu0 + T/2)); % llikehood
llike2 = -km1/2*log(2*pi)  - .5*sum(logVbeta) - 0.5*(betanew'*invVbetanew*betanew + trace(invVbetanew*Kbetabar1)) ; % prior on theta
llike22 = -(iv-1)/2*log(2*pi) - 1/2*sum(log(1./(diag(invVbeta2*invsig2bar)))) - 0.5*(beta2'*invVbeta2*(invsig2bar)*beta2 + trace(invVbeta2*(invsig2bar)*Kbetabar2)) ; % prior on theta
llike3 = linvgammpdf((Sbar/(nu0 + T/2)),nu0,S0(ii)); % prior on the sigma
llike4 = km/2 + km/2*log(2*pi) + 0.5*CholKbeta; % entropy on theta
llike5 = invgamment(nu0 + T/2,Sbar); % entropy on sigma
llike = llike1 + llike2 + llike22 + llike3 - ( llike4 + llike5) + lliketau +  llikelam ;  
llike = llike/(T*n);


Kbetabar_store{ii,1} = Kbetabar;
betabar_store{ii,1} = betabar;
Sbar_store(ii,1) = Sbar;    
llike_store(ii,1) = llike;     
invVbeta_store{ii,1} = invVbeta ;
taubar_store{ii,1} = taubar ;
lambdabar_store{ii,1} =  lambdabar;

diff = abs(llike - pyhat);


if diff < 10^-4
    check = -10;
    store_diff(ii,:) = diff;
else
    pyhat = llike;

   
end

end


end

end

betat = [];
at = [];

for iv = 1:n
    if iv ==1
    indic =  Ridn{iv,1} + 1;     
    km = length(indic);
    theta2 = betabar_store{1,1} ;
%     betat =[betat;theta2];  
    betabar = zeros(k/n,1);
    betabar(indic,:) = theta2(1:length(indic));
    betat =[betat;betabar]; 
    else
    indic =  Ridn{iv,1} + 1;     
    km = length(indic) + iv-1;
    theta2 = betabar_store{iv,1};
    betabar = zeros(k/n,1);
    betabar(indic,:) = theta2(1:length(indic));
    betat =[betat;betabar]; 
    L(iv,1:iv-1) = theta2(length(indic)+1:end)';
    end
    
end
sig2 = Sbar_store./(nu0 + T/2 - 1);
invSig = (L'*diag(1./sig2)*L);


betacoeff = reshape(betat,k/n,n)';
betacoeff = [ L\betacoeff(:,1) L\betacoeff(:,2:n+1)];
[Y_new,munew,Kvar]=Sample_latent_Y_VBapproxcs1test(invSig,betacoeff,A,vecY,M_u,M_o, M_a, Y_con,p,T,n,cserr,idn);

YY = Y_new';
y0 = YY(1:p,:);
y = YY(p+1:end,:);
[T,n] = size(y);
longy = reshape(y',T*n,1);

X = zeros(T,np); 
tempY = [y0; y];
for i=1:p
    X(:,(i-1)*n+1:i*n) = tempY(p-i+1:end-i,:);
end
X = [ones(T,1) X];

Xqlag = zeros(T,nq*p); 
for i=1:p
    Xqlag(:,(i-1)*nq+1:i*nq) = tempY(p-i+1:end-i,1:nq);
end

 diff3 = mean(mean((y(:,nq+1:end)- pyhat0).^2));
% if count == 50 | diff3 < 10^-3
%     chck1 = 10;
% else
%     count = count + 1;
%     pyhat0 = y(:,nq+1:end) ;
% 
% end
%  diff1 = abs(VBlike-pyhat1);
if count>0
    if count == 50 | diff3 < 10^-3
        chck1 = 10;
    end
else
    count = count + 1;
    pyhat0 = y(:,nq+1:end) ;            
end
 

%% CS restriction
err = YY(:,1) - sum(ILT1w.*YY(:,6:17),2);
S0err1 = S0cs + 0.5*sum(err.^2);
cserr = S0err1/(nucs + T/2);

end
 toc
disp(['icount',num2str(icount)]);

Kvar = reshape(full(diag(Kvar\speye(length(Kvar)))),na,T+p)';


for sim = 1:1000
YYann = YY(Tint:end,nq+1:end);
YYannhat = YYann + sqrt(Kvar(Tint:end,:)).*randn(size(YYann,1),na);
%% Forecasting
CholSig = chol(invSig\speye(n),'lower');
betahat = reshape(betacoeff',k,1);
YYhat = [YY(Tint:end,1:nq) YYannhat];
for tt = 1:6
Xhat = kron(speye(n),[1 YYhat(end,:)]);
yhat = Xhat*betahat + CholSig*randn(n,1);
% if FACTOR==1
% if tt==1 && sR==1
%     yhat(6:5+12*nf)=FMatrix(T+1,:);
% end
% end
YYhat = [YYhat;yhat'];
end

yhatann = [];

for tv = 7:size(YYhat,1)
       yhatann = [ yhatann ;.25*YYhat(tv,nq+1:end) + .5*YYhat(tv-1,nq+1:end) + .75*YYhat(tv-2,nq+1:end) + YYhat(tv-3,nq+1:end) + .75*YYhat(tv-4,nq+1:end) + .5*YYhat(tv-5,nq+1:end) + .25*YYhat(tv-6,nq+1:end)];
    
end

fdates = 2013:.25:enddate+.25*6;
fdates = fdates(7:end);

store_nowcastITL1 = [store_nowcastITL1 ;yhatann(find(fdates==ITL1fcast(icount,2)),1:12)];
store_backcastITL1 = [store_backcastITL1;yhatann(find(fdates==ITL1fcast(icount,3)),1:12)];
store_forecastITL1 = [store_forecastITL1;yhatann(find(fdates==ITL1fcast(icount,1)),1:12)];

end

if sR == 1
if icount<5
elseif icount>94
else
  Q1_storeITL1{Q1countITL1,1} = (mean(store_nowcastITL1) - RSTIann(end-2,:)).^2;
  Q1_storeITL1{Q1countITL1,2} = crps(store_nowcastITL1, RSTIann(end-2,:));  
  store_nowcast_error = [store_nowcast_error;mean(store_nowcastITL1 - RSTIann(end-2,:))];
end

if icount<2
elseif icount>91
else
  Q1_storeITL1{Q1countITL1,5} = (mean(store_forecastITL1) - RSTIann(end-1,:)).^2;
  Q1_storeITL1{Q1countITL1,6} = crps(store_forecastITL1, RSTIann(end-1,:));
  store_forecast_error = [store_forecast_error;mean(store_forecastITL1 - RSTIann(end-1,:))];
end
Q1countITL1 =  Q1countITL1 + 1;


elseif sR==2

  if icount<5
elseif icount>94
  else
      Q2_storeITL1{Q2countITL1,1} = (mean(store_nowcastITL1) - RSTIann(end-1,:)).^2;
      Q2_storeITL1{Q2countITL1,2} = crps(store_nowcastITL1, RSTIann(end-1,:));
      store_nowcast_error = [store_nowcast_error;mean(store_nowcastITL1 - RSTIann(end-1,:))];
  end

  if icount<2
elseif icount>91
  else
      Q2_storeITL1{Q2countITL1,5} = (mean(store_forecastITL1) - RSTIann(end,:)).^2;
      Q2_storeITL1{Q2countITL1,6} = crps(store_forecastITL1, RSTIann(end,:));
      store_forecast_error = [store_forecast_error;mean(store_forecastITL1 - RSTIann(end,:))];
  end

  if icount<8
  else
  Q2_storeITL1{Q2countITL1,3} = (mean(store_backcastITL1) - RSTIann(end-2,:)).^2;
  Q2_storeITL1{Q2countITL1,4} = crps(store_backcastITL1, RSTIann(end-2,:));
  store_backcast_error = [store_backcast_error;mean(store_backcastITL1 - RSTIann(end-2,:))];
  end
  Q2countITL1 =  Q2countITL1 + 1;


else

  if icount<5
elseif icount>94
  else
      Q3_storeITL1{Q3countITL1,1} = (mean(store_nowcastITL1) - RSTIann(end-1,:)).^2;
      Q3_storeITL1{Q3countITL1,2} = crps(store_nowcastITL1, RSTIann(end-1,:));
      store_nowcast_error = [store_nowcast_error;mean(store_nowcastITL1 - RSTIann(end-1,:))];
  end

  if icount<2
elseif icount>91
  else
      Q3_storeITL1{Q3countITL1,5} = (mean(store_forecastITL1) - RSTIann(end,:)).^2;
      Q3_storeITL1{Q3countITL1,6} = crps(store_forecastITL1, RSTIann(end,:));
      store_forecast_error = [store_forecast_error;mean(store_forecastITL1 - RSTIann(end,:))];
  end
  if icount<8
  else
  Q3_storeITL1{Q3countITL1,3} = (mean(store_backcastITL1) - RSTIann(end-2,:)).^2;
  Q3_storeITL1{Q3countITL1,4} = crps(store_backcastITL1, RSTIann(end-2,:));
  store_backcast_error = [store_backcast_error;mean(store_backcastITL1 - RSTIann(end-2,:))];
  end
  Q3countITL1 =  Q3countITL1 + 1;


end

 icount = icount  + 1;

end

if icount==96+1

    break

end


end
