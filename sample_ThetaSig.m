% This script obtains the posterior draws under the asymmetric conjugate prior
%
% See:
% Chan, J.C.C. (2019). Asymmetric conjugate priors for large Bayesian VARs,
% CAMA Working Papers 51/2019

function [betabar_store,invsig2_store] = sample_ThetaSig(Y0,Y,Z,p,kappa,Ridn)
[T,n] = size(Y);
sig2 = get_resid_var(Y0,Y);

betabar_store = cell(n,1);
invsig2_store = ones(n,1);
parfor ii = 1:n
    indic =  Ridn{ii,1} + 1;
    yi = Y(:,ii); 
    [mi,Vi,nui,Si] = prior_ACPi(n,p,ii,kappa,sig2);
    Xi = [Z(:,indic) -Y(:,1:ii-1)];
    ki = size(Xi,2);
    Vi = [Vi(indic);Vi(end-ii+2:end)];
    mi = [mi(indic);mi(end-ii+2:end)];
        % compute the parameters of the posterior distribution
    iVi = sparse(1:ki,1:ki,Vi)\speye(ki);
    Kthetai = iVi + Xi'*Xi;
    CKthetai = chol(Kthetai,'lower');    
    thetai_hat = (CKthetai')\(CKthetai\(iVi*mi + Xi'*yi));
    Si_hat = Si + (yi'*yi + mi'*iVi*mi - thetai_hat'*Kthetai*thetai_hat)/2;
        % sample sig and theta
    Sigi = 1./gamrnd(nui+T/2,1./Si_hat,1,1);
    U = randn(1,ki).*repmat(sqrt(Sigi),1,ki);
    Thetai = repmat(thetai_hat',1,1) + U/CKthetai;
    invsig2 =  1./Sigi;
    betabar_store{ii,1} = Thetai;
    invsig2_store(ii,1) = invsig2;
end   
end