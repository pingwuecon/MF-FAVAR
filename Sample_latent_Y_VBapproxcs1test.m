function [Y_new,munew,Knew]=Sample_latent_Y_VBapproxcs1test(invSig,beta,A,vecY,M_u,M_o, M_a, Y_con,p,T,n,cserr,idn)
    Tq = size(M_a,2);
    C=Compute_C(beta(:,2:end),A,p);
    Cu=C*M_u;
    bigK=Cu'*kron(speye(T),invSig);
    K=bigK*Cu; 
%     K(1:12*7,1:12*7) = K(1:12*7,1:12*7) + 100*speye(12*7);  
    K(1:12*7,1:12*7) = K(1:12*7,1:12*7) + 100*speye(12*7);  
    Kmu = bigK*(repmat(beta(:,1),T,1)-C*M_o*vecY);  
    invW = 1e10*ones(Tq,1);
    invW(find(idn==1)) = 1/cserr;
    iW = sparse(1:Tq,1:Tq,invW);
    Knew = K + M_a*iW*M_a';
    munew = Knew\(Kmu + M_a*iW*Y_con);    
     Y_new=reshape(M_o*vecY+M_u*munew,n,T+p);
     mu=K\(bigK*(repmat(beta(:,1),T,1)-C*M_o*vecY));
     CMu= C*M_u;
     err = C*(M_o*vecY+M_u*mu) -   repmat(beta(:,1),T,1);
     err1 =  Y_con - M_a'*munew;
%      C2 = C;
%      bigK2 = C2'*kron(speye(T),invSig);
%      K2=bigK2*C2; 
%      K2(1:12*7,1:12*7) = K2(1:12*7,1:12*7) + 100*speye(12*7);  
%      Kmu2 = (bigK2*(repmat(beta(:,1),T,1)-C*M_o*vecY)); 
%      llike1 = -length(err)/2*log(2*pi) + T*sum(sum(log(diag(chol(invSig,'lower'))))) - .5*(err'*kron(speye(T),invSig)*err);
%      llike2 = -length(err1)/2*log(2*pi) -.5*sum(log((1./diag(iW)))) - .5*(err1'*iW*err1 + trace( M_a'*(Knew\M_a) ));
%      llike3 = length(err1)/2 + length(err1)/2*log(2*pi) - sum(log(diag(chol(Knew ,'lower'))));
%      llike = (llike1 + llike2) - llike3;
%      llike = 0;

%     C=Compute_C(beta(:,2:end),A,p);
%     Cu=C*M_u;
%     bigK=Cu'*sparse(blkdiag(speye(n)*1000,kron(speye(T-1),invSig)));
%     K=bigK*Cu;
%     mu=K\(bigK*(repmat(beta(:,1),T,1)-C*M_o*vecY));
%     diff=Y_con-M_a'*mu;
%     v=K\M_a;
%     munew=mu+v*((M_a'*v)\diff);
%     Y_new=reshape(M_o*vecY+M_u*munew,n,T+p);