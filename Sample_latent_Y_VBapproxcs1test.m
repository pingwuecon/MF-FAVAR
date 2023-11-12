function [Y_new,munew,Knew]=Sample_latent_Y_VBapproxcs1test(invSig,beta,A,vecY,M_u,M_o, M_a, Y_con,p,T,n,cserr,idn)
    Tq = size(M_a,2);
    C=Compute_C(beta(:,2:end),A,p);
    Cu=C*M_u;
    bigK=Cu'*kron(speye(T),invSig);
    K=bigK*Cu; 
    K(1:12*7,1:12*7) = K(1:12*7,1:12*7) + 100*speye(12*7);  
    Kmu = bigK*(repmat(beta(:,1),T,1)-C*M_o*vecY);  
    invW = 1e10*ones(Tq,1);
    invW(find(idn==1)) = 1/cserr;
    iW = sparse(1:Tq,1:Tq,invW);
    Knew = K + M_a*iW*M_a';
    munew = Knew\(Kmu + M_a*iW*Y_con);    
    Y_new=reshape(M_o*vecY+M_u*munew,n,T+p);
