function [F,corrx,FX] = Factor_est_vat(X,idtranlag,nf1,nf2,method,sR,ITL1id,VAT,vatlag)
% method: 1 EMPCA; 2 TW; 3 TP
F = cell(12,1); corrx = cell(12,1);
FX = cell(12,1);
for ivar=1:12
    pred1 = X{ivar};
    idtranlag1 = idtranlag{ivar};

    pred2 = pred1; 
    for j=1:length(idtranlag1)
        if 1<idtranlag1(j,2)/sR && idtranlag1(j,2)/sR<2
            pred2(end,j) = nan;
        elseif idtranlag1(j,2)/sR==2
            pred2(end-1:end,j) = nan;
        end
    end  
    T = length(pred2);

    vat4 = [];
    if ivar<11
        vat1 = VAT(ITL1id==ivar,:)';
        if vatlag==5
            if sR==1
                vat2 = vat1(1:end-2,:);
            else
                vat2 = vat1(1:end-1,:);
            end
        elseif vatlag==3
            if sR==1 || sR==2
                vat2 = vat1(1:end-1,:);
            else
                vat2 = vat1(1:end,:);
            end
        elseif vatlag==1
            vat2 = vat1(1:end,:);
        end
        vat3 = log(vat2(5:end,:)) - log(vat2(4:end-1,:));
        vat4 = [nan*ones(181,size(vat3,2));vat3;nan*ones(T-181-size(vat3,1),size(vat3,2))];
    end

    VARtmps = standardize_miss([pred2 vat4]);  
%     VARtmps = standardize_miss([pred2(:,1) pred2(:,5:6) vat4]);   
    VARtmps(:, all(isnan(VARtmps),1)) = [];
    if method==1
        [F3,FULLX,iter] = empca(VARtmps,nf1,nf2);
    elseif method==2
        [F3,FULLX] = TW3(VARtmps,nf1,1,nf2);
    elseif method==3
        [F3,FULLX] = TP(VARtmps,nf1,nf2); 
    end
    F{ivar} = F3;
%     F{ivar} = F3./100;
%     corrx{ivar} = corr([F3(:,1:3) VARtmps],'rows','complete');
end

