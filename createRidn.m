Ridn1 = (0:5);
for i1=1:5
    Ridn{i1,1} = Ridn1;
%     Ridn{i1,1} = [Ridn1 Ridn1(2:end)+n Ridn1(2:end)+2*n Ridn1(2:end)+3*n Ridn1(2:end)+4*n Ridn1(2:end)+5*n Ridn1(2:end)+6*n];
end

for i1=1:12
    Ridn2 = (0:5);
if nf==0
    Ridn3 = [Ridn2 i1+5];
elseif nf==1
    Ridn3 = [Ridn2 i1+5 i1+5+12];
elseif nf==2
    Ridn3 = [Ridn2 i1+5 i1+5+12 i1+5+12+12];
elseif nf==3
    Ridn3 = [Ridn2 i1+5 i1+5+12 i1+5+12+12 i1+5+12+12+12];    
elseif nf==4
    Ridn3 = [Ridn2 i1+5 i1+5+12 i1+5+12+12 i1+5+12+12+12 i1+5+12+12+12+12]; 
elseif nf==5
    Ridn3 = [Ridn2 i1+5 i1+5+12 i1+5+12+12 i1+5+12+12+12 i1+5+12+12+12+12 i1+5+12+12+12+12+12];
end
if p==1
    Ridn{i1+5,1} = Ridn3;
elseif p==4
    Ridn{i1+5,1} = [Ridn3 Ridn3(2:end)+n Ridn3(2:end)+2*n Ridn3(2:end)+3*n];
elseif p==7
    Ridn{i1+5,1} = [Ridn3 Ridn3(2:end)+n Ridn3(2:end)+2*n Ridn3(2:end)+3*n Ridn3(2:end)+4*n Ridn3(2:end)+5*n Ridn3(2:end)+6*n];
end
end

if FACTOR==1
for i1=1:12
    Ridn2 = (0:5);
if nf==1
    Ridn3 = [Ridn2 i1+5 i1+5+12];
elseif nf==2
    Ridn3 = [Ridn2 i1+5 i1+5+12 i1+5+12+12];
elseif nf==3
    Ridn3 = [Ridn2 i1+5 i1+5+12 i1+5+12+12 i1+5+12+12+12];    
elseif nf==4
    Ridn3 = [Ridn2 i1+5 i1+5+12 i1+5+12+12 i1+5+12+12+12 i1+5+12+12+12+12]; 
elseif nf==5
    Ridn3 = [Ridn2 i1+5 i1+5+12 i1+5+12+12 i1+5+12+12+12 i1+5+12+12+12+12 i1+5+12+12+12+12+12]; 
end
if p==1
    Ridn{i1+5+12,1} = Ridn3;
elseif p==4
    Ridn{i1+5+12,1} = [Ridn3 Ridn3(2:end)+n Ridn3(2:end)+2*n Ridn3(2:end)+3*n];
elseif p==7
    Ridn{i1+5+12,1} = [Ridn3 Ridn3(2:end)+n Ridn3(2:end)+2*n Ridn3(2:end)+3*n Ridn3(2:end)+4*n Ridn3(2:end)+5*n Ridn3(2:end)+6*n];
end
end
% % 
if nf>1
for i1=1:12
    Ridn2 = (0:5);
if nf==1
    Ridn3 = [Ridn2 i1+5 i1+5+12];
elseif nf==2
    Ridn3 = [Ridn2 i1+5 i1+5+12 i1+5+12+12];
elseif nf==3
    Ridn3 = [Ridn2 i1+5 i1+5+12 i1+5+12+12 i1+5+12+12+12];    
elseif nf==4
    Ridn3 = [Ridn2 i1+5 i1+5+12 i1+5+12+12 i1+5+12+12+12 i1+5+12+12+12+12]; 
elseif nf==5
    Ridn3 = [Ridn2 i1+5 i1+5+12 i1+5+12+12 i1+5+12+12+12 i1+5+12+12+12+12 i1+5+12+12+12+12+12];
end
if p==1
    Ridn{i1+5+12*2,1} = Ridn3;
elseif p==4
    Ridn{i1+5+12*2,1} = [Ridn3 Ridn3(2:end)+n Ridn3(2:end)+2*n Ridn3(2:end)+3*n];
elseif p==7
    Ridn{i1+5+12+12,1} = [Ridn3 Ridn3(2:end)+n Ridn3(2:end)+2*n Ridn3(2:end)+3*n Ridn3(2:end)+4*n Ridn3(2:end)+5*n Ridn3(2:end)+6*n];
end
end
end
% % 
if nf>2
for i1=1:12
    Ridn2 = (0:5);
if nf==1
    Ridn3 = [Ridn2 i1+5 i1+5+12];
elseif nf==2
    Ridn3 = [Ridn2 i1+5 i1+5+12 i1+5+12+12];
elseif nf==3
    Ridn3 = [Ridn2 i1+5 i1+5+12 i1+5+12+12 i1+5+12+12+12];    
elseif nf==4
    Ridn3 = [Ridn2 i1+5 i1+5+12 i1+5+12+12 i1+5+12+12+12 i1+5+12+12+12+12]; 
elseif nf==5
    Ridn3 = [Ridn2 i1+5 i1+5+12 i1+5+12+12 i1+5+12+12+12 i1+5+12+12+12+12 i1+5+12+12+12+12+12]; 
end
if p==1
    Ridn{i1+5+12*3,1} = Ridn3;
elseif p==4
    Ridn{i1+5+12*3,1} = [Ridn3 Ridn3(2:end)+n Ridn3(2:end)+2*n Ridn3(2:end)+3*n];
elseif p==7
    Ridn{i1+5+12*3,1} = [Ridn3 Ridn3(2:end)+n Ridn3(2:end)+2*n Ridn3(2:end)+3*n Ridn3(2:end)+4*n Ridn3(2:end)+5*n Ridn3(2:end)+6*n];
end
end
end
% % 
if nf>3
for i1=1:12
    Ridn2 = (0:5);
if nf==1
    Ridn3 = [Ridn2 i1+5 i1+5+12];
elseif nf==2
    Ridn3 = [Ridn2 i1+5 i1+5+12 i1+5+12+12];
elseif nf==3
    Ridn3 = [Ridn2 i1+5 i1+5+12 i1+5+12+12 i1+5+12+12+12];    
elseif nf==4
    Ridn3 = [Ridn2 i1+5 i1+5+12 i1+5+12+12 i1+5+12+12+12 i1+5+12+12+12+12]; 
elseif nf==5
    Ridn3 = [Ridn2 i1+5 i1+5+12 i1+5+12+12 i1+5+12+12+12 i1+5+12+12+12+12 i1+5+12+12+12+12+12];
end
if p==1
    Ridn{i1+5+12*4,1} = Ridn3;
elseif p==4
    Ridn{i1+5+12*4,1} = [Ridn3 Ridn3(2:end)+n Ridn3(2:end)+2*n Ridn3(2:end)+3*n];
elseif p==7
    Ridn{i1+5+12*4,1} = [Ridn3 Ridn3(2:end)+n Ridn3(2:end)+2*n Ridn3(2:end)+3*n Ridn3(2:end)+4*n Ridn3(2:end)+5*n Ridn3(2:end)+6*n];
end
end
end
% % 
if nf>4
for i1=1:12
    Ridn2 = (0:5);
if nf==1
    Ridn3 = [Ridn2 i1+5 i1+5+12];
elseif nf==2
    Ridn3 = [Ridn2 i1+5 i1+5+12 i1+5+12+12];
elseif nf==3
    Ridn3 = [Ridn2 i1+5 i1+5+12 i1+5+12+12 i1+5+12+12+12];    
elseif nf==4
    Ridn3 = [Ridn2 i1+5 i1+5+12 i1+5+12+12 i1+5+12+12+12 i1+5+12+12+12+12]; 
elseif nf==5
    Ridn3 = [Ridn2 i1+5 i1+5+12 i1+5+12+12 i1+5+12+12+12 i1+5+12+12+12+12 i1+5+12+12+12+12+12];
end
if p==1
    Ridn{i1+5+12*5,1} = Ridn3;
elseif p==4
    Ridn{i1+5+12*5,1} = [Ridn3 Ridn3(2:end)+n Ridn3(2:end)+2*n Ridn3(2:end)+3*n];
elseif p==7
    Ridn{i1+5+12*5,1} = [Ridn3 Ridn3(2:end)+n Ridn3(2:end)+2*n Ridn3(2:end)+3*n Ridn3(2:end)+4*n Ridn3(2:end)+5*n Ridn3(2:end)+6*n];
end
end
end
end