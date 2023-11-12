function C=Compute_C(B,A,lag)
n=size(B,1);
C=kron(A{1},speye(n))-kron(A{2},B(:,1:n));
end  

