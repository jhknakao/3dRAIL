function X = TKR(A,B)
% Transpose Khatri-Rao product between matrices A and B.
% A = size N x r1
% B = size N x r2
% TKR(A,B) = size N x r1r2
[N1,r1] = size(A);
[N2,r2] = size(B);
if N1~=N2
    disp('A and B must have the same number of rows!')
    return
else
    N = N1;
end
X = zeros(N,r1*r2);
for k = 1:N
    X(k,:) = kron(A(k,:),B(k,:));
end
end
