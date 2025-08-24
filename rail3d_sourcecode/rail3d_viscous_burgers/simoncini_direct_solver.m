function [xfinal,X] = simoncini_direct_solver(A1,A2,A3,M1,M2,H1,H3,B,check_res)
% 
% V. Simoncini, "Numerical solution of a class of third order tensor linear
% equations" Bollettino dell'Unione Matematica Italiana (BUMI), 2020.
% Theorem 2.1, Algorithm T3-sylv, extended to righthand sides in Tucker format.
% 
% Make sure Tensorlab Toolbox is in the path.
% 
% 
% Assume B (rhs = f = vec(B)) is in a general NxNxN tensor format.
% Solve G * xfinal -f=0
% G = kron(kron(M1,A1),H1)+kron(kron(A2,M2),H1)+kron(kron(H3,M2),A3) ]
% f = vec(B)
%
% All coeff matrices are nxn
%
% X = tensor(xfinal,n,n,n)
%
% Uses MATLAB's lyap function, contained in Control Systems Toolbox

[Q,R]=schur(full( A3'/H1' ));

bb = tens2mat(B,1); % bb = B_(1), unfolding of B in first dimension
%vecB = real(bb(:)); %commented so it does not display
%nrmrhs=norm(vecB); %commented so it does not display

g = (bb'/H1')*Q; % size r2*r3 x ng
[~,ng] = size(g);
[n1,n2]=size(A1);

for k=1:ng
    G = reshape(g(:,k),n2,n2); % g_k = vec(G_k)
    F = M2\G/M1';
    if k>1
        W = reshape(zz(:,1:k-1)*R(1:k-1,k),n1,n2);
        F = F - (W*H3')/M1';
    end
    Zhat = lyap(M2\A1, (R(k,k)*H3' + A2')/M1', -F);
    %Zhat = sylvester(M2\A1, (R(k,k)*H3' + A2')/M1', F);
    zz(:,k)=reshape(Zhat,n1*n2,1);
end

xx=Q*zz';
xfinal=real(xx(:));
X=reshape(real(xx),n1,n2,n2);

if check_res
    %L=kron(kron(M1,A1),H1)+kron(kron(A2,M2),H1)+kron(kron(H3,M2),A3); %commented so it does not display
    %normres=norm(L*xfinal- vecB)/nrmrhs; %commented so it does not display
    %fprintf('final full relative residual: %d\n',normres); %commented so it does not display
end
end












