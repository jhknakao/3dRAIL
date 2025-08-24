function [U_nn,G_nn,MLR_nn] = IMEX111(U_n,G_n,MLR_n,A,B,C,P,tn,dtn,Dx,Dy,Dz,Dxx,Dyy,Dzz,Nx,Ny,Nz,tol,w1,w2,w3,dx,dy,dz,rhoM,JxM,JyM,JzM,kM,c,cc,xvals,yvals,zvals)

% Compute flux terms (at time t^(0) = t^{n})
E1_n = {{TKR(U_n{1},U_n{1}),TKR(U_n{2},U_n{2}),TKR(U_n{3},U_n{3})},tkron(0.5*G_n,G_n)}; %E1 = A*U
E2_n = {{TKR(U_n{1},U_n{1}),TKR(U_n{2},U_n{2}),TKR(U_n{3},U_n{3})},tkron(0.5*G_n,G_n)}; %E2 = B*U
E3_n = {{TKR(U_n{1},U_n{1}),TKR(U_n{2},U_n{2}),TKR(U_n{3},U_n{3})},tkron(0.5*G_n,G_n)}; %E3 = C*U
% Compute source term (at time t^(1)=t^{n+1})
P_nn = P(tn+dtn);
% Pull out bases
Vx_n = U_n{1};
Vy_n = U_n{2};
Vz_n = U_n{3};
S_n = G_n;
r1_n = MLR_n(1);
r2_n = MLR_n(2);
r3_n = MLR_n(3);
Vx_star = Vx_n;
Vy_star = Vy_n;
Vz_star = Vz_n;

% K1 step (x basis update)
K1 = sylvester( full(speye(Nx,Nx) - dtn*Dxx) ,...
    full( -dtn*(kron( (Dzz*Vz_star)'*Vz_star , speye(r2_n,r2_n) ) + kron( speye(r3_n,r3_n) , (Dyy*Vy_star)'*Vy_star )) ) ,...
    Vx_n*tens2mat(S_n,1)*kron( Vz_n'*Vz_star , Vy_n'*Vy_star )...
    - dtn*(  ( (Dx*E1_n{1}{1})*tens2mat(E1_n{2},1)*kron( E1_n{1}{3}'*Vz_star , E1_n{1}{2}'*Vy_star ) )...
            +( E2_n{1}{1}*tens2mat(E2_n{2},1)*kron( E2_n{1}{3}'*Vz_star , (Dy*E2_n{1}{2})'*Vy_star ) )...
            +( E3_n{1}{1}*tens2mat(E3_n{2},1)*kron( (Dz*E3_n{1}{3})'*Vz_star , E3_n{1}{2}'*Vy_star ) )   )...
    + dtn*( P_nn{1}{1}*tens2mat(P_nn{2},1)*kron( P_nn{1}{3}'*Vz_star , P_nn{1}{2}'*Vy_star ) ) );
[Vx_ddagger,~] = qr(K1,0);

% K2 step (y basis update)
K2 = sylvester( full(speye(Ny,Ny) - dtn*Dyy) ,...
    full( -dtn*(kron( (Dzz*Vz_star)'*Vz_star , speye(r1_n,r1_n) ) + kron( speye(r3_n,r3_n) , (Dxx*Vx_star)'*Vx_star )) ) ,...
    Vy_n*tens2mat(S_n,2)*kron( Vz_n'*Vz_star , Vx_n'*Vx_star )...
    - dtn*(  ( E1_n{1}{2}*tens2mat(E1_n{2},2)*kron( E1_n{1}{3}'*Vz_star , (Dx*E1_n{1}{1})'*Vx_star ) )...
            +( (Dy*E2_n{1}{2})*tens2mat(E2_n{2},2)*kron( E2_n{1}{3}'*Vz_star , E2_n{1}{1}'*Vx_star ) )...
            +( E3_n{1}{2}*tens2mat(E3_n{2},2)*kron( (Dz*E3_n{1}{3})'*Vz_star , E3_n{1}{1}'*Vx_star ) )   )...
    + dtn*( P_nn{1}{2}*tens2mat(P_nn{2},2)*kron( P_nn{1}{3}'*Vz_star , P_nn{1}{1}'*Vx_star ) ) );
[Vy_ddagger,~] = qr(K2,0);

% K3 step (z basis update)
K3 = sylvester( full(speye(Nz,Nz) - dtn*Dzz) ,...
    full( -dtn*(kron( (Dyy*Vy_star)'*Vy_star , speye(r1_n,r1_n) ) + kron( speye(r2_n,r2_n) , (Dxx*Vx_star)'*Vx_star )) ) ,...
    Vz_n*tens2mat(S_n,3)*kron( Vy_n'*Vy_star , Vx_n'*Vx_star )...
    - dtn*(  ( E1_n{1}{3}*tens2mat(E1_n{2},3)*kron( E1_n{1}{2}'*Vy_star , (Dx*E1_n{1}{1})'*Vx_star ) )...
            +( E2_n{1}{3}*tens2mat(E2_n{2},3)*kron( (Dy*E2_n{1}{2})'*Vy_star , E2_n{1}{1}'*Vx_star ) )...
            +( (Dz*E3_n{1}{3})*tens2mat(E3_n{2},3)*kron( E3_n{1}{2}'*Vy_star , E3_n{1}{1}'*Vx_star ) )   )...
    + dtn*( P_nn{1}{3}*tens2mat(P_nn{2},3)*kron( P_nn{1}{2}'*Vy_star , P_nn{1}{1}'*Vx_star ) ) );
[Vz_ddagger,~] = qr(K3,0);

% S step (update transfer tensor)
% Using Simoncini's solver requires square coefficient matrices. So, we
% must make the dimensions r1=r2=r3=R for the S step. These dimensions
% could be different after truncation. Note that Simoncini's solver will
% output a "square" order-3 tensor. Hence, RxRxR is okay if R>1. However,
% if R=1, then instead of 1x1x1, MATLAB will instead output a 1x1 array.
% So, we must force a 1x1x1 tensor in this situation after solving for
% S_nn using Simoncini's solver.
[Vx_nn,Vy_nn,Vz_nn,R] = red_aug_S(Vx_ddagger,Vx_n,Vy_ddagger,Vy_n,Vz_ddagger,Vz_n);
% r1 = r2 = r3 = R after reduced augmentation. Could be different after
% truncation.

A1 = -dtn*Vy_nn'*(Dyy*Vy_nn);
A2 = speye(R,R)-dtn*Vz_nn'*(Dzz*Vz_nn);
A3 = -dtn*Vx_nn'*(Dxx*Vx_nn);
M1 = speye(R,R); %r3xr3
M2 = speye(R,R); %r2xr2
H1 = speye(R,R); %r1xr1
H3 = speye(R,R); %r3xr3
check_res = 1;

B1 = Vx_nn'*Vx_n;
B2 = Vy_nn'*Vy_n;
B3 = Vz_nn'*Vz_n;
S = S_n;
B = lmlragen({B1,B2,B3},S);

B1 = Vx_nn'*(Dx*E1_n{1}{1});
B2 = Vy_nn'*E1_n{1}{2};
B3 = Vz_nn'*E1_n{1}{3};
S = -dtn*E1_n{2};
B = B + lmlragen({B1,B2,B3},S);

B1 = Vx_nn'*E2_n{1}{1};
B2 = Vy_nn'*(Dy*E2_n{1}{2});
B3 = Vz_nn'*E2_n{1}{3};
S = -dtn*E2_n{2};
B = B + lmlragen({B1,B2,B3},S);

B1 = Vx_nn'*E3_n{1}{1};
B2 = Vy_nn'*E3_n{1}{2};
B3 = Vz_nn'*(Dz*E3_n{1}{3});
S = -dtn*E3_n{2};
B = B + lmlragen({B1,B2,B3},S);

B1 = Vx_nn'*P_nn{1}{1};
B2 = Vy_nn'*P_nn{1}{2};
B3 = Vz_nn'*P_nn{1}{3};
S = dtn*P_nn{2};
B = B + lmlragen({B1,B2,B3},S);

[~,S_nn] = simoncini_direct_solver(A1,A2,A3,M1,M2,H1,H3,B,check_res); %outputs double data type

% Truncation
% [U_nn,G_nn,MLR_nn] = nonconstrun({Vx_nn,Vy_nn,Vz_nn},S_nn,tol);
[U_nn,G_nn,MLR_nn] = lomac_0({Vx_nn,Vy_nn,Vz_nn},S_nn,tol,w1,w2,w3,Nx,Ny,Nz,dx,dy,dz,rhoM);
% [U_nn,G_nn,MLR_nn] = lomac_01({Vx_nn,Vy_nn,Vz_nn},S_nn,tol,w1,w2,w3,Nx,Ny,Nz,dx,dy,dz,rhoM,JxM,JyM,JzM,xvals,yvals,zvals);
% [U_nn,G_nn,MLR_nn] = lomac_012({Vx_nn,Vy_nn,Vz_nn},S_nn,tol,w1,w2,w3,Nx,Ny,Nz,dx,dy,dz,rhoM,JxM,JyM,JzM,kM,c,cc,xvals,yvals,zvals);
end



