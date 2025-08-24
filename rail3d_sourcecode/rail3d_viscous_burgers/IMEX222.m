function [U_nn,G_nn,MLR_nn] = IMEX222(U_n,G_n,MLR_n,A,B,C,P,tn,dtn,Dx,Dy,Dz,Dxx,Dyy,Dzz,Nx,Ny,Nz,tol,w1,w2,w3,dx,dy,dz,rhoM,JxM,JyM,JzM,kM,c,cc,xvals,yvals,zvals)
gamma = 1-1/sqrt(2);
delta = 1-1/(2*gamma);

% Compute flux terms (at time t^(0) = t^{n})
E1_n = {{TKR(U_n{1},U_n{1}),TKR(U_n{2},U_n{2}),TKR(U_n{3},U_n{3})},tkron(0.5*G_n,G_n)}; %E1 = A*U
E2_n = {{TKR(U_n{1},U_n{1}),TKR(U_n{2},U_n{2}),TKR(U_n{3},U_n{3})},tkron(0.5*G_n,G_n)}; %E2 = B*U
E3_n = {{TKR(U_n{1},U_n{1}),TKR(U_n{2},U_n{2}),TKR(U_n{3},U_n{3})},tkron(0.5*G_n,G_n)}; %E3 = C*U
% Compute source term (at time t^(1) = t^n + gamma*dt)
P_1 = P(tn+gamma*dtn);
% Pull out bases
Vx_n = U_n{1};
Vy_n = U_n{2};
Vz_n = U_n{3};
S_n = G_n;

% % % STAGE 1, t^(1) = t^n + gamma*dt % % %
[U_1,G_1,~] = IMEX111(U_n,G_n,MLR_n,A,B,C,P,tn,gamma*dtn,Dx,Dy,Dz,Dxx,Dyy,Dzz,Nx,Ny,Nz,tol,w1,w2,w3,dx,dy,dz,rhoM,JxM,JyM,JzM,kM,c,cc,xvals,yvals,zvals);
% Compute flux terms (at time t^(1) = t^{n} + gamma*dt)
E1_1 = {{TKR(U_1{1},U_1{1}),TKR(U_1{2},U_1{2}),TKR(U_1{3},U_1{3})},tkron(0.5*G_1,G_1)}; %E1 = A*U
E2_1 = {{TKR(U_1{1},U_1{1}),TKR(U_1{2},U_1{2}),TKR(U_1{3},U_1{3})},tkron(0.5*G_1,G_1)}; %E2 = B*U
E3_1 = {{TKR(U_1{1},U_1{1}),TKR(U_1{2},U_1{2}),TKR(U_1{3},U_1{3})},tkron(0.5*G_1,G_1)}; %E3 = C*U
% Compute source term (at time t^(2) = t^{n+1})
P_nn = P(tn+dtn);
% Pull out bases
Vx_1 = U_1{1};
Vy_1 = U_1{2};
Vz_1 = U_1{3};
S_1 = G_1;


% % % STAGE 2, t^(2) = t^{n+1} % % %
% K Steps
[U_dagger,~,~] = IMEX111(U_n,G_n,MLR_n,A,B,C,P,tn,dtn,Dx,Dy,Dz,Dxx,Dyy,Dzz,Nx,Ny,Nz,tol,w1,w2,w3,dx,dy,dz,rhoM,JxM,JyM,JzM,kM,c,cc,xvals,yvals,zvals);
[Vx_star,Vy_star,Vz_star,r1,r2,r3] = red_aug_K(U_dagger{1},[Vx_1,Vx_n],U_dagger{2},[Vy_1,Vy_n],U_dagger{3},[Vz_1,Vz_n]);

% K1 step (x basis update)
K1 = sylvester( full(speye(Nx,Nx) - (gamma*dtn)*Dxx) ,...
    full( -(gamma*dtn)*(kron( (Dzz*Vz_star)'*Vz_star , speye(r2,r2) ) + kron( speye(r3,r3) , (Dyy*Vy_star)'*Vy_star )) ) ,...
    Vx_n*tens2mat(S_n,1)*kron( Vz_n'*Vz_star , Vy_n'*Vy_star )...
    + ((1-gamma)*dtn)*(  ((Dxx*Vx_1)*tens2mat(S_1,1)*kron( Vz_1'*Vz_star , Vy_1'*Vy_star ))...
                        +(Vx_1*tens2mat(S_1,1)*kron( Vz_1'*Vz_star , (Dyy*Vy_1)'*Vy_star ))...
                        +(Vx_1*tens2mat(S_1,1)*kron( (Dzz*Vz_1)'*Vz_star , Vy_1'*Vy_star ))...
                        +(P_1{1}{1}*tens2mat(P_1{2},1)*kron( P_1{1}{3}'*Vz_star , P_1{1}{2}'*Vy_star ))  )...
    - delta*dtn*(  ( (Dx*E1_n{1}{1})*tens2mat(E1_n{2},1)*kron( E1_n{1}{3}'*Vz_star , E1_n{1}{2}'*Vy_star ) )...
                  +( E2_n{1}{1}*tens2mat(E2_n{2},1)*kron( E2_n{1}{3}'*Vz_star , (Dy*E2_n{1}{2})'*Vy_star ) )...
                  +( E3_n{1}{1}*tens2mat(E3_n{2},1)*kron( (Dz*E3_n{1}{3})'*Vz_star , E3_n{1}{2}'*Vy_star ) )  )...
    - (1-delta)*dtn*(  ( (Dx*E1_1{1}{1})*tens2mat(E1_1{2},1)*kron( E1_1{1}{3}'*Vz_star , E1_1{1}{2}'*Vy_star ) )...
                      +( E2_1{1}{1}*tens2mat(E2_1{2},1)*kron( E2_1{1}{3}'*Vz_star , (Dy*E2_1{1}{2})'*Vy_star ) )...
                      +( E3_1{1}{1}*tens2mat(E3_1{2},1)*kron( (Dz*E3_1{1}{3})'*Vz_star , E3_1{1}{2}'*Vy_star ) )  )...
    + gamma*dtn*( P_nn{1}{1}*tens2mat(P_nn{2},1)*kron( P_nn{1}{3}'*Vz_star , P_nn{1}{2}'*Vy_star ) ) );
[Vx_ddagger,~] = qr(K1,0);

% K2 step (y basis update)
K2 = sylvester( full(speye(Ny,Ny) - (gamma*dtn)*Dyy) ,...
    full( -(gamma*dtn)*(kron( (Dzz*Vz_star)'*Vz_star , speye(r1,r1) ) + kron( speye(r3,r3) , (Dxx*Vx_star)'*Vx_star )) ) ,...
    Vy_n*tens2mat(S_n,2)*kron( Vz_n'*Vz_star , Vx_n'*Vx_star )...
    + ((1-gamma)*dtn)*(  (Vy_1*tens2mat(S_1,2)*kron( Vz_1'*Vz_star , (Dxx*Vx_1)'*Vx_star ))...
                        +((Dyy*Vy_1)*tens2mat(S_1,2)*kron( Vz_1'*Vz_star , Vx_1'*Vx_star ))...
                        +(Vy_1*tens2mat(S_1,2)*kron( (Dzz*Vz_1)'*Vz_star , Vx_1'*Vx_star ))...
                        +(P_1{1}{2}*tens2mat(P_1{2},2)*kron( P_1{1}{3}'*Vz_star , P_1{1}{1}'*Vx_star ))  )...
    - delta*dtn*(  ( E1_n{1}{2}*tens2mat(E1_n{2},2)*kron( E1_n{1}{3}'*Vz_star , (Dx*E1_n{1}{1})'*Vx_star ) )...
                  +( (Dy*E2_n{1}{2})*tens2mat(E2_n{2},2)*kron( E2_n{1}{3}'*Vz_star , E2_n{1}{1}'*Vx_star ) )...
                  +( E3_n{1}{2}*tens2mat(E3_n{2},2)*kron( (Dz*E3_n{1}{3})'*Vz_star , E3_n{1}{1}'*Vx_star ) )  )...
    - (1-delta)*dtn*(  ( E1_1{1}{2}*tens2mat(E1_1{2},2)*kron( E1_1{1}{3}'*Vz_star , (Dx*E1_1{1}{1})'*Vx_star ) )...
                      +( (Dy*E2_1{1}{2})*tens2mat(E2_1{2},2)*kron( E2_1{1}{3}'*Vz_star , E2_1{1}{1}'*Vx_star ) )...
                      +( E3_1{1}{2}*tens2mat(E3_1{2},2)*kron( (Dz*E3_1{1}{3})'*Vz_star , E3_1{1}{1}'*Vx_star ) )  )...
    + gamma*dtn*( P_nn{1}{2}*tens2mat(P_nn{2},2)*kron( P_nn{1}{3}'*Vz_star , P_nn{1}{1}'*Vx_star ) ) );
[Vy_ddagger,~] = qr(K2,0);

% K3 step (z basis update)
K3 = sylvester( full(speye(Nz,Nz) - (gamma*dtn)*Dzz) ,...
    full( -(gamma*dtn)*(kron( (Dyy*Vy_star)'*Vy_star , speye(r1,r1) ) + kron( speye(r2,r2) , (Dxx*Vx_star)'*Vx_star )) ) ,...
    Vz_n*tens2mat(S_n,3)*kron( Vy_n'*Vy_star , Vx_n'*Vx_star )...
    + ((1-gamma)*dtn)*(  (Vz_1*tens2mat(S_1,3)*kron( Vy_1'*Vy_star , (Dxx*Vx_1)'*Vx_star ))...
                        +(Vz_1*tens2mat(S_1,3)*kron( (Dyy*Vy_1)'*Vy_star , Vx_1'*Vx_star ))...
                        +((Dzz*Vz_1)*tens2mat(S_1,3)*kron( Vy_1'*Vy_star , Vx_1'*Vx_star ))...
                        +(P_1{1}{3}*tens2mat(P_1{2},3)*kron( P_1{1}{2}'*Vy_star , P_1{1}{1}'*Vx_star ))  )...
    - delta*dtn*(  ( E1_n{1}{3}*tens2mat(E1_n{2},3)*kron( E1_n{1}{2}'*Vy_star , (Dx*E1_n{1}{1})'*Vx_star ) )...
                  +( E2_n{1}{3}*tens2mat(E2_n{2},3)*kron( (Dy*E2_n{1}{2})'*Vy_star , E2_n{1}{1}'*Vx_star ) )...
                  +( (Dz*E3_n{1}{3})*tens2mat(E3_n{2},3)*kron( E3_n{1}{2}'*Vy_star , E3_n{1}{1}'*Vx_star ) )  )...
    - (1-delta)*dtn*(  ( E1_1{1}{3}*tens2mat(E1_1{2},3)*kron( E1_1{1}{2}'*Vy_star , (Dx*E1_1{1}{1})'*Vx_star ) )...
                      +( E2_1{1}{3}*tens2mat(E2_1{2},3)*kron( (Dy*E2_1{1}{2})'*Vy_star , E2_1{1}{1}'*Vx_star ) )...
                      +( (Dz*E3_1{1}{3})*tens2mat(E3_1{2},3)*kron( E3_1{1}{2}'*Vy_star , E3_1{1}{1}'*Vx_star ) )  )...
    + gamma*dtn*( P_nn{1}{3}*tens2mat(P_nn{2},3)*kron( P_nn{1}{2}'*Vy_star , P_nn{1}{1}'*Vx_star ) ) );
[Vz_ddagger,~] = qr(K3,0);


% S step (update transfer tensor)
% Using Simoncini's solver requires square coefficient matrices. So, we
% must make the dimensions r1=r2=r3=R for the S step. These dimensions
% could be different after truncation. Note that Simoncini's solver will
% output a "square" order-3 tensor. Hence, RxRxR is okay if R>1. However,
% if R=1, then instead of 1x1x1, MATLAB will instead output a 1x1 array.
% So, we must force a 1x1x1 tensor in this situation after solving for
% S_nn using Simoncini's solver.
[Vx_nn,Vy_nn,Vz_nn,R] = red_aug_S(Vx_ddagger,[Vx_1,Vx_n],Vy_ddagger,[Vy_1,Vy_n],Vz_ddagger,[Vz_1,Vz_n]);
% r1 = r2 = r3 = R after reduced augmentation. Could be different after
% truncation.

A1 = -(gamma*dtn)*Vy_nn'*(Dyy*Vy_nn);
A2 = speye(R,R)-(gamma*dtn)*Vz_nn'*(Dzz*Vz_nn);
A3 = -(gamma*dtn)*Vx_nn'*(Dxx*Vx_nn);
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

B1 = Vx_nn'*(Dxx*Vx_1);
B2 = Vy_nn'*Vy_1;
B3 = Vz_nn'*Vz_1;
S = ((1-gamma)*dtn)*S_1;
B = B + lmlragen({B1,B2,B3},S);

B1 = Vx_nn'*Vx_1;
B2 = Vy_nn'*(Dyy*Vy_1);
B3 = Vz_nn'*Vz_1;
S = ((1-gamma)*dtn)*S_1;
B = B + lmlragen({B1,B2,B3},S);

B1 = Vx_nn'*Vx_1;
B2 = Vy_nn'*Vy_1;
B3 = Vz_nn'*(Dzz*Vz_1);
S = ((1-gamma)*dtn)*S_1;
B = B + lmlragen({B1,B2,B3},S);

B1 = Vx_nn'*P_1{1}{1};
B2 = Vy_nn'*P_1{1}{2};
B3 = Vz_nn'*P_1{1}{3};
S = ((1-gamma)*dtn)*P_1{2};
B = B + lmlragen({B1,B2,B3},S);

B1 = Vx_nn'*(Dx*E1_n{1}{1});
B2 = Vy_nn'*E1_n{1}{2};
B3 = Vz_nn'*E1_n{1}{3};
S = -(delta*dtn)*E1_n{2};
B = B + lmlragen({B1,B2,B3},S);

B1 = Vx_nn'*E2_n{1}{1};
B2 = Vy_nn'*(Dy*E2_n{1}{2});
B3 = Vz_nn'*E2_n{1}{3};
S = -(delta*dtn)*E2_n{2};
B = B + lmlragen({B1,B2,B3},S);

B1 = Vx_nn'*E3_n{1}{1};
B2 = Vy_nn'*E3_n{1}{2};
B3 = Vz_nn'*(Dz*E3_n{1}{3});
S = -(delta*dtn)*E3_n{2};
B = B + lmlragen({B1,B2,B3},S);

B1 = Vx_nn'*(Dx*E1_1{1}{1});
B2 = Vy_nn'*E1_1{1}{2};
B3 = Vz_nn'*E1_1{1}{3};
S = -((1-delta)*dtn)*E1_1{2};
B = B + lmlragen({B1,B2,B3},S);

B1 = Vx_nn'*E2_1{1}{1};
B2 = Vy_nn'*(Dy*E2_1{1}{2});
B3 = Vz_nn'*E2_1{1}{3};
S = -((1-delta)*dtn)*E2_1{2};
B = B + lmlragen({B1,B2,B3},S);

B1 = Vx_nn'*E3_1{1}{1};
B2 = Vy_nn'*E3_1{1}{2};
B3 = Vz_nn'*(Dz*E3_1{1}{3});
S = -((1-delta)*dtn)*E3_1{2};
B = B + lmlragen({B1,B2,B3},S);

B1 = Vx_nn'*P_nn{1}{1};
B2 = Vy_nn'*P_nn{1}{2};
B3 = Vz_nn'*P_nn{1}{3};
S = (gamma*dtn)*P_nn{2};
B = B + lmlragen({B1,B2,B3},S);


[~,S_nn] = simoncini_direct_solver(A1,A2,A3,M1,M2,H1,H3,B,check_res); %outputs double data type

% Truncation
% [U_nn,G_nn,MLR_nn] = nonconstrun({Vx_nn,Vy_nn,Vz_nn},S_nn,tol);
[U_nn,G_nn,MLR_nn] = lomac_0({Vx_nn,Vy_nn,Vz_nn},S_nn,tol,w1,w2,w3,Nx,Ny,Nz,dx,dy,dz,rhoM);
% [U_nn,G_nn,MLR_nn] = lomac_01({Vx_nn,Vy_nn,Vz_nn},S_nn,tol,w1,w2,w3,Nx,Ny,Nz,dx,dy,dz,rhoM,JxM,JyM,JzM,xvals,yvals,zvals);
% [U_nn,G_nn,MLR_nn] = lomac_012({Vx_nn,Vy_nn,Vz_nn},S_nn,tol,w1,w2,w3,Nx,Ny,Nz,dx,dy,dz,rhoM,JxM,JyM,JzM,kM,c,cc,xvals,yvals,zvals);
end











