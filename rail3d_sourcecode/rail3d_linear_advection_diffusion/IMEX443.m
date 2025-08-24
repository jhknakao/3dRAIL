function [U_nn,G_nn,MLR_nn] = IMEX443(U_n,G_n,MLR_n,A,B,C,P,tn,dtn,Dx,Dy,Dz,Dxx,Dyy,Dzz,Nx,Ny,Nz,tol,w1,w2,w3,dx,dy,dz,rhoM,JxM,JyM,JzM,kM,c,cc,xvals,yvals,zvals)

% Compute flux terms (at time t^(0) = t^{n})
A_n = A(tn); %A(tn)
B_n = B(tn); %B(tn)
C_n = C(tn); %C(tn)
E1_n = {{TKR(A_n{1}{1},U_n{1}),TKR(A_n{1}{2},U_n{2}),TKR(A_n{1}{3},U_n{3})},tkron(A_n{2},G_n)}; %E1 = A*U
E2_n = {{TKR(B_n{1}{1},U_n{1}),TKR(B_n{1}{2},U_n{2}),TKR(B_n{1}{3},U_n{3})},tkron(B_n{2},G_n)}; %E2 = B*U
E3_n = {{TKR(C_n{1}{1},U_n{1}),TKR(C_n{1}{2},U_n{2}),TKR(C_n{1}{3},U_n{3})},tkron(C_n{2},G_n)}; %E3 = C*U
% Compute source term (at time t^(1) = t^n + dt/2)
P_1 = P(tn+dtn/2);
% Pull out bases
Vx_n = U_n{1};
Vy_n = U_n{2};
Vz_n = U_n{3};
S_n = G_n;

%%
% % % STAGE 1, t^(1) = t^n + dt/2 % % %
[U_1,G_1,~] = IMEX111(U_n,G_n,MLR_n,A,B,C,P,tn,dtn/2,Dx,Dy,Dz,Dxx,Dyy,Dzz,Nx,Ny,Nz,tol,w1,w2,w3,dx,dy,dz,rhoM,JxM,JyM,JzM,kM,c,cc,xvals,yvals,zvals);
% Compute flux terms (at time t^(1) = t^n + dt/2)
A_1 = A(tn+dtn/2); %A(t^(1))
B_1 = B(tn+dtn/2); %B(t^(1))
C_1 = C(tn+dtn/2); %C(t^(1))
E1_1 = {{TKR(A_1{1}{1},U_1{1}),TKR(A_1{1}{2},U_1{2}),TKR(A_1{1}{3},U_1{3})},tkron(A_1{2},G_1)}; %E1 = A*U
E2_1 = {{TKR(B_1{1}{1},U_1{1}),TKR(B_1{1}{2},U_1{2}),TKR(B_1{1}{3},U_1{3})},tkron(B_1{2},G_1)}; %E2 = B*U
E3_1 = {{TKR(C_1{1}{1},U_1{1}),TKR(C_1{1}{2},U_1{2}),TKR(C_1{1}{3},U_1{3})},tkron(C_1{2},G_1)}; %E3 = C*U
% Compute source term (at time t^(2) = t^n + 2dt/3)
P_2 = P(tn+(2/3)*dtn);
% Pull out bases
Vx_1 = U_1{1};
Vy_1 = U_1{2};
Vz_1 = U_1{3};
S_1 = G_1;


%%
% % % STAGE 2, t^(2) = t^n + 2dt/3 % % %
% K Steps
[U_dagger,~,~] = IMEX111(U_n,G_n,MLR_n,A,B,C,P,tn,(2/3)*dtn,Dx,Dy,Dz,Dxx,Dyy,Dzz,Nx,Ny,Nz,tol,w1,w2,w3,dx,dy,dz,rhoM,JxM,JyM,JzM,kM,c,cc,xvals,yvals,zvals);
[Vx_star,Vy_star,Vz_star,r1,r2,r3] = red_aug_K(U_dagger{1},[Vx_1,Vx_n],U_dagger{2},[Vy_1,Vy_n],U_dagger{3},[Vz_1,Vz_n]);

% K1 step (x basis update)
K1 = sylvester( full(speye(Nx,Nx) - (dtn/2)*Dxx) ,...
    full( -(dtn/2)*(kron( (Dzz*Vz_star)'*Vz_star , speye(r2,r2) ) + kron( speye(r3,r3) , (Dyy*Vy_star)'*Vy_star )) ) ,...
    Vx_n*tens2mat(S_n,1)*kron( Vz_n'*Vz_star , Vy_n'*Vy_star )...
    + (dtn/6)*(  ((Dxx*Vx_1)*tens2mat(S_1,1)*kron( Vz_1'*Vz_star , Vy_1'*Vy_star ))...
                        +(Vx_1*tens2mat(S_1,1)*kron( Vz_1'*Vz_star , (Dyy*Vy_1)'*Vy_star ))...
                        +(Vx_1*tens2mat(S_1,1)*kron( (Dzz*Vz_1)'*Vz_star , Vy_1'*Vy_star ))...
                        +(P_1{1}{1}*tens2mat(P_1{2},1)*kron( P_1{1}{3}'*Vz_star , P_1{1}{2}'*Vy_star ))  )...
    - (11/18)*dtn*(  ( (Dx*E1_n{1}{1})*tens2mat(E1_n{2},1)*kron( E1_n{1}{3}'*Vz_star , E1_n{1}{2}'*Vy_star ) )...
                  +( E2_n{1}{1}*tens2mat(E2_n{2},1)*kron( E2_n{1}{3}'*Vz_star , (Dy*E2_n{1}{2})'*Vy_star ) )...
                  +( E3_n{1}{1}*tens2mat(E3_n{2},1)*kron( (Dz*E3_n{1}{3})'*Vz_star , E3_n{1}{2}'*Vy_star ) )  )...
    - (1/18)*dtn*(  ( (Dx*E1_1{1}{1})*tens2mat(E1_1{2},1)*kron( E1_1{1}{3}'*Vz_star , E1_1{1}{2}'*Vy_star ) )...
                      +( E2_1{1}{1}*tens2mat(E2_1{2},1)*kron( E2_1{1}{3}'*Vz_star , (Dy*E2_1{1}{2})'*Vy_star ) )...
                      +( E3_1{1}{1}*tens2mat(E3_1{2},1)*kron( (Dz*E3_1{1}{3})'*Vz_star , E3_1{1}{2}'*Vy_star ) )  )...
    + (dtn/2)*( P_2{1}{1}*tens2mat(P_2{2},1)*kron( P_2{1}{3}'*Vz_star , P_2{1}{2}'*Vy_star ) ) );
[Vx_ddagger,~] = qr(K1,0);

% K2 step (y basis update)
K2 = sylvester( full(speye(Ny,Ny) - (dtn/2)*Dyy) ,...
    full( -(dtn/2)*(kron( (Dzz*Vz_star)'*Vz_star , speye(r1,r1) ) + kron( speye(r3,r3) , (Dxx*Vx_star)'*Vx_star )) ) ,...
    Vy_n*tens2mat(S_n,2)*kron( Vz_n'*Vz_star , Vx_n'*Vx_star )...
    + (dtn/6)*(  (Vy_1*tens2mat(S_1,2)*kron( Vz_1'*Vz_star , (Dxx*Vx_1)'*Vx_star ))...
                        +((Dyy*Vy_1)*tens2mat(S_1,2)*kron( Vz_1'*Vz_star , Vx_1'*Vx_star ))...
                        +(Vy_1*tens2mat(S_1,2)*kron( (Dzz*Vz_1)'*Vz_star , Vx_1'*Vx_star ))...
                        +(P_1{1}{2}*tens2mat(P_1{2},2)*kron( P_1{1}{3}'*Vz_star , P_1{1}{1}'*Vx_star ))  )...
    - (11/18)*dtn*(  ( E1_n{1}{2}*tens2mat(E1_n{2},2)*kron( E1_n{1}{3}'*Vz_star , (Dx*E1_n{1}{1})'*Vx_star ) )...
                  +( (Dy*E2_n{1}{2})*tens2mat(E2_n{2},2)*kron( E2_n{1}{3}'*Vz_star , E2_n{1}{1}'*Vx_star ) )...
                  +( E3_n{1}{2}*tens2mat(E3_n{2},2)*kron( (Dz*E3_n{1}{3})'*Vz_star , E3_n{1}{1}'*Vx_star ) )  )...
    - (1/18)*dtn*(  ( E1_1{1}{2}*tens2mat(E1_1{2},2)*kron( E1_1{1}{3}'*Vz_star , (Dx*E1_1{1}{1})'*Vx_star ) )...
                      +( (Dy*E2_1{1}{2})*tens2mat(E2_1{2},2)*kron( E2_1{1}{3}'*Vz_star , E2_1{1}{1}'*Vx_star ) )...
                      +( E3_1{1}{2}*tens2mat(E3_1{2},2)*kron( (Dz*E3_1{1}{3})'*Vz_star , E3_1{1}{1}'*Vx_star ) )  )...
    + (dtn/2)*( P_2{1}{2}*tens2mat(P_2{2},2)*kron( P_2{1}{3}'*Vz_star , P_2{1}{1}'*Vx_star ) ) );
[Vy_ddagger,~] = qr(K2,0);

% K3 step (z basis update)
K3 = sylvester( full(speye(Nz,Nz) - (dtn/2)*Dzz) ,...
    full( -(dtn/2)*(kron( (Dyy*Vy_star)'*Vy_star , speye(r1,r1) ) + kron( speye(r2,r2) , (Dxx*Vx_star)'*Vx_star )) ) ,...
    Vz_n*tens2mat(S_n,3)*kron( Vy_n'*Vy_star , Vx_n'*Vx_star )...
    + (dtn/6)*(  (Vz_1*tens2mat(S_1,3)*kron( Vy_1'*Vy_star , (Dxx*Vx_1)'*Vx_star ))...
                        +(Vz_1*tens2mat(S_1,3)*kron( (Dyy*Vy_1)'*Vy_star , Vx_1'*Vx_star ))...
                        +((Dzz*Vz_1)*tens2mat(S_1,3)*kron( Vy_1'*Vy_star , Vx_1'*Vx_star ))...
                        +(P_1{1}{3}*tens2mat(P_1{2},3)*kron( P_1{1}{2}'*Vy_star , P_1{1}{1}'*Vx_star ))  )...
    - (11/18)*dtn*(  ( E1_n{1}{3}*tens2mat(E1_n{2},3)*kron( E1_n{1}{2}'*Vy_star , (Dx*E1_n{1}{1})'*Vx_star ) )...
                  +( E2_n{1}{3}*tens2mat(E2_n{2},3)*kron( (Dy*E2_n{1}{2})'*Vy_star , E2_n{1}{1}'*Vx_star ) )...
                  +( (Dz*E3_n{1}{3})*tens2mat(E3_n{2},3)*kron( E3_n{1}{2}'*Vy_star , E3_n{1}{1}'*Vx_star ) )  )...
    - (1/18)*dtn*(  ( E1_1{1}{3}*tens2mat(E1_1{2},3)*kron( E1_1{1}{2}'*Vy_star , (Dx*E1_1{1}{1})'*Vx_star ) )...
                      +( E2_1{1}{3}*tens2mat(E2_1{2},3)*kron( (Dy*E2_1{1}{2})'*Vy_star , E2_1{1}{1}'*Vx_star ) )...
                      +( (Dz*E3_1{1}{3})*tens2mat(E3_1{2},3)*kron( E3_1{1}{2}'*Vy_star , E3_1{1}{1}'*Vx_star ) )  )...
    + (dtn/2)*( P_2{1}{3}*tens2mat(P_2{2},3)*kron( P_2{1}{2}'*Vy_star , P_2{1}{1}'*Vx_star ) ) );
[Vz_ddagger,~] = qr(K3,0);


% S step (update transfer tensor)
% Using Simoncini's solver requires square coefficient matrices. So, we
% must make the dimensions r1=r2=r3=R for the S step. These dimensions
% could be different after truncation. Note that Simoncini's solver will
% output a "square" order-3 tensor. Hence, RxRxR is okay if R>1. However,
% if R=1, then instead of 1x1x1, MATLAB will instead output a 1x1 array.
% So, we must force a 1x1x1 tensor in this situation after solving for
% S_nn using Simoncini's solver.
[Vx_2,Vy_2,Vz_2,R] = red_aug_S(Vx_ddagger,[Vx_1,Vx_n],Vy_ddagger,[Vy_1,Vy_n],Vz_ddagger,[Vz_1,Vz_n]);
% r1 = r2 = r3 = R after reduced augmentation. Could be different after
% truncation.

A1 = -(dtn/2)*Vy_2'*(Dyy*Vy_2);
A2 = speye(R,R)-(dtn/2)*Vz_2'*(Dzz*Vz_2);
A3 = -(dtn/2)*Vx_2'*(Dxx*Vx_2);
M1 = speye(R,R); %r3xr3
M2 = speye(R,R); %r2xr2
H1 = speye(R,R); %r1xr1
H3 = speye(R,R); %r3xr3
check_res = 1;


O1 = Vx_2'*Vx_n;
O2 = Vy_2'*Vy_n;
O3 = Vz_2'*Vz_n;
S = S_n;
O = lmlragen({O1,O2,O3},S);

O1 = Vx_2'*(Dxx*Vx_1);
O2 = Vy_2'*Vy_1;
O3 = Vz_2'*Vz_1;
S = (dtn/6)*S_1;
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_2'*Vx_1;
O2 = Vy_2'*(Dyy*Vy_1);
O3 = Vz_2'*Vz_1;
S = (dtn/6)*S_1;
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_2'*Vx_1;
O2 = Vy_2'*Vy_1;
O3 = Vz_2'*(Dzz*Vz_1);
S = (dtn/6)*S_1;
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_2'*P_1{1}{1};
O2 = Vy_2'*P_1{1}{2};
O3 = Vz_2'*P_1{1}{3};
S = (dtn/6)*P_1{2};
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_2'*(Dx*E1_n{1}{1});
O2 = Vy_2'*E1_n{1}{2};
O3 = Vz_2'*E1_n{1}{3};
S = -((11/18)*dtn)*E1_n{2};
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_2'*E2_n{1}{1};
O2 = Vy_2'*(Dy*E2_n{1}{2});
O3 = Vz_2'*E2_n{1}{3};
S = -((11/18)*dtn)*E2_n{2};
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_2'*E3_n{1}{1};
O2 = Vy_2'*E3_n{1}{2};
O3 = Vz_2'*(Dz*E3_n{1}{3});
S = -((11/18)*dtn)*E3_n{2};
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_2'*(Dx*E1_1{1}{1});
O2 = Vy_2'*E1_1{1}{2};
O3 = Vz_2'*E1_1{1}{3};
S = -((1/18)*dtn)*E1_1{2};
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_2'*E2_1{1}{1};
O2 = Vy_2'*(Dy*E2_1{1}{2});
O3 = Vz_2'*E2_1{1}{3};
S = -((1/18)*dtn)*E2_1{2};
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_2'*E3_1{1}{1};
O2 = Vy_2'*E3_1{1}{2};
O3 = Vz_2'*(Dz*E3_1{1}{3});
S = -((1/18)*dtn)*E3_1{2};
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_2'*P_2{1}{1};
O2 = Vy_2'*P_2{1}{2};
O3 = Vz_2'*P_2{1}{3};
S = (dtn/2)*P_2{2};
O = O + lmlragen({O1,O2,O3},S);

[~,S_2] = simoncini_direct_solver(A1,A2,A3,M1,M2,H1,H3,O,check_res); %outputs double data type

% Truncation
% [U_2,G_2,~] = nonconstrun({Vx_2,Vy_2,Vz_2},S_2,tol);
% [U_2,G_2,~] = lomac_0({Vx_2,Vy_2,Vz_2},S_2,tol,w1,w2,w3,Nx,Ny,Nz,dx,dy,dz,rhoM);
% [U_2,G_2,~] = lomac_01({Vx_2,Vy_2,Vz_2},S_2,tol,w1,w2,w3,Nx,Ny,Nz,dx,dy,dz,rhoM,JxM,JyM,JzM,xvals,yvals,zvals);
[U_2,G_2,~] = lomac_012({Vx_2,Vy_2,Vz_2},S_2,tol,w1,w2,w3,Nx,Ny,Nz,dx,dy,dz,rhoM,JxM,JyM,JzM,kM,c,cc,xvals,yvals,zvals);

% Compute flux terms (at time t^(2) = t^n + 2dt/3)
A_2 = A(tn+(2/3)*dtn); %A(t^(2))
B_2 = B(tn+(2/3)*dtn); %B(t^(2))
C_2 = C(tn+(2/3)*dtn); %C(t^(2))
E1_2 = {{TKR(A_2{1}{1},U_2{1}),TKR(A_2{1}{2},U_2{2}),TKR(A_2{1}{3},U_2{3})},tkron(A_2{2},G_2)}; %E1 = A*U
E2_2 = {{TKR(B_2{1}{1},U_2{1}),TKR(B_2{1}{2},U_2{2}),TKR(B_2{1}{3},U_2{3})},tkron(B_2{2},G_2)}; %E2 = B*U
E3_2 = {{TKR(C_2{1}{1},U_2{1}),TKR(C_2{1}{2},U_2{2}),TKR(C_2{1}{3},U_2{3})},tkron(C_2{2},G_2)}; %E3 = C*U
% Compute source term (at time t^(3) = t^n + dt/2)
P_3 = P_1;
% Pull out bases
Vx_2 = U_2{1};
Vy_2 = U_2{2};
Vz_2 = U_2{3};
S_2 = G_2;


%%
% % % STAGE 3, t^(3) = t^n + dt/2 % % %
% K Steps
[U_dagger,~,~] = IMEX111(U_n,G_n,MLR_n,A,B,C,P,tn,dtn/2,Dx,Dy,Dz,Dxx,Dyy,Dzz,Nx,Ny,Nz,tol,w1,w2,w3,dx,dy,dz,rhoM,JxM,JyM,JzM,kM,c,cc,xvals,yvals,zvals);
[Vx_star,Vy_star,Vz_star,r1,r2,r3] = red_aug_K(U_dagger{1},[Vx_2,Vx_1,Vx_n],U_dagger{2},[Vy_2,Vy_1,Vy_n],U_dagger{3},[Vz_2,Vz_1,Vz_n]);

% K1 step (x basis update)
K1 = sylvester( full(speye(Nx,Nx) - (dtn/2)*Dxx) ,...
    full( -(dtn/2)*(kron( (Dzz*Vz_star)'*Vz_star , speye(r2,r2) ) + kron( speye(r3,r3) , (Dyy*Vy_star)'*Vy_star )) ) ,...
    Vx_n*tens2mat(S_n,1)*kron( Vz_n'*Vz_star , Vy_n'*Vy_star )...
    + (-dtn/2)*(  ((Dxx*Vx_1)*tens2mat(S_1,1)*kron( Vz_1'*Vz_star , Vy_1'*Vy_star ))...
                        +(Vx_1*tens2mat(S_1,1)*kron( Vz_1'*Vz_star , (Dyy*Vy_1)'*Vy_star ))...
                        +(Vx_1*tens2mat(S_1,1)*kron( (Dzz*Vz_1)'*Vz_star , Vy_1'*Vy_star ))...
                        +(P_1{1}{1}*tens2mat(P_1{2},1)*kron( P_1{1}{3}'*Vz_star , P_1{1}{2}'*Vy_star ))  )...
    + (dtn/2)*(  ((Dxx*Vx_2)*tens2mat(S_2,1)*kron( Vz_2'*Vz_star , Vy_2'*Vy_star ))...
                        +(Vx_2*tens2mat(S_2,1)*kron( Vz_2'*Vz_star , (Dyy*Vy_2)'*Vy_star ))...
                        +(Vx_2*tens2mat(S_2,1)*kron( (Dzz*Vz_2)'*Vz_star , Vy_2'*Vy_star ))...
                        +(P_2{1}{1}*tens2mat(P_2{2},1)*kron( P_2{1}{3}'*Vz_star , P_2{1}{2}'*Vy_star ))  )...
    - (5/6)*dtn*(  ( (Dx*E1_n{1}{1})*tens2mat(E1_n{2},1)*kron( E1_n{1}{3}'*Vz_star , E1_n{1}{2}'*Vy_star ) )...
                  +( E2_n{1}{1}*tens2mat(E2_n{2},1)*kron( E2_n{1}{3}'*Vz_star , (Dy*E2_n{1}{2})'*Vy_star ) )...
                  +( E3_n{1}{1}*tens2mat(E3_n{2},1)*kron( (Dz*E3_n{1}{3})'*Vz_star , E3_n{1}{2}'*Vy_star ) )  )...
    - (-5/6)*dtn*(  ( (Dx*E1_1{1}{1})*tens2mat(E1_1{2},1)*kron( E1_1{1}{3}'*Vz_star , E1_1{1}{2}'*Vy_star ) )...
                      +( E2_1{1}{1}*tens2mat(E2_1{2},1)*kron( E2_1{1}{3}'*Vz_star , (Dy*E2_1{1}{2})'*Vy_star ) )...
                      +( E3_1{1}{1}*tens2mat(E3_1{2},1)*kron( (Dz*E3_1{1}{3})'*Vz_star , E3_1{1}{2}'*Vy_star ) )  )...
    - (1/2)*dtn*(  ( (Dx*E1_2{1}{1})*tens2mat(E1_2{2},1)*kron( E1_2{1}{3}'*Vz_star , E1_2{1}{2}'*Vy_star ) )...
                      +( E2_2{1}{1}*tens2mat(E2_2{2},1)*kron( E2_2{1}{3}'*Vz_star , (Dy*E2_2{1}{2})'*Vy_star ) )...
                      +( E3_2{1}{1}*tens2mat(E3_2{2},1)*kron( (Dz*E3_2{1}{3})'*Vz_star , E3_2{1}{2}'*Vy_star ) )  )...
    + (dtn/2)*( P_3{1}{1}*tens2mat(P_3{2},1)*kron( P_3{1}{3}'*Vz_star , P_3{1}{2}'*Vy_star ) ) );
[Vx_ddagger,~] = qr(K1,0);

% K2 step (y basis update)
K2 = sylvester( full(speye(Ny,Ny) - (dtn/2)*Dyy) ,...
    full( -(dtn/2)*(kron( (Dzz*Vz_star)'*Vz_star , speye(r1,r1) ) + kron( speye(r3,r3) , (Dxx*Vx_star)'*Vx_star )) ) ,...
    Vy_n*tens2mat(S_n,2)*kron( Vz_n'*Vz_star , Vx_n'*Vx_star )...
    + (-dtn/2)*(  (Vy_1*tens2mat(S_1,2)*kron( Vz_1'*Vz_star , (Dxx*Vx_1)'*Vx_star ))...
                        +((Dyy*Vy_1)*tens2mat(S_1,2)*kron( Vz_1'*Vz_star , Vx_1'*Vx_star ))...
                        +(Vy_1*tens2mat(S_1,2)*kron( (Dzz*Vz_1)'*Vz_star , Vx_1'*Vx_star ))...
                        +(P_1{1}{2}*tens2mat(P_1{2},2)*kron( P_1{1}{3}'*Vz_star , P_1{1}{1}'*Vx_star ))  )...
    + (dtn/2)*(  (Vy_2*tens2mat(S_2,2)*kron( Vz_2'*Vz_star , (Dxx*Vx_2)'*Vx_star ))...
                        +((Dyy*Vy_2)*tens2mat(S_2,2)*kron( Vz_2'*Vz_star , Vx_2'*Vx_star ))...
                        +(Vy_2*tens2mat(S_2,2)*kron( (Dzz*Vz_2)'*Vz_star , Vx_2'*Vx_star ))...
                        +(P_2{1}{2}*tens2mat(P_2{2},2)*kron( P_2{1}{3}'*Vz_star , P_2{1}{1}'*Vx_star ))  )...
    - (5/6)*dtn*(  ( E1_n{1}{2}*tens2mat(E1_n{2},2)*kron( E1_n{1}{3}'*Vz_star , (Dx*E1_n{1}{1})'*Vx_star ) )...
                  +( (Dy*E2_n{1}{2})*tens2mat(E2_n{2},2)*kron( E2_n{1}{3}'*Vz_star , E2_n{1}{1}'*Vx_star ) )...
                  +( E3_n{1}{2}*tens2mat(E3_n{2},2)*kron( (Dz*E3_n{1}{3})'*Vz_star , E3_n{1}{1}'*Vx_star ) )  )...
    - (-5/6)*dtn*(  ( E1_1{1}{2}*tens2mat(E1_1{2},2)*kron( E1_1{1}{3}'*Vz_star , (Dx*E1_1{1}{1})'*Vx_star ) )...
                      +( (Dy*E2_1{1}{2})*tens2mat(E2_1{2},2)*kron( E2_1{1}{3}'*Vz_star , E2_1{1}{1}'*Vx_star ) )...
                      +( E3_1{1}{2}*tens2mat(E3_1{2},2)*kron( (Dz*E3_1{1}{3})'*Vz_star , E3_1{1}{1}'*Vx_star ) )  )...
    - (1/2)*dtn*(  ( E1_2{1}{2}*tens2mat(E1_2{2},2)*kron( E1_2{1}{3}'*Vz_star , (Dx*E1_2{1}{1})'*Vx_star ) )...
                      +( (Dy*E2_2{1}{2})*tens2mat(E2_2{2},2)*kron( E2_2{1}{3}'*Vz_star , E2_2{1}{1}'*Vx_star ) )...
                      +( E3_2{1}{2}*tens2mat(E3_2{2},2)*kron( (Dz*E3_2{1}{3})'*Vz_star , E3_2{1}{1}'*Vx_star ) )  )...
    + (dtn/2)*( P_3{1}{2}*tens2mat(P_3{2},2)*kron( P_3{1}{3}'*Vz_star , P_3{1}{1}'*Vx_star ) ) );
[Vy_ddagger,~] = qr(K2,0);

% K3 step (z basis update)
K3 = sylvester( full(speye(Nz,Nz) - (dtn/2)*Dzz) ,...
    full( -(dtn/2)*(kron( (Dyy*Vy_star)'*Vy_star , speye(r1,r1) ) + kron( speye(r2,r2) , (Dxx*Vx_star)'*Vx_star )) ) ,...
    Vz_n*tens2mat(S_n,3)*kron( Vy_n'*Vy_star , Vx_n'*Vx_star )...
    + (-dtn/2)*(  (Vz_1*tens2mat(S_1,3)*kron( Vy_1'*Vy_star , (Dxx*Vx_1)'*Vx_star ))...
                        +(Vz_1*tens2mat(S_1,3)*kron( (Dyy*Vy_1)'*Vy_star , Vx_1'*Vx_star ))...
                        +((Dzz*Vz_1)*tens2mat(S_1,3)*kron( Vy_1'*Vy_star , Vx_1'*Vx_star ))...
                        +(P_1{1}{3}*tens2mat(P_1{2},3)*kron( P_1{1}{2}'*Vy_star , P_1{1}{1}'*Vx_star ))  )...
    + (dtn/2)*(  (Vz_2*tens2mat(S_2,3)*kron( Vy_2'*Vy_star , (Dxx*Vx_2)'*Vx_star ))...
                        +(Vz_2*tens2mat(S_2,3)*kron( (Dyy*Vy_2)'*Vy_star , Vx_2'*Vx_star ))...
                        +((Dzz*Vz_2)*tens2mat(S_2,3)*kron( Vy_2'*Vy_star , Vx_2'*Vx_star ))...
                        +(P_2{1}{3}*tens2mat(P_2{2},3)*kron( P_2{1}{2}'*Vy_star , P_2{1}{1}'*Vx_star ))  )...
    - (5/6)*dtn*(  ( E1_n{1}{3}*tens2mat(E1_n{2},3)*kron( E1_n{1}{2}'*Vy_star , (Dx*E1_n{1}{1})'*Vx_star ) )...
                  +( E2_n{1}{3}*tens2mat(E2_n{2},3)*kron( (Dy*E2_n{1}{2})'*Vy_star , E2_n{1}{1}'*Vx_star ) )...
                  +( (Dz*E3_n{1}{3})*tens2mat(E3_n{2},3)*kron( E3_n{1}{2}'*Vy_star , E3_n{1}{1}'*Vx_star ) )  )...
    - (-5/6)*dtn*(  ( E1_1{1}{3}*tens2mat(E1_1{2},3)*kron( E1_1{1}{2}'*Vy_star , (Dx*E1_1{1}{1})'*Vx_star ) )...
                      +( E2_1{1}{3}*tens2mat(E2_1{2},3)*kron( (Dy*E2_1{1}{2})'*Vy_star , E2_1{1}{1}'*Vx_star ) )...
                      +( (Dz*E3_1{1}{3})*tens2mat(E3_1{2},3)*kron( E3_1{1}{2}'*Vy_star , E3_1{1}{1}'*Vx_star ) )  )...
    - (1/2)*dtn*(  ( E1_2{1}{3}*tens2mat(E1_2{2},3)*kron( E1_2{1}{2}'*Vy_star , (Dx*E1_2{1}{1})'*Vx_star ) )...
                      +( E2_2{1}{3}*tens2mat(E2_2{2},3)*kron( (Dy*E2_2{1}{2})'*Vy_star , E2_2{1}{1}'*Vx_star ) )...
                      +( (Dz*E3_2{1}{3})*tens2mat(E3_2{2},3)*kron( E3_2{1}{2}'*Vy_star , E3_2{1}{1}'*Vx_star ) )  )...
    + (dtn/2)*( P_3{1}{3}*tens2mat(P_3{2},3)*kron( P_3{1}{2}'*Vy_star , P_3{1}{1}'*Vx_star ) ) );
[Vz_ddagger,~] = qr(K3,0);


% S step (update transfer tensor)
% Using Simoncini's solver requires square coefficient matrices. So, we
% must make the dimensions r1=r2=r3=R for the S step. These dimensions
% could be different after truncation. Note that Simoncini's solver will
% output a "square" order-3 tensor. Hence, RxRxR is okay if R>1. However,
% if R=1, then instead of 1x1x1, MATLAB will instead output a 1x1 array.
% So, we must force a 1x1x1 tensor in this situation after solving for
% S_nn using Simoncini's solver.
[Vx_3,Vy_3,Vz_3,R] = red_aug_S(Vx_ddagger,[Vx_2,Vx_1,Vx_n],Vy_ddagger,[Vy_2,Vy_1,Vy_n],Vz_ddagger,[Vz_2,Vz_1,Vz_n]);
% r1 = r2 = r3 = R after reduced augmentation. Could be different after
% truncation.

A1 = -(dtn/2)*Vy_3'*(Dyy*Vy_3);
A2 = speye(R,R)-(dtn/2)*Vz_3'*(Dzz*Vz_3);
A3 = -(dtn/2)*Vx_3'*(Dxx*Vx_3);
M1 = speye(R,R); %r3xr3
M2 = speye(R,R); %r2xr2
H1 = speye(R,R); %r1xr1
H3 = speye(R,R); %r3xr3
check_res = 1;

O1 = Vx_3'*Vx_n;
O2 = Vy_3'*Vy_n;
O3 = Vz_3'*Vz_n;
S = S_n;
O = lmlragen({O1,O2,O3},S);

O1 = Vx_3'*(Dxx*Vx_1);
O2 = Vy_3'*Vy_1;
O3 = Vz_3'*Vz_1;
S = (-dtn/2)*S_1;
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_3'*Vx_1;
O2 = Vy_3'*(Dyy*Vy_1);
O3 = Vz_3'*Vz_1;
S = (-dtn/2)*S_1;
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_3'*Vx_1;
O2 = Vy_3'*Vy_1;
O3 = Vz_3'*(Dzz*Vz_1);
S = (-dtn/2)*S_1;
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_3'*P_1{1}{1};
O2 = Vy_3'*P_1{1}{2};
O3 = Vz_3'*P_1{1}{3};
S = (-dtn/2)*P_1{2};
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_3'*(Dxx*Vx_2);
O2 = Vy_3'*Vy_2;
O3 = Vz_3'*Vz_2;
S = (dtn/2)*S_2;
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_3'*Vx_2;
O2 = Vy_3'*(Dyy*Vy_2);
O3 = Vz_3'*Vz_2;
S = (dtn/2)*S_2;
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_3'*Vx_2;
O2 = Vy_3'*Vy_2;
O3 = Vz_3'*(Dzz*Vz_2);
S = (dtn/2)*S_2;
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_3'*P_2{1}{1};
O2 = Vy_3'*P_2{1}{2};
O3 = Vz_3'*P_2{1}{3};
S = (dtn/2)*P_2{2};
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_3'*(Dx*E1_n{1}{1});
O2 = Vy_3'*E1_n{1}{2};
O3 = Vz_3'*E1_n{1}{3};
S = -((5/6)*dtn)*E1_n{2};
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_3'*E2_n{1}{1};
O2 = Vy_3'*(Dy*E2_n{1}{2});
O3 = Vz_3'*E2_n{1}{3};
S = -((5/6)*dtn)*E2_n{2};
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_3'*E3_n{1}{1};
O2 = Vy_3'*E3_n{1}{2};
O3 = Vz_3'*(Dz*E3_n{1}{3});
S = -((5/6)*dtn)*E3_n{2};
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_3'*(Dx*E1_1{1}{1});
O2 = Vy_3'*E1_1{1}{2};
O3 = Vz_3'*E1_1{1}{3};
S = -((-5/6)*dtn)*E1_1{2};
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_3'*E2_1{1}{1};
O2 = Vy_3'*(Dy*E2_1{1}{2});
O3 = Vz_3'*E2_1{1}{3};
S = -((-5/6)*dtn)*E2_1{2};
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_3'*E3_1{1}{1};
O2 = Vy_3'*E3_1{1}{2};
O3 = Vz_3'*(Dz*E3_1{1}{3});
S = -((-5/6)*dtn)*E3_1{2};
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_3'*(Dx*E1_2{1}{1});
O2 = Vy_3'*E1_2{1}{2};
O3 = Vz_3'*E1_2{1}{3};
S = -((1/2)*dtn)*E1_2{2};
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_3'*E2_2{1}{1};
O2 = Vy_3'*(Dy*E2_2{1}{2});
O3 = Vz_3'*E2_2{1}{3};
S = -((1/2)*dtn)*E2_2{2};
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_3'*E3_2{1}{1};
O2 = Vy_3'*E3_2{1}{2};
O3 = Vz_3'*(Dz*E3_2{1}{3});
S = -((1/2)*dtn)*E3_2{2};
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_3'*P_3{1}{1};
O2 = Vy_3'*P_3{1}{2};
O3 = Vz_3'*P_3{1}{3};
S = (dtn/2)*P_3{2};
O = O + lmlragen({O1,O2,O3},S);

[~,S_3] = simoncini_direct_solver(A1,A2,A3,M1,M2,H1,H3,O,check_res); %outputs double data type

% Truncation
% [U_3,G_3,~] = nonconstrun({Vx_3,Vy_3,Vz_3},S_3,tol);
% [U_3,G_3,~] = lomac_0({Vx_3,Vy_3,Vz_3},S_3,tol,w1,w2,w3,Nx,Ny,Nz,dx,dy,dz,rhoM);
% [U_3,G_3,~] = lomac_01({Vx_3,Vy_3,Vz_3},S_3,tol,w1,w2,w3,Nx,Ny,Nz,dx,dy,dz,rhoM,JxM,JyM,JzM,xvals,yvals,zvals);
[U_3,G_3,~] = lomac_012({Vx_3,Vy_3,Vz_3},S_3,tol,w1,w2,w3,Nx,Ny,Nz,dx,dy,dz,rhoM,JxM,JyM,JzM,kM,c,cc,xvals,yvals,zvals);

% Compute flux terms (at time t^(3) = t^n + dt/2)
A_3 = A_1; %A(t^(3))
B_3 = B_1; %B(t^(3))
C_3 = C_1; %C(t^(3))
E1_3 = {{TKR(A_3{1}{1},U_3{1}),TKR(A_3{1}{2},U_3{2}),TKR(A_3{1}{3},U_3{3})},tkron(A_3{2},G_3)}; %E1 = A*U
E2_3 = {{TKR(B_3{1}{1},U_3{1}),TKR(B_3{1}{2},U_3{2}),TKR(B_3{1}{3},U_3{3})},tkron(B_3{2},G_3)}; %E2 = B*U
E3_3 = {{TKR(C_3{1}{1},U_3{1}),TKR(C_3{1}{2},U_3{2}),TKR(C_3{1}{3},U_3{3})},tkron(C_3{2},G_3)}; %E3 = C*U
% Compute source term (at time t^(4) = t^{n+1})
P_4 = P(tn+dtn);% Pull out bases
Vx_3 = U_3{1};
Vy_3 = U_3{2};
Vz_3 = U_3{3};
S_3 = G_3;


%%
% % % STAGE 4, t^(4) = t^{n+1} % % %
% K Steps
[U_dagger,~,~] = IMEX111(U_n,G_n,MLR_n,A,B,C,P,tn,dtn,Dx,Dy,Dz,Dxx,Dyy,Dzz,Nx,Ny,Nz,tol,w1,w2,w3,dx,dy,dz,rhoM,JxM,JyM,JzM,kM,c,cc,xvals,yvals,zvals);
[Vx_star,Vy_star,Vz_star,r1,r2,r3] = red_aug_K(U_dagger{1},[Vx_3,Vx_2,Vx_1,Vx_n],U_dagger{2},[Vy_3,Vy_2,Vy_1,Vy_n],U_dagger{3},[Vz_3,Vz_2,Vz_1,Vz_n]);

% K1 step (x basis update)
K1 = sylvester( full(speye(Nx,Nx) - (dtn/2)*Dxx) ,...
    full( -(dtn/2)*(kron( (Dzz*Vz_star)'*Vz_star , speye(r2,r2) ) + kron( speye(r3,r3) , (Dyy*Vy_star)'*Vy_star )) ) ,...
    Vx_n*tens2mat(S_n,1)*kron( Vz_n'*Vz_star , Vy_n'*Vy_star )...
    + ((3/2)*dtn)*(  ((Dxx*Vx_1)*tens2mat(S_1,1)*kron( Vz_1'*Vz_star , Vy_1'*Vy_star ))...
                        +(Vx_1*tens2mat(S_1,1)*kron( Vz_1'*Vz_star , (Dyy*Vy_1)'*Vy_star ))...
                        +(Vx_1*tens2mat(S_1,1)*kron( (Dzz*Vz_1)'*Vz_star , Vy_1'*Vy_star ))...
                        +(P_1{1}{1}*tens2mat(P_1{2},1)*kron( P_1{1}{3}'*Vz_star , P_1{1}{2}'*Vy_star ))  )...
    + (-(3/2)*dtn)*(  ((Dxx*Vx_2)*tens2mat(S_2,1)*kron( Vz_2'*Vz_star , Vy_2'*Vy_star ))...
                        +(Vx_2*tens2mat(S_2,1)*kron( Vz_2'*Vz_star , (Dyy*Vy_2)'*Vy_star ))...
                        +(Vx_2*tens2mat(S_2,1)*kron( (Dzz*Vz_2)'*Vz_star , Vy_2'*Vy_star ))...
                        +(P_2{1}{1}*tens2mat(P_2{2},1)*kron( P_2{1}{3}'*Vz_star , P_2{1}{2}'*Vy_star ))  )...
    + (dtn/2)*(  ((Dxx*Vx_3)*tens2mat(S_3,1)*kron( Vz_3'*Vz_star , Vy_3'*Vy_star ))...
                        +(Vx_3*tens2mat(S_3,1)*kron( Vz_3'*Vz_star , (Dyy*Vy_3)'*Vy_star ))...
                        +(Vx_3*tens2mat(S_3,1)*kron( (Dzz*Vz_3)'*Vz_star , Vy_3'*Vy_star ))...
                        +(P_3{1}{1}*tens2mat(P_3{2},1)*kron( P_3{1}{3}'*Vz_star , P_3{1}{2}'*Vy_star ))  )...
    - (dtn/4)*(  ( (Dx*E1_n{1}{1})*tens2mat(E1_n{2},1)*kron( E1_n{1}{3}'*Vz_star , E1_n{1}{2}'*Vy_star ) )...
                  +( E2_n{1}{1}*tens2mat(E2_n{2},1)*kron( E2_n{1}{3}'*Vz_star , (Dy*E2_n{1}{2})'*Vy_star ) )...
                  +( E3_n{1}{1}*tens2mat(E3_n{2},1)*kron( (Dz*E3_n{1}{3})'*Vz_star , E3_n{1}{2}'*Vy_star ) )  )...
    - (7*dtn/4)*(  ( (Dx*E1_1{1}{1})*tens2mat(E1_1{2},1)*kron( E1_1{1}{3}'*Vz_star , E1_1{1}{2}'*Vy_star ) )...
                      +( E2_1{1}{1}*tens2mat(E2_1{2},1)*kron( E2_1{1}{3}'*Vz_star , (Dy*E2_1{1}{2})'*Vy_star ) )...
                      +( E3_1{1}{1}*tens2mat(E3_1{2},1)*kron( (Dz*E3_1{1}{3})'*Vz_star , E3_1{1}{2}'*Vy_star ) )  )...
    - ((3/4)*dtn)*(  ( (Dx*E1_2{1}{1})*tens2mat(E1_2{2},1)*kron( E1_2{1}{3}'*Vz_star , E1_2{1}{2}'*Vy_star ) )...
                      +( E2_2{1}{1}*tens2mat(E2_2{2},1)*kron( E2_2{1}{3}'*Vz_star , (Dy*E2_2{1}{2})'*Vy_star ) )...
                      +( E3_2{1}{1}*tens2mat(E3_2{2},1)*kron( (Dz*E3_2{1}{3})'*Vz_star , E3_2{1}{2}'*Vy_star ) )  )...
    - (-(7/4)*dtn)*(  ( (Dx*E1_3{1}{1})*tens2mat(E1_3{2},1)*kron( E1_3{1}{3}'*Vz_star , E1_3{1}{2}'*Vy_star ) )...
                      +( E2_3{1}{1}*tens2mat(E2_3{2},1)*kron( E2_3{1}{3}'*Vz_star , (Dy*E2_3{1}{2})'*Vy_star ) )...
                      +( E3_3{1}{1}*tens2mat(E3_3{2},1)*kron( (Dz*E3_3{1}{3})'*Vz_star , E3_3{1}{2}'*Vy_star ) )  )...
    + (dtn/2)*( P_4{1}{1}*tens2mat(P_4{2},1)*kron( P_4{1}{3}'*Vz_star , P_4{1}{2}'*Vy_star ) ) );
[Vx_ddagger,~] = qr(K1,0);

% K2 step (y basis update)
K2 = sylvester( full(speye(Ny,Ny) - (dtn/2)*Dyy) ,...
    full( -(dtn/2)*(kron( (Dzz*Vz_star)'*Vz_star , speye(r1,r1) ) + kron( speye(r3,r3) , (Dxx*Vx_star)'*Vx_star )) ) ,...
    Vy_n*tens2mat(S_n,2)*kron( Vz_n'*Vz_star , Vx_n'*Vx_star )...
    + ((3/2)*dtn)*(  (Vy_1*tens2mat(S_1,2)*kron( Vz_1'*Vz_star , (Dxx*Vx_1)'*Vx_star ))...
                        +((Dyy*Vy_1)*tens2mat(S_1,2)*kron( Vz_1'*Vz_star , Vx_1'*Vx_star ))...
                        +(Vy_1*tens2mat(S_1,2)*kron( (Dzz*Vz_1)'*Vz_star , Vx_1'*Vx_star ))...
                        +(P_1{1}{2}*tens2mat(P_1{2},2)*kron( P_1{1}{3}'*Vz_star , P_1{1}{1}'*Vx_star ))  )...
    + (-(3/2)*dtn)*(  (Vy_2*tens2mat(S_2,2)*kron( Vz_2'*Vz_star , (Dxx*Vx_2)'*Vx_star ))...
                        +((Dyy*Vy_2)*tens2mat(S_2,2)*kron( Vz_2'*Vz_star , Vx_2'*Vx_star ))...
                        +(Vy_2*tens2mat(S_2,2)*kron( (Dzz*Vz_2)'*Vz_star , Vx_2'*Vx_star ))...
                        +(P_2{1}{2}*tens2mat(P_2{2},2)*kron( P_2{1}{3}'*Vz_star , P_2{1}{1}'*Vx_star ))  )...
    + (dtn/2)*(  (Vy_3*tens2mat(S_3,2)*kron( Vz_3'*Vz_star , (Dxx*Vx_3)'*Vx_star ))...
                        +((Dyy*Vy_3)*tens2mat(S_3,2)*kron( Vz_3'*Vz_star , Vx_3'*Vx_star ))...
                        +(Vy_3*tens2mat(S_3,2)*kron( (Dzz*Vz_3)'*Vz_star , Vx_3'*Vx_star ))...
                        +(P_3{1}{2}*tens2mat(P_3{2},2)*kron( P_3{1}{3}'*Vz_star , P_3{1}{1}'*Vx_star ))  )...
    - (dtn/4)*(  ( E1_n{1}{2}*tens2mat(E1_n{2},2)*kron( E1_n{1}{3}'*Vz_star , (Dx*E1_n{1}{1})'*Vx_star ) )...
                  +( (Dy*E2_n{1}{2})*tens2mat(E2_n{2},2)*kron( E2_n{1}{3}'*Vz_star , E2_n{1}{1}'*Vx_star ) )...
                  +( E3_n{1}{2}*tens2mat(E3_n{2},2)*kron( (Dz*E3_n{1}{3})'*Vz_star , E3_n{1}{1}'*Vx_star ) )  )...
    - (7*dtn/4)*(  ( E1_1{1}{2}*tens2mat(E1_1{2},2)*kron( E1_1{1}{3}'*Vz_star , (Dx*E1_1{1}{1})'*Vx_star ) )...
                      +( (Dy*E2_1{1}{2})*tens2mat(E2_1{2},2)*kron( E2_1{1}{3}'*Vz_star , E2_1{1}{1}'*Vx_star ) )...
                      +( E3_1{1}{2}*tens2mat(E3_1{2},2)*kron( (Dz*E3_1{1}{3})'*Vz_star , E3_1{1}{1}'*Vx_star ) )  )...
    - (3*dtn/4)*(  ( E1_2{1}{2}*tens2mat(E1_2{2},2)*kron( E1_2{1}{3}'*Vz_star , (Dx*E1_2{1}{1})'*Vx_star ) )...
                      +( (Dy*E2_2{1}{2})*tens2mat(E2_2{2},2)*kron( E2_2{1}{3}'*Vz_star , E2_2{1}{1}'*Vx_star ) )...
                      +( E3_2{1}{2}*tens2mat(E3_2{2},2)*kron( (Dz*E3_2{1}{3})'*Vz_star , E3_2{1}{1}'*Vx_star ) )  )...
    - (-7*dtn/4)*(  ( E1_3{1}{2}*tens2mat(E1_3{2},2)*kron( E1_3{1}{3}'*Vz_star , (Dx*E1_3{1}{1})'*Vx_star ) )...
                      +( (Dy*E2_3{1}{2})*tens2mat(E2_3{2},2)*kron( E2_3{1}{3}'*Vz_star , E2_3{1}{1}'*Vx_star ) )...
                      +( E3_3{1}{2}*tens2mat(E3_3{2},2)*kron( (Dz*E3_3{1}{3})'*Vz_star , E3_3{1}{1}'*Vx_star ) )  )...
    + (dtn/2)*( P_4{1}{2}*tens2mat(P_4{2},2)*kron( P_4{1}{3}'*Vz_star , P_4{1}{1}'*Vx_star ) ) );
[Vy_ddagger,~] = qr(K2,0);

% K3 step (z basis update)
K3 = sylvester( full(speye(Nz,Nz) - (dtn/2)*Dzz) ,...
    full( -(dtn/2)*(kron( (Dyy*Vy_star)'*Vy_star , speye(r1,r1) ) + kron( speye(r2,r2) , (Dxx*Vx_star)'*Vx_star )) ) ,...
    Vz_n*tens2mat(S_n,3)*kron( Vy_n'*Vy_star , Vx_n'*Vx_star )...
    + ((3/2)*dtn)*(  (Vz_1*tens2mat(S_1,3)*kron( Vy_1'*Vy_star , (Dxx*Vx_1)'*Vx_star ))...
                        +(Vz_1*tens2mat(S_1,3)*kron( (Dyy*Vy_1)'*Vy_star , Vx_1'*Vx_star ))...
                        +((Dzz*Vz_1)*tens2mat(S_1,3)*kron( Vy_1'*Vy_star , Vx_1'*Vx_star ))...
                        +(P_1{1}{3}*tens2mat(P_1{2},3)*kron( P_1{1}{2}'*Vy_star , P_1{1}{1}'*Vx_star ))  )...
    + (-(3/2)*dtn)*(  (Vz_2*tens2mat(S_2,3)*kron( Vy_2'*Vy_star , (Dxx*Vx_2)'*Vx_star ))...
                        +(Vz_2*tens2mat(S_2,3)*kron( (Dyy*Vy_2)'*Vy_star , Vx_2'*Vx_star ))...
                        +((Dzz*Vz_2)*tens2mat(S_2,3)*kron( Vy_2'*Vy_star , Vx_2'*Vx_star ))...
                        +(P_2{1}{3}*tens2mat(P_2{2},3)*kron( P_2{1}{2}'*Vy_star , P_2{1}{1}'*Vx_star ))  )...
    + (dtn/2)*(  (Vz_3*tens2mat(S_3,3)*kron( Vy_3'*Vy_star , (Dxx*Vx_3)'*Vx_star ))...
                        +(Vz_3*tens2mat(S_3,3)*kron( (Dyy*Vy_3)'*Vy_star , Vx_3'*Vx_star ))...
                        +((Dzz*Vz_3)*tens2mat(S_3,3)*kron( Vy_3'*Vy_star , Vx_3'*Vx_star ))...
                        +(P_3{1}{3}*tens2mat(P_3{2},3)*kron( P_3{1}{2}'*Vy_star , P_3{1}{1}'*Vx_star ))  )...
    - (dtn/4)*(  ( E1_n{1}{3}*tens2mat(E1_n{2},3)*kron( E1_n{1}{2}'*Vy_star , (Dx*E1_n{1}{1})'*Vx_star ) )...
                  +( E2_n{1}{3}*tens2mat(E2_n{2},3)*kron( (Dy*E2_n{1}{2})'*Vy_star , E2_n{1}{1}'*Vx_star ) )...
                  +( (Dz*E3_n{1}{3})*tens2mat(E3_n{2},3)*kron( E3_n{1}{2}'*Vy_star , E3_n{1}{1}'*Vx_star ) )  )...
    - (7*dtn/4)*(  ( E1_1{1}{3}*tens2mat(E1_1{2},3)*kron( E1_1{1}{2}'*Vy_star , (Dx*E1_1{1}{1})'*Vx_star ) )...
                      +( E2_1{1}{3}*tens2mat(E2_1{2},3)*kron( (Dy*E2_1{1}{2})'*Vy_star , E2_1{1}{1}'*Vx_star ) )...
                      +( (Dz*E3_1{1}{3})*tens2mat(E3_1{2},3)*kron( E3_1{1}{2}'*Vy_star , E3_1{1}{1}'*Vx_star ) )  )...
    - (3*dtn/4)*(  ( E1_2{1}{3}*tens2mat(E1_2{2},3)*kron( E1_2{1}{2}'*Vy_star , (Dx*E1_2{1}{1})'*Vx_star ) )...
                      +( E2_2{1}{3}*tens2mat(E2_2{2},3)*kron( (Dy*E2_2{1}{2})'*Vy_star , E2_2{1}{1}'*Vx_star ) )...
                      +( (Dz*E3_2{1}{3})*tens2mat(E3_2{2},3)*kron( E3_2{1}{2}'*Vy_star , E3_2{1}{1}'*Vx_star ) )  )...
    - (-7*dtn/4)*(  ( E1_3{1}{3}*tens2mat(E1_3{2},3)*kron( E1_3{1}{2}'*Vy_star , (Dx*E1_3{1}{1})'*Vx_star ) )...
                      +( E2_3{1}{3}*tens2mat(E2_3{2},3)*kron( (Dy*E2_3{1}{2})'*Vy_star , E2_3{1}{1}'*Vx_star ) )...
                      +( (Dz*E3_3{1}{3})*tens2mat(E3_3{2},3)*kron( E3_3{1}{2}'*Vy_star , E3_3{1}{1}'*Vx_star ) )  )...
    + (dtn/2)*( P_4{1}{3}*tens2mat(P_4{2},3)*kron( P_4{1}{2}'*Vy_star , P_4{1}{1}'*Vx_star ) ) );
[Vz_ddagger,~] = qr(K3,0);


% S step (update transfer tensor)
% Using Simoncini's solver requires square coefficient matrices. So, we
% must make the dimensions r1=r2=r3=R for the S step. These dimensions
% could be different after truncation. Note that Simoncini's solver will
% output a "square" order-3 tensor. Hence, RxRxR is okay if R>1. However,
% if R=1, then instead of 1x1x1, MATLAB will instead output a 1x1 array.
% So, we must force a 1x1x1 tensor in this situation after solving for
% S_nn using Simoncini's solver.
[Vx_nn,Vy_nn,Vz_nn,R] = red_aug_S(Vx_ddagger,[Vx_3,Vx_2,Vx_1,Vx_n],Vy_ddagger,[Vy_3,Vy_2,Vy_1,Vy_n],Vz_ddagger,[Vz_3,Vz_2,Vz_1,Vz_n]);
% r1 = r2 = r3 = R after reduced augmentation. Could be different after
% truncation.

A1 = -(dtn/2)*Vy_nn'*(Dyy*Vy_nn);
A2 = speye(R,R)-(dtn/2)*Vz_nn'*(Dzz*Vz_nn);
A3 = -(dtn/2)*Vx_nn'*(Dxx*Vx_nn);
M1 = speye(R,R); %r3xr3
M2 = speye(R,R); %r2xr2
H1 = speye(R,R); %r1xr1
H3 = speye(R,R); %r3xr3
check_res = 1;

O1 = Vx_nn'*Vx_n;
O2 = Vy_nn'*Vy_n;
O3 = Vz_nn'*Vz_n;
S = S_n;
O = lmlragen({O1,O2,O3},S);

O1 = Vx_nn'*(Dxx*Vx_1);
O2 = Vy_nn'*Vy_1;
O3 = Vz_nn'*Vz_1;
S = (3*dtn/2)*S_1;
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_nn'*Vx_1;
O2 = Vy_nn'*(Dyy*Vy_1);
O3 = Vz_nn'*Vz_1;
S = (3*dtn/2)*S_1;
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_nn'*Vx_1;
O2 = Vy_nn'*Vy_1;
O3 = Vz_nn'*(Dzz*Vz_1);
S = (3*dtn/2)*S_1;
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_nn'*P_1{1}{1};
O2 = Vy_nn'*P_1{1}{2};
O3 = Vz_nn'*P_1{1}{3};
S = (3*dtn/2)*P_1{2};
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_nn'*(Dxx*Vx_2);
O2 = Vy_nn'*Vy_2;
O3 = Vz_nn'*Vz_2;
S = (-3*dtn/2)*S_2;
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_nn'*Vx_2;
O2 = Vy_nn'*(Dyy*Vy_2);
O3 = Vz_nn'*Vz_2;
S = (-3*dtn/2)*S_2;
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_nn'*Vx_2;
O2 = Vy_nn'*Vy_2;
O3 = Vz_nn'*(Dzz*Vz_2);
S = (-3*dtn/2)*S_2;
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_nn'*P_2{1}{1};
O2 = Vy_nn'*P_2{1}{2};
O3 = Vz_nn'*P_2{1}{3};
S = (-3*dtn/2)*P_2{2};
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_nn'*(Dxx*Vx_3);
O2 = Vy_nn'*Vy_3;
O3 = Vz_nn'*Vz_3;
S = (dtn/2)*S_3;
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_nn'*Vx_3;
O2 = Vy_nn'*(Dyy*Vy_3);
O3 = Vz_nn'*Vz_3;
S = (dtn/2)*S_3;
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_nn'*Vx_3;
O2 = Vy_nn'*Vy_3;
O3 = Vz_nn'*(Dzz*Vz_3);
S = (dtn/2)*S_3;
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_nn'*P_3{1}{1};
O2 = Vy_nn'*P_3{1}{2};
O3 = Vz_nn'*P_3{1}{3};
S = (dtn/2)*P_3{2};
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_nn'*(Dx*E1_n{1}{1});
O2 = Vy_nn'*E1_n{1}{2};
O3 = Vz_nn'*E1_n{1}{3};
S = -((1/4)*dtn)*E1_n{2};
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_nn'*E2_n{1}{1};
O2 = Vy_nn'*(Dy*E2_n{1}{2});
O3 = Vz_nn'*E2_n{1}{3};
S = -((1/4)*dtn)*E2_n{2};
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_nn'*E3_n{1}{1};
O2 = Vy_nn'*E3_n{1}{2};
O3 = Vz_nn'*(Dz*E3_n{1}{3});
S = -((1/4)*dtn)*E3_n{2};
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_nn'*(Dx*E1_1{1}{1});
O2 = Vy_nn'*E1_1{1}{2};
O3 = Vz_nn'*E1_1{1}{3};
S = -((7/4)*dtn)*E1_1{2};
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_nn'*E2_1{1}{1};
O2 = Vy_nn'*(Dy*E2_1{1}{2});
O3 = Vz_nn'*E2_1{1}{3};
S = -((7/4)*dtn)*E2_1{2};
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_nn'*E3_1{1}{1};
O2 = Vy_nn'*E3_1{1}{2};
O3 = Vz_nn'*(Dz*E3_1{1}{3});
S = -((7/4)*dtn)*E3_1{2};
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_nn'*(Dx*E1_2{1}{1});
O2 = Vy_nn'*E1_2{1}{2};
O3 = Vz_nn'*E1_2{1}{3};
S = -((3/4)*dtn)*E1_2{2};
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_nn'*E2_2{1}{1};
O2 = Vy_nn'*(Dy*E2_2{1}{2});
O3 = Vz_nn'*E2_2{1}{3};
S = -((3/4)*dtn)*E2_2{2};
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_nn'*E3_2{1}{1};
O2 = Vy_nn'*E3_2{1}{2};
O3 = Vz_nn'*(Dz*E3_2{1}{3});
S = -((3/4)*dtn)*E3_2{2};
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_nn'*(Dx*E1_3{1}{1});
O2 = Vy_nn'*E1_3{1}{2};
O3 = Vz_nn'*E1_3{1}{3};
S = -(-(7/4)*dtn)*E1_3{2};
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_nn'*E2_3{1}{1};
O2 = Vy_nn'*(Dy*E2_3{1}{2});
O3 = Vz_nn'*E2_3{1}{3};
S = -(-(7/4)*dtn)*E2_3{2};
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_nn'*E3_3{1}{1};
O2 = Vy_nn'*E3_3{1}{2};
O3 = Vz_nn'*(Dz*E3_3{1}{3});
S = -(-(7/4)*dtn)*E3_3{2};
O = O + lmlragen({O1,O2,O3},S);

O1 = Vx_nn'*P_4{1}{1};
O2 = Vy_nn'*P_4{1}{2};
O3 = Vz_nn'*P_4{1}{3};
S = (dtn/2)*P_4{2};
O = O + lmlragen({O1,O2,O3},S);

[~,S_nn] = simoncini_direct_solver(A1,A2,A3,M1,M2,H1,H3,O,check_res); %outputs double data type

% Truncation
% [U_nn,G_nn,MLR_nn] = nonconstrun({Vx_nn,Vy_nn,Vz_nn},S_nn,tol);
% [U_nn,G_nn,MLR_nn] = lomac_0({Vx_nn,Vy_nn,Vz_nn},S_nn,tol,w1,w2,w3,Nx,Ny,Nz,dx,dy,dz,rhoM);
% [U_nn,G_nn,MLR_nn] = lomac_01({Vx_nn,Vy_nn,Vz_nn},S_nn,tol,w1,w2,w3,Nx,Ny,Nz,dx,dy,dz,rhoM,JxM,JyM,JzM,xvals,yvals,zvals);
[U_nn,G_nn,MLR_nn] = lomac_012({Vx_nn,Vy_nn,Vz_nn},S_nn,tol,w1,w2,w3,Nx,Ny,Nz,dx,dy,dz,rhoM,JxM,JyM,JzM,kM,c,cc,xvals,yvals,zvals);
end











