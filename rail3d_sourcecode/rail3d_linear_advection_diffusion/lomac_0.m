function [Uout,Gout,rank] = lomac_0(Uin,Gin,tol,w1,w2,w3,Nx,Ny,Nz,dx,dy,dz,rhoM)
% LoMaC truncation. Conserves zeroth moment (i.e., mass).
rhoH = dx*dy*dz*sum(Uin{1},1)*tens2mat(Gin,1)*kron(sum(Uin{3},1)',sum(Uin{2},1)');

% Compute f1
Vx_f1 = w1.*[ones(Nx,1)];
Vy_f1 = w2.*[ones(Ny,1)];
Vz_f1 = w3.*[ones(Nz,1)];
S_f1 = [rhoH/(dx*dy*dz*sum(w1)*sum(w2)*sum(w3))];
% Compute and truncate f2=f-f1
[U2,G2] = tensoraddition(Uin,Gin,{Vx_f1,Vy_f1,Vz_f1},-S_f1);
[G2_U,G2_G,~] = mlsvd(G2,tol); %truncate wrt tol
SU_f2 = cell(1,3);
SU_f2{1} = U2{1}*G2_U{1};
SU_f2{2} = U2{2}*G2_U{2};
SU_f2{3} = U2{3}*G2_U{3};
SG_f2 = G2_G;
if numel(SG_f2)==1
    SU_f2{3} = [1];
end
if numel(SU_f2)==2 %if resulting hosvd/mlsvd is, e.g., size 1x1x1 (output as matrix)
    SU_f2 = {SU_f2{1},SU_f2{2},[1]};
end

% Compute PN(trun(f2)), to subtract later to ensure zero moments.
rho2 = dx*dy*dz*sum(SU_f2{1},1)*tens2mat(SG_f2,1)*kron(sum(SU_f2{3},1)',sum(SU_f2{2},1)');
S_f2P = [rho2/(dx*dy*dz*sum(w1)*sum(w2)*sum(w3))];
S_fM = [rhoM/(dx*dy*dz*sum(w1)*sum(w2)*sum(w3))];
% Since they share the same O.N. basis, only add(subtract) core matrices.
% Redefine fM = fM-PN(trun(f2))

% Final solution (make 1d bases orthonormal)
[U,G] = tensoraddition({Vx_f1,Vy_f1,Vz_f1},S_fM-S_f2P,SU_f2,SG_f2);
[Vx,Rx] = qr(U{1},0);
[Vy,Ry] = qr(U{2},0);
[Vz,Rz] = qr(U{3},0);
[U,G] = mlsvd(lmlragen({Rx,Ry,Rz},G));
U{1} = Vx*U{1};
U{2} = Vy*U{2};
U{3} = Vz*U{3};
if numel(U)==2 %if resulting hosvd/mlsvd is, e.g., size 1x1x1 (output as matrix)
    U = {U{1},U{2},[1]};
end
Uout = U;
Gout = G;
sz = size(Gout);
if numel(sz)==2
    rank = [sz(1),sz(2),1];
else
    rank = [sz(1),sz(2),sz(3)];
end
end