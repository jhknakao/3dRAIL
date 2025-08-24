function [Vx_star,Vy_star,Vz_star,rx,ry,rz] = red_aug_K(Vx_dagger,Vx_RK,Vy_dagger,Vy_RK,Vz_dagger,Vz_RK)
% Reduced augmentation
% V_dagger is backward Euler prediction
% Vx_RK, Vy_RK, Vz_RK hold the bases from the previous RK stages

[Qx,Rx] = qr([Vx_dagger,Vx_RK],0);
[Qy,Ry] = qr([Vy_dagger,Vy_RK],0);
[Qz,Rz] = qr([Vz_dagger,Vz_RK],0);
[Vx_temp,Sx_temp,~] = svd(Rx,0);
[Vy_temp,Sy_temp,~] = svd(Ry,0);
[Vz_temp,Sz_temp,~] = svd(Rz,0);
rx = find(diag(Sx_temp)>1.0e-12,1,'last');
ry = find(diag(Sy_temp)>1.0e-12,1,'last');
rz = find(diag(Sz_temp)>1.0e-12,1,'last');

Vx_star = Qx*Vx_temp(:,1:rx);
Vy_star = Qy*Vy_temp(:,1:ry);
Vz_star = Qz*Vz_temp(:,1:rz);
end