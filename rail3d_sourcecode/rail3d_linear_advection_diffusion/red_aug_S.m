function [Vx_nn,Vy_nn,Vz_nn,R] = red_aug_S(Vx_ddagger,Vx_RK,Vy_ddagger,Vy_RK,Vz_ddagger,Vz_RK)
% Reduced augmentation
% V_ddagger from K steps
% Vx_RK, Vy_RK, Vz_RK hold the bases from the previous RK stages

[Qx,Rx] = qr([Vx_ddagger,Vx_RK],0);
[Qy,Ry] = qr([Vy_ddagger,Vy_RK],0);
[Qz,Rz] = qr([Vz_ddagger,Vz_RK],0);
[Vx_temp,Sx_temp,~] = svd(Rx,0);
[Vy_temp,Sy_temp,~] = svd(Ry,0);
[Vz_temp,Sz_temp,~] = svd(Rz,0);
rx = find(diag(Sx_temp)>1.0e-12,1,'last');
ry = find(diag(Sy_temp)>1.0e-12,1,'last');
rz = find(diag(Sz_temp)>1.0e-12,1,'last');

R = min(max([rx,ry,rz]),min([length(Rx),length(Ry),length(Rz)]));
%R = max([rx,ry,rz]);
Vx_nn = Qx*Vx_temp(:,1:R);
Vy_nn = Qy*Vy_temp(:,1:R);
Vz_nn = Qz*Vz_temp(:,1:R);
end