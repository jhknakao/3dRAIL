function finf = QCM(xvals,yvals,zvals,rhoM,JxM,JyM,JzM,kM)
% Quadrature Corrected Moment Matching
% "An equilibrium-preserving discretization for the nonlinear
% Rosenbluth-Fokker-Planck operator in arbitrary multi-dimensional
% geometry" by WT Taitano, L Chacon, AN Simakov, JSC 2024.

dx = xvals(2)-xvals(1);
dy = yvals(2)-yvals(1);
dz = zvals(2)-zvals(1);

[x,y,z] = meshgrid(xvals,yvals,zvals);
% Meshgrid ordering is y-x-z, so reorder to x-y-z
x = permute(x,[2,1,3]);
y = permute(y,[2,1,3]);
z = permute(z,[2,1,3]); %not technically necessary since third component is not affected
tol2 = 1.0e-13;
tol3 = 1.0e-15;
R = 1/6;
Rk_norm = 1;
Mk = [pi^(3/2);0;0;0;3]; %theoretically exact moments
n_k = Mk(1); ux_k = Mk(2); uy_k = Mk(3); uz_k = Mk(4); T_k = Mk(5);
Rk = zeros(5,1);
Rk(1,1) = rhoM - dx*dy*dz*sum(sum(sum( (n_k/(2*pi*R*T_k)^1.5)*exp(-((x-ux_k).^2+(y-uy_k).^2+(z-uz_k).^2)/(2*R*T_k)) )));
Rk(2,1) = JxM - dx*dy*dz*sum(sum(sum( (n_k/(2*pi*R*T_k)^1.5)*exp(-((x-ux_k).^2+(y-uy_k).^2+(z-uz_k).^2)/(2*R*T_k)).*x )));
Rk(3,1) = JyM - dx*dy*dz*sum(sum(sum( (n_k/(2*pi*R*T_k)^1.5)*exp(-((x-ux_k).^2+(y-uy_k).^2+(z-uz_k).^2)/(2*R*T_k)).*y )));
Rk(4,1) = JzM - dx*dy*dz*sum(sum(sum( (n_k/(2*pi*R*T_k)^1.5)*exp(-((x-ux_k).^2+(y-uy_k).^2+(z-uz_k).^2)/(2*R*T_k)).*z )));
Rk(5,1) = kM - dx*dy*dz*sum(sum(sum( (n_k/(2*pi*R*T_k)^1.5)*exp(-((x-ux_k).^2+(y-uy_k).^2+(z-uz_k).^2)/(2*R*T_k)).*(x.^2+y.^2+z.^2)/2 )));
R0_norm = norm(Rk,2);
count = 0;
while Rk_norm > tol2*R0_norm + tol3
Jk = zeros(5,5);
Jk(1,1) = - dx*dy*dz*sum(sum(sum( (1/(2*pi*R*T_k)^1.5)*exp(-((x-ux_k).^2+(y-uy_k).^2+(z-uz_k).^2)/(2*R*T_k)) )));
Jk(1,2) = - dx*dy*dz*sum(sum(sum( (n_k/(2*pi*R*T_k)^1.5)*((x-ux_k)/(R*T_k)).*exp(-((x-ux_k).^2+(y-uy_k).^2+(z-uz_k).^2)/(2*R*T_k)) )));
Jk(1,3) = - dx*dy*dz*sum(sum(sum( (n_k/(2*pi*R*T_k)^1.5)*((y-uy_k)/(R*T_k)).*exp(-((x-ux_k).^2+(y-uy_k).^2+(z-uz_k).^2)/(2*R*T_k)) )));
Jk(1,4) = - dx*dy*dz*sum(sum(sum( (n_k/(2*pi*R*T_k)^1.5)*((z-uz_k)/(R*T_k)).*exp(-((x-ux_k).^2+(y-uy_k).^2+(z-uz_k).^2)/(2*R*T_k)) )));
Jk(1,5) = - dx*dy*dz*sum(sum(sum( (n_k/(2*pi*R)^1.5)*exp(-((x-ux_k).^2+(y-uy_k).^2+(z-uz_k).^2)/(2*R*T_k)).*(-3/(2*T_k^2.5) + ((x-ux_k).^2+(y-uy_k).^2+(z-uz_k).^2)/(2*R*T_k^3.5)) )));

Jk(2,1) = - dx*dy*dz*sum(sum(sum( (1/(2*pi*R*T_k)^1.5)*exp(-((x-ux_k).^2+(y-uy_k).^2+(z-uz_k).^2)/(2*R*T_k)).*x )));
Jk(2,2) = - dx*dy*dz*sum(sum(sum( (n_k/(2*pi*R*T_k)^1.5)*((x-ux_k)/(R*T_k)).*exp(-((x-ux_k).^2+(y-uy_k).^2+(z-uz_k).^2)/(2*R*T_k)).*x )));
Jk(2,3) = - dx*dy*dz*sum(sum(sum( (n_k/(2*pi*R*T_k)^1.5)*((y-uy_k)/(R*T_k)).*exp(-((x-ux_k).^2+(y-uy_k).^2+(z-uz_k).^2)/(2*R*T_k)).*x )));
Jk(2,4) = - dx*dy*dz*sum(sum(sum( (n_k/(2*pi*R*T_k)^1.5)*((z-uz_k)/(R*T_k)).*exp(-((x-ux_k).^2+(y-uy_k).^2+(z-uz_k).^2)/(2*R*T_k)).*x )));
Jk(2,5) = - dx*dy*dz*sum(sum(sum( (n_k/(2*pi*R)^1.5)*exp(-((x-ux_k).^2+(y-uy_k).^2+(z-uz_k).^2)/(2*R*T_k)).*(-3/(2*T_k^2.5) + ((x-ux_k).^2+(y-uy_k).^2+(z-uz_k).^2)/(2*R*T_k^3.5)).*x )));

Jk(3,1) = - dx*dy*dz*sum(sum(sum( (1/(2*pi*R*T_k)^1.5)*exp(-((x-ux_k).^2+(y-uy_k).^2+(z-uz_k).^2)/(2*R*T_k)).*y )));
Jk(3,2) = - dx*dy*dz*sum(sum(sum( (n_k/(2*pi*R*T_k)^1.5)*((x-ux_k)/(R*T_k)).*exp(-((x-ux_k).^2+(y-uy_k).^2+(z-uz_k).^2)/(2*R*T_k)).*y )));
Jk(3,3) = - dx*dy*dz*sum(sum(sum( (n_k/(2*pi*R*T_k)^1.5)*((y-uy_k)/(R*T_k)).*exp(-((x-ux_k).^2+(y-uy_k).^2+(z-uz_k).^2)/(2*R*T_k)).*y )));
Jk(3,4) = - dx*dy*dz*sum(sum(sum( (n_k/(2*pi*R*T_k)^1.5)*((z-uz_k)/(R*T_k)).*exp(-((x-ux_k).^2+(y-uy_k).^2+(z-uz_k).^2)/(2*R*T_k)).*y )));
Jk(3,5) = - dx*dy*dz*sum(sum(sum( (n_k/(2*pi*R)^1.5)*exp(-((x-ux_k).^2+(y-uy_k).^2+(z-uz_k).^2)/(2*R*T_k)).*(-3/(2*T_k^2.5) + ((x-ux_k).^2+(y-uy_k).^2+(z-uz_k).^2)/(2*R*T_k^3.5)).*y )));

Jk(4,1) = - dx*dy*dz*sum(sum(sum( (1/(2*pi*R*T_k)^1.5)*exp(-((x-ux_k).^2+(y-uy_k).^2+(z-uz_k).^2)/(2*R*T_k)).*z )));
Jk(4,2) = - dx*dy*dz*sum(sum(sum( (n_k/(2*pi*R*T_k)^1.5)*((x-ux_k)/(R*T_k)).*exp(-((x-ux_k).^2+(y-uy_k).^2+(z-uz_k).^2)/(2*R*T_k)).*z )));
Jk(4,3) = - dx*dy*dz*sum(sum(sum( (n_k/(2*pi*R*T_k)^1.5)*((y-uy_k)/(R*T_k)).*exp(-((x-ux_k).^2+(y-uy_k).^2+(z-uz_k).^2)/(2*R*T_k)).*z )));
Jk(4,4) = - dx*dy*dz*sum(sum(sum( (n_k/(2*pi*R*T_k)^1.5)*((z-uz_k)/(R*T_k)).*exp(-((x-ux_k).^2+(y-uy_k).^2+(z-uz_k).^2)/(2*R*T_k)).*z )));
Jk(4,5) = - dx*dy*dz*sum(sum(sum( (n_k/(2*pi*R)^1.5)*exp(-((x-ux_k).^2+(y-uy_k).^2+(z-uz_k).^2)/(2*R*T_k)).*(-3/(2*T_k^2.5) + ((x-ux_k).^2+(y-uy_k).^2+(z-uz_k).^2)/(2*R*T_k^3.5)).*z )));

Jk(5,1) = - dx*dy*dz*sum(sum(sum( (1/(2*pi*R*T_k)^1.5)*exp(-((x-ux_k).^2+(y-uy_k).^2+(z-uz_k).^2)/(2*R*T_k)).*(((x-ux_k).^2+(y-uy_k).^2+(z-uz_k).^2)/2) )));
Jk(5,2) = - dx*dy*dz*sum(sum(sum( (n_k/(2*pi*R*T_k)^1.5)*((x-ux_k)/(R*T_k)).*exp(-((x-ux_k).^2+(y-uy_k).^2+(z-uz_k).^2)/(2*R*T_k)).*(((x-ux_k).^2+(y-uy_k).^2+(z-uz_k).^2)/2) )));
Jk(5,3) = - dx*dy*dz*sum(sum(sum( (n_k/(2*pi*R*T_k)^1.5)*((y-uy_k)/(R*T_k)).*exp(-((x-ux_k).^2+(y-uy_k).^2+(z-uz_k).^2)/(2*R*T_k)).*(((x-ux_k).^2+(y-uy_k).^2+(z-uz_k).^2)/2) )));
Jk(5,4) = - dx*dy*dz*sum(sum(sum( (n_k/(2*pi*R*T_k)^1.5)*((z-uz_k)/(R*T_k)).*exp(-((x-ux_k).^2+(y-uy_k).^2+(z-uz_k).^2)/(2*R*T_k)).*(((x-ux_k).^2+(y-uy_k).^2+(z-uz_k).^2)/2) )));
Jk(5,5) = - dx*dy*dz*sum(sum(sum( (n_k/(2*pi*R)^1.5)*exp(-((x-ux_k).^2+(y-uy_k).^2+(z-uz_k).^2)/(2*R*T_k)).*(-3/(2*T_k^2.5) + ((x-ux_k).^2+(y-uy_k).^2+(z-uz_k).^2)/(2*R*T_k^3.5)).*(((x-ux_k).^2+(y-uy_k).^2+(z-uz_k).^2)/2) )));

dMk = -Jk\Rk;
Mk = Mk + dMk;

n_k = Mk(1); ux_k = Mk(2); uy_k = Mk(3); uz_k = Mk(4); T_k = Mk(5);
Rk = zeros(5,1);
Rk(1,1) = rhoM - dx*dy*dz*sum(sum(sum( (n_k/(2*pi*R*T_k)^1.5)*exp(-((x-ux_k).^2+(y-uy_k).^2+(z-uz_k).^2)/(2*R*T_k)) )));
Rk(2,1) = JxM - dx*dy*dz*sum(sum(sum( (n_k/(2*pi*R*T_k)^1.5)*exp(-((x-ux_k).^2+(y-uy_k).^2+(z-uz_k).^2)/(2*R*T_k)).*x )));
Rk(3,1) = JyM - dx*dy*dz*sum(sum(sum( (n_k/(2*pi*R*T_k)^1.5)*exp(-((x-ux_k).^2+(y-uy_k).^2+(z-uz_k).^2)/(2*R*T_k)).*y )));
Rk(4,1) = JzM - dx*dy*dz*sum(sum(sum( (n_k/(2*pi*R*T_k)^1.5)*exp(-((x-ux_k).^2+(y-uy_k).^2+(z-uz_k).^2)/(2*R*T_k)).*z )));
Rk(5,1) = kM - dx*dy*dz*sum(sum(sum( (n_k/(2*pi*R*T_k)^1.5)*exp(-((x-ux_k).^2+(y-uy_k).^2+(z-uz_k).^2)/(2*R*T_k)).*(x.^2+y.^2+z.^2)/2 )));

Rk_norm = norm(Rk,2);
count = count + 1;
if count > 10
    disp(['Rk norm = ',num2str(Rk_norm)]);
    return
end
end

finf = (n_k/(2*pi*R*T_k)^1.5)*exp(-((x-ux_k).^2+(y-uy_k).^2+(z-uz_k).^2)/(2*R*T_k));

end


