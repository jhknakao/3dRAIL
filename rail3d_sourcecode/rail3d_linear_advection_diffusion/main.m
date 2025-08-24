% Author: J. Nakao
% Date: August 2025
% 
% u_t + (a1*u)_x + (a2*u)_y + (a3*u)_z = d1 u_xx + d2 u_yy + d3 u_zz + c(x,y,z,t)
% Periodic BCs
% 
% 
clear all; close all;
addpath('./tensorlab_2016-03-28/'); %tensor toolbox

%% Define parameters
Tf = 10; %final time
tol = 1.0e-4; %truncation tolerance
Nx = 100; Ny = 100; Nz = 100; %number of cells

Lambdavals = 0.9;           %single run: CFL condition, dt = lambda/(max|f1'(u)|/dx + max|f2'(u)|/dy + max|f3'(u)|/dz)
% Lambdavals = 0.5:0.02:2;    %testing order of accuracy: CFL condition, dt = lambda/(max|f1'(u)|/dx + max|f2'(u)|/dy + max|f3'(u)|/dz)

testnumber = 8;
[L,CFLconstraints,diffcoefs,U0,G0,u_exact,A,B,C,P,xvals,yvals,zvals,dx,dy,dz,compute_error] = test_parameters(testnumber,Tf,Nx,Ny,Nz);
d1 = diffcoefs(1); d2 = diffcoefs(2); d3 = diffcoefs(3);

% Spectrally accurate differentiation matrices for periodic BCs [Trefethen, 2000]
columnx = [0 0.5*(-1).^(1:Nx-1).*cot((1:Nx-1)*(2*pi*dx/(2*L)))];
Dx = (2*pi/L)*toeplitz(columnx,columnx([1 Nx:-1:2]));
columny = [0 0.5*(-1).^(1:Ny-1).*cot((1:Ny-1)*(2*pi*dy/(2*L)))];
Dy = (2*pi/L)*toeplitz(columny,columny([1 Ny:-1:2]));
columnz = [0 0.5*(-1).^(1:Nz-1).*cot((1:Nz-1)*(2*pi*dz/(2*L)))];
Dz = (2*pi/L)*toeplitz(columnz,columnz([1 Nz:-1:2]));
Dxx = (2*pi/L)^2*toeplitz([-1/(3*(2*dx/L)^2)-1/6 ...
    .5*(-1).^(2:Nx)./sin((2*pi*dx/L)*(1:Nx-1)/2).^2]);
Dyy = (2*pi/L)^2*toeplitz([-1/(3*(2*dy/L)^2)-1/6 ...
    .5*(-1).^(2:Ny)./sin((2*pi*dy/L)*(1:Ny-1)/2).^2]);
Dzz = (2*pi/L)^2*toeplitz([-1/(3*(2*dz/L)^2)-1/6 ...
    .5*(-1).^(2:Nz)./sin((2*pi*dz/L)*(1:Nz-1)/2).^2]);
Dxx = d1*Dxx; Dyy = d2*Dyy; Dzz = d3*Dzz;

errvals = zeros(numel(Lambdavals),1); %L1 error for accuracy test

tic
for k = 1:numel(Lambdavals)
disp(['Starting lambda = ',num2str(Lambdavals(k))]);
lambda = Lambdavals(k);
dt = lambda/(CFLconstraints(1)/dx + CFLconstraints(2)/dy + CFLconstraints(3)/dz);
tvals = 0:dt:Tf;
    if tvals(end)~=Tf
        tvals = [tvals,Tf];
    end
Nt = numel(tvals);

U = U0; G = G0; %initial condition

rankvals = zeros(Nt,3); %multilinear rank [r1,r2,r3]
numden = zeros(Nt,1);
bulkvelx = zeros(Nt,1);
bulkvely = zeros(Nt,1);
bulkvelz = zeros(Nt,1);
energy = zeros(Nt,1);

numden(1) = dx*dy*dz*sum(U{1},1)*tens2mat(G,1)*kron(sum(U{3},1)',sum(U{2},1)');
bulkvelx(1) = dx*dy*dz*sum(xvals'.*U{1},1)*tens2mat(G,1)*kron(sum(U{3},1)',sum(U{2},1)');
bulkvely(1) = dx*dy*dz*sum(U{1},1)*tens2mat(G,1)*kron(sum(U{3},1)',sum(yvals'.*U{2},1)');
bulkvelz(1) = dx*dy*dz*sum(U{1},1)*tens2mat(G,1)*kron(sum(zvals'.*U{3},1)',sum(U{2},1)');
energy(1) = dx*dy*dz*(sum(xvals'.^2.*U{1}/2,1)*tens2mat(G,1)*kron(sum(U{3},1)',sum(U{2},1)')...
    +sum(U{1},1)*tens2mat(G,1)*kron(sum(U{3},1)',sum(yvals'.^2.*U{2}/2,1)')...
    +sum(U{1},1)*tens2mat(G,1)*kron(sum(zvals'.^2.*U{3}/2,1)',sum(U{2},1)'));

% Macroscopic quantities
rhoM = numden(1);
JxM = bulkvelx(1);
JyM = bulkvely(1);
JzM = bulkvelz(1);
kM = energy(1);

% Weight function for LoMaC [Guo and Qiu, JSC 2024]
s = 1; %weight function should sufficiently decay at the domain boundary
w1 = exp(-s*xvals'.^2); w2 = exp(-s*yvals'.^2); w3 = exp(-s*zvals'.^2); %w=w1*w2*w3
c = (sum(xvals'.^2.*w1)*sum(w2)*sum(w3) + sum(w1)*sum(yvals'.^2.*w2)*sum(w3) + sum(w1)*sum(w2)*sum(zvals'.^2.*w3))/(sum(w1)*sum(w2)*sum(w3));
cc = dx*dy*dz*sum(w1.*[ones(Nx,1),xvals'.^2,xvals'.^4],1)*tens2mat(cat(3,[c^2,-2*c,1;-2*c,2,0;1,0,0],[-2*c,2,0;2,0,0;0,0,0],[1,0,0;0,0,0;0,0,0]),1)*kron(sum(w3.*[ones(Nz,1),zvals'.^2,zvals'.^4],1)',sum(w2.*[ones(Ny,1),yvals'.^2,yvals'.^4],1)');

[r1_n,r2_n,r3_n] = size(G);
MLR = [r1_n,r2_n,r3_n];
rankvals(1,:) = MLR;

if testnumber == 8 %Dougherty-Fokker-Planck example
% Only for full verification of equilibrium-preservation and relative entropy dissipation.
    finf = QCM(xvals,yvals,zvals,rhoM,JxM,JyM,JzM,kM); %discrete moments of equilibrium distribution equal discrete moments of IC, and then used to define Maxwellian distribution
    equilvals = zeros(Nt,1); %L1 decay to equilibrium
    relentvals = zeros(Nt,1); %relative entropy
    u_approx = lmlragen(U,G);
    equilvals(1) = dx*dy*dz*sum(sum(sum(abs(u_approx-finf))));
    relentvals(1) = dx*dy*dz*sum(sum(sum(u_approx.*log(u_approx./finf))));
end

for n = 2:Nt
disp(['Starting t = ',num2str(tvals(n))]);
dtn = tvals(n)-tvals(n-1);
tn = tvals(n-1); %current time
[U,G,MLR] = IMEX111(U,G,MLR,A,B,C,P,tn,dtn,Dx,Dy,Dz,Dxx,Dyy,Dzz,Nx,Ny,Nz,tol,w1,w2,w3,dx,dy,dz,rhoM,JxM,JyM,JzM,kM,c,cc,xvals',yvals',zvals');
rankvals(n,:) = MLR;
% disp(MLR);
if max(MLR)>25
    disp(['Multilinear rank is (',num2str(MLR),').']);
    disp('Rank is getting large; terminate script?');
    yn = input('Enter 1 for yes, 0 for no :');
    if yn == 1
        return
    end
end
numden(n) = dx*dy*dz*sum(U{1},1)*tens2mat(G,1)*kron(sum(U{3},1)',sum(U{2},1)');
bulkvelx(n) = dx*dy*dz*sum(xvals'.*U{1},1)*tens2mat(G,1)*kron(sum(U{3},1)',sum(U{2},1)');
bulkvely(n) = dx*dy*dz*sum(U{1},1)*tens2mat(G,1)*kron(sum(U{3},1)',sum(yvals'.*U{2},1)');
bulkvelz(n) = dx*dy*dz*sum(U{1},1)*tens2mat(G,1)*kron(sum(zvals'.*U{3},1)',sum(U{2},1)');
energy(n) = dx*dy*dz*(sum(xvals'.^2.*U{1}/2,1)*tens2mat(G,1)*kron(sum(U{3},1)',sum(U{2},1)')...
    +sum(U{1},1)*tens2mat(G,1)*kron(sum(U{3},1)',sum(yvals'.^2.*U{2}/2,1)')...
    +sum(U{1},1)*tens2mat(G,1)*kron(sum(zvals'.^2.*U{3}/2,1)',sum(U{2},1)'));

if testnumber == 8 %Dougherty-Fokker-Planck example
% Only for full verification of equilibrium-preservation and relative entropy dissipation.
    u_approx = lmlragen(U,G);
    equilvals(n) = dx*dy*dz*sum(sum(sum(abs(u_approx-finf))));
    relentvals(n) = dx*dy*dz*sum(sum(sum(u_approx.*log(u_approx./finf))));
end
end

if compute_error == 1
% Only for full verification of the order of accuracy.
    u_approx = lmlragen(U,G);
    errvals(k) = dx*dy*dz*sum(sum(sum(abs(u_approx - u_exact)))); %L1
end
end
toc



%% Plots

% Rank plot
figure(1);clf;plot(tvals,rankvals(:,1),'b-','linewidth',1.5);
hold on;plot(tvals,rankvals(:,2),'g-','linewidth',1.5);
hold on;plot(tvals,rankvals(:,3),'m-','linewidth',1.5);
hold on;plot(tvals,sum(rankvals,2)/3,'k-.','linewidth',1.5);
xlabel('t');ylabel('rank');axis([0,Tf,0,max(rankvals(:))+1]);
legend('r_1','r_2','r_3','(r_1+r_2+r_3)/3');

% Error plot (N/A for single runs)
figure(2);clf;loglog(Lambdavals,errvals,'k-','linewidth',1.5);
%hold on;loglog(Lambdavals(ceil(0.35*end):ceil(0.75*end)),0.008*Lambdavals(ceil(0.35*end):ceil(0.75*end)).^1,'k-.','linewidth',1.5);
hold on;loglog(Lambdavals(ceil(0.25*end):ceil(1.0*end)),0.000071*Lambdavals(ceil(0.25*end):ceil(1.0*end)).^2,'b-.','linewidth',1.5);
%hold on;loglog(Lambdavals(ceil(0.25*end):ceil(1.0*end)),0.000018*Lambdavals(ceil(0.25*end):ceil(1.0*end)).^3,'m-.','linewidth',1.5);
xlabel('\lambda');ylabel('Error');
legend('L^{\infty}','L^1','L^2','Order 3','location','northwest');
title(['N_x=N_y=N_z=',num2str(Nx)]);

% Mass
figure(3);clf;semilogy(tvals,abs(numden-numden(1))/numden(1),'k-','linewidth',1.5);
xlabel('t');ylabel('relative number density');%axis([0,Tf,0,max(rankvals(:))+1]);
%axis([0,Tf,0,1.3*max(rankvals(:))]);

% Momentum
figure(4);clf;semilogy(tvals,abs(bulkvelx-bulkvelx(1)),'k-','linewidth',1.5);
hold on;semilogy(tvals,abs(bulkvely-bulkvely(1)),'g-','linewidth',1.5);
hold on;semilogy(tvals,abs(bulkvelz-bulkvelz(1)),'m-','linewidth',1.5);
xlabel('t');ylabel('absolute bulk velocity');
legend('$|\overline{v}_x-\overline{v}_x^0|$','$|\overline{v}_y-\overline{v}_y^0|$','$|\overline{v}_z-\overline{v}_z^0|$','Interpreter','Latex','Location','southeast');

% Energy
figure(5);clf;semilogy(tvals,abs(energy-energy(1))/energy(1),'k-','linewidth',1.5);
xlabel('t');ylabel('relative energy');%axis([0,Tf,0,max(rankvals(:))+1]);
%axis([0,Tf,0,1.3*max(rankvals(:))]);

if testnumber==8
    % L1 decay, equilibrium preservation
    figure(10);clf;semilogy(tvals,equilvals,'k-','linewidth',1.5);
    xlabel('t');ylabel('||f-f^{\infty}||_1');
    % Relative entropy dissipation
    figure(11);clf;semilogy(tvals,relentvals,'k-','linewidth',1.5);
    xlabel('t');ylabel('relative entropy');
end


