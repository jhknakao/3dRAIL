% Author: J. Nakao
% Date: July 31, 2025
% 
% Parameters for each numerical test (linear advection-diffusion models)
% 
% 
function [L,CFLconstraints,diffcoefs,U,G,u_exact,A,B,C,P,xvals,yvals,zvals,dx,dy,dz,compute_error] = test_parameters(testnumber,Tf,Nx,Ny,Nz)

    % % % % % % % % % % % % % % % % % % % % % % % OUTPUTS % % % % % % % % % % % % % % % % % % % % % % % % %
    %    L = length of domain, [-L/2,L/2]
    %    CFLconstraints = bounds for CFL condition: max|f1'(u)|, max|f2'(u)|, max|f3'(u)|
    %    diffcoefs = diffusion coefficients for symmetric diffusion tensor: d1*u_xx + d2*u_yy + d3*u_zz
    %    U = initial condition: frames for one-dimensional bases (Tucker decomposition)
    %    G = initial condition: core tensor (Tucker decomposition)
    %    u_exact = exact solution for computing L1 error
    %           *Full tensor, only for verifying order of accuracy; not needed otherwise.
    %    A,B,C = flow fields: div(<A,B,C>u)
    %    P = source term: phi(x,y,z,t), stored as Tucker tensor Q x1 P1 x2 P2 x3 P3
    %    compute_error = whether or not to compute the error, i.e., store the exact solution. (0=no, 1=yes)
    % % % % % % % % % % % % % % % % % % % % % % % OUTPUTS % % % % % % % % % % % % % % % % % % % % % % % % %

if testnumber==1
    % Constant coefficient
    % u_t + u_x + u_y + u_z = (1/6)*(u_xx + u_yy + u_zz)
    % Accuracy test (rank-1)
    L = 2*pi;
    xvals = linspace(-L/2,L/2,Nx+1); dx = xvals(2)-xvals(1); xvals = xvals(1:Nx)+dx/2; %cell centers
    yvals = linspace(-L/2,L/2,Ny+1); dy = yvals(2)-yvals(1); yvals = yvals(1:Ny)+dy/2; %cell centers
    zvals = linspace(-L/2,L/2,Nz+1); dz = zvals(2)-zvals(1); zvals = zvals(1:Nz)+dz/2; %cell centers
    CFLconstraints = [1,1,1];
    d1 = 1/6; d2 = 1/6; d3 = 1/6;
    diffcoefs = [d1,d2,d3];
    U = {sin(xvals'),sin(yvals'),sin(zvals')};
    G = [1];
    [x,y,z] = meshgrid(xvals,yvals,zvals);
    % Meshgrid ordering is y-x-z, so reorder to x-y-z
    x = permute(x,[2,1,3]);
    y = permute(y,[2,1,3]);
    z = permute(z,[2,1,3]); %not technically necessary since third component is not affected
    u_exact = exp(-3*d1*Tf)*sin(x-Tf).*sin(y-Tf).*sin(z-Tf);
    compute_error = 1;
    A = @(t) {{ones(Nx,1),ones(Ny,1),ones(Nz,1)},[1]};
    B = @(t) {{ones(Nx,1),ones(Ny,1),ones(Nz,1)},[1]};
    C = @(t) {{ones(Nx,1),ones(Ny,1),ones(Nz,1)},[1]};
    Q = @(t) [0];
    P1 = @(t) zeros(Nx,1);
    P2 = @(t) zeros(Ny,1);
    P3 = @(t) zeros(Nz,1);
    P = @(t) {{P1(t),P2(t),P3(t)},Q(t)};
elseif testnumber==2
    % Constant coefficient
    % u_t + u_x + u_y + u_z = (1/6)*(u_xx + u_yy + u_zz)
    % Accuracy test (rank-2)
    L = 2*pi;
    xvals = linspace(-L/2,L/2,Nx+1); dx = xvals(2)-xvals(1); xvals = xvals(1:Nx)+dx/2; %cell centers
    yvals = linspace(-L/2,L/2,Ny+1); dy = yvals(2)-yvals(1); yvals = yvals(1:Ny)+dy/2; %cell centers
    zvals = linspace(-L/2,L/2,Nz+1); dz = zvals(2)-zvals(1); zvals = zvals(1:Nz)+dz/2; %cell centers
    CFLconstraints = [1,1,1];
    d1 = 1/6; d2 = 1/6; d3 = 1/6;
    diffcoefs = [d1,d2,d3];
    U = {[ones(Nx,1),sin(xvals'),sin(2*xvals')],[ones(Ny,1),sin(yvals'),sin(2*yvals')],[ones(Nz,1),sin(zvals'),sin(2*zvals')]};
    G = cat(3,[1,0,0;0,0,0;0,0,0],[0,0,0;0,1,0;0,0,0],[0,0,0;0,0,0;0,0,1]);
    [x,y,z] = meshgrid(xvals,yvals,zvals);
    % Meshgrid ordering is y-x-z, so reorder to x-y-z
    x = permute(x,[2,1,3]);
    y = permute(y,[2,1,3]);
    z = permute(z,[2,1,3]); %not technically necessary since third component is not affected
    u_exact = 1+exp(-3*d1*Tf)*sin(x-Tf).*sin(y-Tf).*sin(z-Tf) + exp(-3*d1*4*Tf)*sin(2*(x-Tf)).*sin(2*(y-Tf)).*sin(2*(z-Tf));
    compute_error = 1;
    A = @(t) {{ones(Nx,1),ones(Ny,1),ones(Nz,1)},[1]};
    B = @(t) {{ones(Nx,1),ones(Ny,1),ones(Nz,1)},[1]};
    C = @(t) {{ones(Nx,1),ones(Ny,1),ones(Nz,1)},[1]};
    Q = @(t) [0];
    P1 = @(t) zeros(Nx,1);
    P2 = @(t) zeros(Ny,1);
    P3 = @(t) zeros(Nz,1);
    P = @(t) {{P1(t),P2(t),P3(t)},Q(t)};
elseif testnumber==3
    % Rigid body rotation with diffusion: Rotate about vector <0,0,1>
    % u_t - yu_x + xu_y = (1/3)*(u_xx + u_yy + u_zz)
    % Accuracy test
    L = 4*pi;
    xvals = linspace(-L/2,L/2,Nx+1); dx = xvals(2)-xvals(1); xvals = xvals(1:Nx)+dx/2; %cell centers
    yvals = linspace(-L/2,L/2,Ny+1); dy = yvals(2)-yvals(1); yvals = yvals(1:Ny)+dy/2; %cell centers
    zvals = linspace(-L/2,L/2,Nz+1); dz = zvals(2)-zvals(1); zvals = zvals(1:Nz)+dz/2; %cell centers
    CFLconstraints = [L/2,L/2,0];
    d1 = 1/3; d2 = 1/3; d3 = 1/3;
    diffcoefs = [d1,d2,d3];
    U = {exp(-(xvals').^2),exp(-2*(yvals').^2),exp(-3*(zvals').^2)};
    G = [1];
    [x,y,z] = meshgrid(xvals,yvals,zvals);
    % Meshgrid ordering is y-x-z, so reorder to x-y-z
    x = permute(x,[2,1,3]);
    y = permute(y,[2,1,3]);
    z = permute(z,[2,1,3]); %not technically necessary since third component is not affected
    u_exact = exp(-(x.^2 + 2*y.^2 + 3*z.^2 + 3*d1*Tf));
    compute_error = 1;
    A = @(t) {{ones(Nx,1),-yvals',ones(Nz,1)},[1]};
    B = @(t) {{xvals',ones(Ny,1),ones(Nz,1)},[1]};
    C = @(t) {{zeros(Nx,1),zeros(Ny,1),zeros(Nz,1)},[0]};
    Q = @(t) exp(-3*d1*t)*cat(3,[9*d1,0,-16*d1;0,-2,0;-4*d1,0,0],[-36*d1,0,0;0,0,0;0,0,0]);
    P1 = @(t) [exp(-xvals'.^2),(xvals').*exp(-xvals'.^2),(xvals'.^2).*exp(-xvals'.^2)];
    P2 = @(t) [exp(-2*yvals'.^2),(yvals').*exp(-2*yvals'.^2),(yvals'.^2).*exp(-2*yvals'.^2)];
    P3 = @(t) [exp(-3*zvals'.^2),(zvals'.^2).*exp(-3*zvals'.^2)];
    P = @(t) {{P1(t),P2(t),P3(t)},Q(t)};
elseif testnumber==4
    % Rigid body rotation with diffusion: Rotate about vector <0,0,1>
    % u_t - 2yu_x + 2xu_y = (1/12)*(u_xx + u_yy + u_zz)
    % Rank test #1
    L = 4*pi;
    xvals = linspace(-L/2,L/2,Nx+1); dx = xvals(2)-xvals(1); xvals = xvals(1:Nx)+dx/2; %cell centers
    yvals = linspace(-L/2,L/2,Ny+1); dy = yvals(2)-yvals(1); yvals = yvals(1:Ny)+dy/2; %cell centers
    zvals = linspace(-L/2,L/2,Nz+1); dz = zvals(2)-zvals(1); zvals = zvals(1:Nz)+dz/2; %cell centers
    CFLconstraints = [L,L,0];
    d1 = 1/12; d2 = 1/12; d3 = 1/12;
    diffcoefs = [d1,d2,d3];
    U = {exp(-(xvals').^2),exp(-9*(yvals').^2),exp(-(zvals').^2)};
    G = [1];
    u_exact = 0; %no exact solution needed
    compute_error = 0;
    A = @(t) {{ones(Nx,1),-yvals',ones(Nz,1)},[2]};
    B = @(t) {{xvals',ones(Ny,1),ones(Nz,1)},[2]};
    C = @(t) {{zeros(Nx,1),zeros(Ny,1),zeros(Nz,1)},[0]};
    Q = @(t) [0];
    P1 = @(t) zeros(Nx,1);
    P2 = @(t) zeros(Ny,1);
    P3 = @(t) zeros(Nz,1);
    P = @(t) {{P1(t),P2(t),P3(t)},Q(t)};
elseif testnumber==5
    % Rigid body rotation with diffusion: Rotate about vector <1,1,1>
    % u_t + (-y+z)u_x + (x-z)u_y + (-x+y)u_z = (1/12)*(u_xx + u_yy + u_zz)
    % Rank test #2
    L = 4*pi;
    xvals = linspace(-L/2,L/2,Nx+1); dx = xvals(2)-xvals(1); xvals = xvals(1:Nx)+dx/2; %cell centers
    yvals = linspace(-L/2,L/2,Ny+1); dy = yvals(2)-yvals(1); yvals = yvals(1:Ny)+dy/2; %cell centers
    zvals = linspace(-L/2,L/2,Nz+1); dz = zvals(2)-zvals(1); zvals = zvals(1:Nz)+dz/2; %cell centers
    CFLconstraints = [L,L,L];
    d1 = 1/12; d2 = 1/12; d3 = 1/12;
    diffcoefs = [d1,d2,d3];
    U = {exp(-2*(xvals'-pi/2).^2),exp(-2*(yvals'+pi/2).^2),exp(-2*(zvals').^2)};
    G = [1];
    u_exact = 0; %no exact solution needed
    compute_error = 0;
    A = @(t) {{ones(Nx,1),[ones(Ny,1),yvals'],[ones(Nz,1),zvals']},cat(3,[0,-1],[1,0])};
    B = @(t) {{[ones(Nx,1),xvals'],ones(Ny,1),[ones(Nz,1),zvals']},cat(3,[0;1],[-1;0])};
    C = @(t) {{[ones(Nx,1),xvals'],[ones(Ny,1),yvals'],ones(Nz,1)},[0 1;-1 0]};
    Q = @(t) [0];
    P1 = @(t) zeros(Nx,1);
    P2 = @(t) zeros(Ny,1);
    P3 = @(t) zeros(Nz,1);
    P = @(t) {{P1(t),P2(t),P3(t)},Q(t)};
elseif testnumber==6
    % Rigid body rotation with diffusion (time-dependent flow field): Rotate about vector <0,0,1>
    % u_t - tyu_x + txu_y = (1/12)*(u_xx + u_yy + u_zz)
    % Accuracy test
    L = 4*pi;
    xvals = linspace(-L/2,L/2,Nx+1); dx = xvals(2)-xvals(1); xvals = xvals(1:Nx)+dx/2; %cell centers
    yvals = linspace(-L/2,L/2,Ny+1); dy = yvals(2)-yvals(1); yvals = yvals(1:Ny)+dy/2; %cell centers
    zvals = linspace(-L/2,L/2,Nz+1); dz = zvals(2)-zvals(1); zvals = zvals(1:Nz)+dz/2; %cell centers
    CFLconstraints = [Tf*L/2,Tf*L/2,0];
    d1 = 1/3; d2 = 1/3; d3 = 1/3;
    diffcoefs = [d1,d2,d3];
    U = {exp(-(xvals').^2),exp(-2*(yvals').^2),exp(-3*(zvals').^2)};
    G = [1];
    [x,y,z] = meshgrid(xvals,yvals,zvals);
    % Meshgrid ordering is y-x-z, so reorder to x-y-z
    x = permute(x,[2,1,3]);
    y = permute(y,[2,1,3]);
    z = permute(z,[2,1,3]); %not technically necessary since third component is not affected
    u_exact = exp(-(x.^2 + 2*y.^2 + 3*z.^2 + 3*d1*Tf));
    compute_error = 1;
    A = @(t) {{ones(Nx,1),-yvals',ones(Nz,1)},[t]};
    B = @(t) {{xvals',ones(Ny,1),ones(Nz,1)},[t]};
    C = @(t) {{zeros(Nx,1),zeros(Ny,1),zeros(Nz,1)},[0]};
    Q = @(t) exp(-3*d1*t)*cat(3,[9*d1,0,-16*d1;0,-2*t,0;-4*d1,0,0],[-36*d1,0,0;0,0,0;0,0,0]);
    P1 = @(t) [exp(-xvals'.^2),(xvals').*exp(-xvals'.^2),(xvals'.^2).*exp(-xvals'.^2)];
    P2 = @(t) [exp(-2*yvals'.^2),(yvals').*exp(-2*yvals'.^2),(yvals'.^2).*exp(-2*yvals'.^2)];
    P3 = @(t) [exp(-3*zvals'.^2),(zvals'.^2).*exp(-3*zvals'.^2)];
    P = @(t) {{P1(t),P2(t),P3(t)},Q(t)};
elseif testnumber==7
    % Rigid body rotation with diffusion (time-dependent flow field): Rotate about vector <0,0,1>
    % u_t - 2tyu_x + 2txu_y = (1/12)*(u_xx + u_yy + u_zz)
    % Rank test
    L = 4*pi;
    xvals = linspace(-L/2,L/2,Nx+1); dx = xvals(2)-xvals(1); xvals = xvals(1:Nx)+dx/2; %cell centers
    yvals = linspace(-L/2,L/2,Ny+1); dy = yvals(2)-yvals(1); yvals = yvals(1:Ny)+dy/2; %cell centers
    zvals = linspace(-L/2,L/2,Nz+1); dz = zvals(2)-zvals(1); zvals = zvals(1:Nz)+dz/2; %cell centers
    CFLconstraints = [Tf*L,Tf*L,0];
    d1 = 1/12; d2 = 1/12; d3 = 1/12;
    diffcoefs = [d1,d2,d3];
    U = {exp(-(xvals').^2),exp(-9*(yvals').^2),exp(-(zvals').^2)};
    G = [1];
    u_exact = 0; %no exact solution needed
    compute_error = 0;
    A = @(t) {{ones(Nx,1),-yvals',ones(Nz,1)},[2*t]};
    B = @(t) {{xvals',ones(Ny,1),ones(Nz,1)},[2*t]};
    C = @(t) {{zeros(Nx,1),zeros(Ny,1),zeros(Nz,1)},[0]};
    Q = @(t) [0];
    P1 = @(t) zeros(Nx,1);
    P2 = @(t) zeros(Ny,1);
    P3 = @(t) zeros(Nz,1);
    P = @(t) {{P1(t),P2(t),P3(t)},Q(t)};
elseif testnumber==8
    % Dougherty-Fokker-Planck equation
    % u_t - (xu)_x - (yu)_y - (zu)_z = (1/2)*(u_xx + u_yy + u_zz)
    % Rank test
    L = 16;
    xvals = linspace(-L/2,L/2,Nx+1); dx = xvals(2)-xvals(1); xvals = xvals(1:Nx)+dx/2; %cell centers
    yvals = linspace(-L/2,L/2,Ny+1); dy = yvals(2)-yvals(1); yvals = yvals(1:Ny)+dy/2; %cell centers
    zvals = linspace(-L/2,L/2,Nz+1); dz = zvals(2)-zvals(1); zvals = zvals(1:Nz)+dz/2; %cell centers
    CFLconstraints = [L/2,L/2,L/2];
    d1 = 1/2; d2 = 1/2; d3 = 1/2;
    diffcoefs = [d1,d2,d3];
    R = 1/6; %gas constant
    n1 = 0.8037121822811545;
    u1x = -0.3403147128006618;
    u1y = 0;
    u1z = 0;
    T1 = 0.1033754314349305;
    n2 = 4.764615814550553;
    u2x = 0.05740548475117823;
    u2y = 0;
    u2z = 0;
    T2 = 3.442950196134546;
    U = {[exp(-(xvals'-u1x).^2/(2*R*T1)),exp(-(xvals'-u2x).^2/(2*R*T2))],[exp(-(yvals'-u1y).^2/(2*R*T1)),exp(-(yvals'-u2y).^2/(2*R*T2))],[exp(-(zvals'-u1z).^2/(2*R*T1)),exp(-(zvals'-u2z).^2/(2*R*T2))]};
    G = cat(3,[n1*(2*pi*R*T1)^(-3/2),0;0,0],[0,0;0,n2*(2*pi*R*T2)^(-3/2)]);
    u_exact = 0; %no exact solution needed
        % Macroscopic quantities: R = 1/6, T = 3, u = 0, n = pi^(3/2)
        % Equilibrium distribution: exp(-x^2-y^2-z^2)
    compute_error = 0;
    A = @(t) {{-xvals',ones(Ny,1),ones(Nz,1)},[1]};
    B = @(t) {{ones(Nx,1),-yvals',ones(Nz,1)},[1]};
    C = @(t) {{ones(Nx,1),ones(Ny,1),-zvals'},[1]};
    Q = @(t) [0];
    P1 = @(t) zeros(Nx,1);
    P2 = @(t) zeros(Ny,1);
    P3 = @(t) zeros(Nz,1);
    P = @(t) {{P1(t),P2(t),P3(t)},Q(t)};
end

































end