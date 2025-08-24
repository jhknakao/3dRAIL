% Author: J. Nakao
% Date: July 31, 2025
% 
% Parameters for each numerical test (viscous Burgers' equation)
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
    %    P = source term: phi(x,y,z,t), stored as Tucker tensor Q x1 P1 x2 P2 x3 P3
    %    compute_error = whether or not to compute the error, i.e., store the exact solution. (0=no, 1=yes)
    % % % % % % % % % % % % % % % % % % % % % % % OUTPUTS % % % % % % % % % % % % % % % % % % % % % % % % %

if testnumber==1
    % Viscous Burgers' equation
    % Accuracy test (multilinear-rank (2,2,2))
    L = 2*pi;
    xvals = linspace(-L/2,L/2,Nx+1); dx = xvals(2)-xvals(1); xvals = xvals(1:Nx)+dx/2; %cell centers
    yvals = linspace(-L/2,L/2,Ny+1); dy = yvals(2)-yvals(1); yvals = yvals(1:Ny)+dy/2; %cell centers
    zvals = linspace(-L/2,L/2,Nz+1); dz = zvals(2)-zvals(1); zvals = zvals(1:Nz)+dz/2; %cell centers
    CFLconstraints = [1,1,1];
    d1 = 1/2; d2 = 1/2; d3 = 1/2;
    diffcoefs = [d1,d2,d3];
    U = {[sin(xvals'),cos(xvals')],[sin(yvals'),cos(yvals')],[sin(zvals'),cos(zvals')]};
    G = cat(3,[-1,0;0,1],[0,1;1,0]);
    [x,y,z] = meshgrid(xvals,yvals,zvals);
    % Meshgrid ordering is y-x-z, so reorder to x-y-z
    x = permute(x,[2,1,3]);
    y = permute(y,[2,1,3]);
    z = permute(z,[2,1,3]); %not technically necessary since third component is not affected
    u_exact = exp(-3*d1*Tf)*sin(x+y+z);
    compute_error = 1;
    A = @(t) 0; %not needed for this problem
    B = @(t) 0; %not needed for this problem
    C = @(t) 0; %not needed for this problem
    Q = @(t) (3/2)*exp(-6*d1*t)*cat(3,[-1,0;0,1],[0,1;1,0]);
    P1 = @(t) [sin(2*xvals'),cos(2*xvals')];
    P2 = @(t) [sin(2*yvals'),cos(2*yvals')];
    P3 = @(t) [sin(2*zvals'),cos(2*zvals')];
    P = @(t) {{P1(t),P2(t),P3(t)},Q(t)};
elseif testnumber==2
    % Viscous Burgers' equation
    % Rank test (vary diffusion coefficient, d)
    L = 2*pi;
    xvals = linspace(-L/2,L/2,Nx+1); dx = xvals(2)-xvals(1); xvals = xvals(1:Nx)+dx/2; %cell centers
    yvals = linspace(-L/2,L/2,Ny+1); dy = yvals(2)-yvals(1); yvals = yvals(1:Ny)+dy/2; %cell centers
    zvals = linspace(-L/2,L/2,Nz+1); dz = zvals(2)-zvals(1); zvals = zvals(1:Nz)+dz/2; %cell centers
    CFLconstraints = [1,1,1];
    d1 = 1/2; d2 = 1/2; d3 = 1/2;
    diffcoefs = [d1,d2,d3];
    U = {[sin(xvals'),cos(xvals')],[sin(yvals'),cos(yvals')],[sin(zvals'),cos(zvals')]};
    G = cat(3,[-1,0;0,1],[0,1;1,0]);
    u_exact = 0; %no exact solution needed
    compute_error = 0;
    A = @(t) 0; %not needed for this problem
    B = @(t) 0; %not needed for this problem
    C = @(t) 0; %not needed for this problem
    Q = @(t) [0];
    P1 = @(t) zeros(Nx,1);
    P2 = @(t) zeros(Ny,1);
    P3 = @(t) zeros(Nz,1);
    P = @(t) {{P1(t),P2(t),P3(t)},Q(t)};
end

































end