%% References:
% [1]: M. Farazmand, G. Haller, Computing Lagrangian coherent structures
% from their variational theory. Chaos. (2012)
%%
% function [lambda_1,lambda_2,xi1,xi2,C11,C12,C22] = cgTensor(xi,yi,tspan,rho,options)

% Input arguments:
%   xi: x-component of grid points defining the input scalar field - matrix
%   yi: y-component of grid points defining the input scalar field - matrix
%   tspan: time span for advecting particles 
%   rho: auxiliary distance for computing the gradient of the flow map
%   options: options structure for ordinary differential equation solvers

% Output arguments:
%   lambda_1,lambda_2: eigenvalues of CG, 0<lambda_1<lambda_2
%   xi1,xi2: First and second eigenvectors of CG corresponding to lambda_1 and lambda_2 respectively 
%   C11,C12,C22: Cauchy-Green strain tensor blocks, C = [C11,C12;C12,C22]
%--------------------------------------------------------------------------
% Author: Mohammad Farazmand, mfaraz@mit.edu
% https://sites.google.com/site/mfarazmand/
%--------------------------------------------------------------------------
function [lambda_1,lambda_2,xi1,xi2,C11,C12,C22] = cgTensor(xi,yi,tspan,rho,mode,options)

[m, n]=size(xi);
Nrad = 4;
xt=zeros(m,n,Nrad);
yt=zeros(m,n,Nrad);
for k=1:Nrad
    xt(:,:,k) = xi + rho.x*cos( (k-1)*pi/2 );
    yt(:,:,k) = yi + rho.y*sin( (k-1)*pi/2 );
end

%% Time integration
[xt,yt] = Integrator(xt(:),yt(:),tspan,mode,options);

xt = reshape(xt, m, n, Nrad);
yt = reshape(yt, m, n, Nrad);
%% computation of eigen-values and eigen-vectors

F11 = (xt(:,:,1)-xt(:,:,3))/(2*rho.x);
F12 = (xt(:,:,2)-xt(:,:,4))/(2*rho.x);
F21 = (yt(:,:,1)-yt(:,:,3))/(2*rho.y);
F22 = (yt(:,:,2)-yt(:,:,4))/(2*rho.y);

C11 = F11.^2+F21.^2;
C12 = F11.*F12+F22.*F21;
C22 = F22.^2+F12.^2;

trC  = C11+C22;
detC = C11.*C22-C12.^2;

lambda_1 = 0.5*trC-sqrt((0.5*trC).^2-detC);
lambda_2 = 0.5*trC+sqrt((0.5*trC).^2-detC);

xi1 = zeros(m,n,2);
xi2 = zeros(m,n,2);

xi2(:,:,1) = -C12./sqrt(C12.^2+(C11-lambda_2).^2);
xi2(:,:,2) = (C11-lambda_2)./sqrt(C12.^2+(C11-lambda_2).^2);
xi1(:,:,1) =  xi2(:,:,2);
xi1(:,:,2) = -xi2(:,:,1);

end
