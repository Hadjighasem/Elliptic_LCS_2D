%% References:
% [1]: G. Haller, FJ. Beron-Vera, Coherent Lagrangian vortices: the black
% holes of turbulence. Journal of Fluid Mechanics. (2013)
% [2]: M. Farazmand, G. Haller, Computing Lagrangian coherent structures
% from their variational theory. Chaos. (2012)
% [3]: A. Hadjighasem & G. Haller, Geodesic Transport barriers in Jupiter’s
% Atmosphere: A Video-Based Analysis, SIAM Review, in press (2015).
%%
% [pxt,pyt] = LambdaLine(x0,y0,shear_length,xi,yi,l1,l2,xi1,xi2,sgn,L,mode,options)

% Input arguments:
%   x0: x-component
%   y0: x-component of grid points - vector
%   ArcLength: % parameterization of \lambda lines: linspace(0,arclength,NumPointsOnCurve)
%   xi: x-component of the meshgrid - matrix
%   yi: y-component of the meshgrid - matrix
%   l1,l2,xi1,xi2: First and second eigenvalues and eigenvectors of Cauchy Green
%   sgn: Sign of the \lambda vector field
%   L: Stretching parameter
%   options: ODE solver options

% Output arguments:
%   pxt: x-coordinates of \lambda lines - size: [#times,#particles]
%   pyt: y-coordinates of \lambda lines - size: [#times,#particles]
%--------------------------------------------------------------------------
% Author: Alireza Hadjighasem  alirezah@ethz.ch
% http://www.zfm.ethz.ch/~hadjighasem/index.html
%--------------------------------------------------------------------------
function [pxt,pyt] = eta_tracing(x0,y0,ArcLength,xi,yi,l1,l2,xi1,xi2,sgn,L,options)
[m,n] = size(xi);
Np = numel(x0);
nn=0;
x = xi(1,:)';     dx = abs( x(2)-x(1) );
y = yi(:,1);      dy = abs( y(2)-y(1) );

Xmin = min(x);    Ymin = min(y);
Xmax = max(x);    Ymax = max(y);

%% -- Constructing eta vector field
eta_1 = sqrt((l2-L^2)./(l2-l1)).*xi1(:,:,1)+sgn*sqrt((L^2-l1)./(l2-l1)).*xi2(:,:,1);
eta_2 = sqrt((l2-L^2)./(l2-l1)).*xi1(:,:,2)+sgn*sqrt((L^2-l1)./(l2-l1)).*xi2(:,:,2);

eta_1( imag(eta_1)~=0 ) = NaN;
eta_2( imag(eta_2)~=0 ) = NaN;

[xi1_1b,xi1_2b] = SmoothVectorField(x0,y0);
%-- make sure that all tensor lines will launch in the same direction 
sgn_0  = sign( xi1_1b(1).*xi1_1b+xi1_2b(1).*xi1_2b );
xi1_1b = sgn_0.*xi1_1b;   xi1_1b = sgn_0.*xi1_1b;

[~,F] = ode113(@fun,ArcLength,[x0(:);y0(:)],options);
pxt = F(:,1:end/2);
pyt = F(:,end/2+1:end);

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
    function dy = fun(~,y)
        nn=nn+1;
        %- Freez the particles at the boundary
        y(1:Np,1)      = min( y(1:Np,1)      ,Xmax );
        y(1:Np,1)      = max( y(1:Np,1)      ,Xmin );
        y(Np+1:2*Np,1) = min( y(Np+1:2*Np,1) ,Ymax );
        y(Np+1:2*Np,1) = max( y(Np+1:2*Np,1) ,Ymin );

        [xi1_1,xi1_2] = SmoothVectorField(y(1:Np),y(Np+1:2*Np));

        sgn_2 = sign( xi1_1.*xi1_1b+xi1_2.*xi1_2b );

        xi1_1 = sgn_2.*xi1_1;
        xi1_2 = sgn_2.*xi1_2;

        dy = zeros(2*Np,1);     % a column vector
        dy(1:Np)      = xi1_1;
        dy(Np+1:2*Np) = xi1_2;
        
        dy(isnan(dy)) = 0;

        xi1_1b = xi1_1;
        xi1_2b = xi1_2;
    end
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
    function [v1_0,v2_0] = SmoothVectorField(x0,y0)
        %-- Get indices for 4 neighbors
        id1_UL = floor( (y0-Ymin)/dy ) + 1;   
        id2_UL = floor( (x0-Xmin)/dx ) + 1;
        ind_UL = safe_sub2ind([m,n],id1_UL,id2_UL);

        id1_UR = id1_UL;  
        id2_UR = id2_UL+1;
        ind_UR = safe_sub2ind([m,n],id1_UR,id2_UR);

        id1_DL = id1_UL+1;   
        id2_DL = id2_UL;
        ind_DL = safe_sub2ind([m,n],id1_DL,id2_DL);

        id1_DR = id1_UL+1;   
        id2_DR = id2_UL+1;
        ind_DR = safe_sub2ind([m,n],id1_DR,id2_DR);

        %-- Get vector field values for 4 neighbors
        v1_UL = eta_1(ind_UL);
        v1_UR = eta_1(ind_UR);
        v1_DL = eta_1(ind_DL);
        v1_DR = eta_1(ind_DR);

        v2_UL = eta_2(ind_UL);
        v2_UR = eta_2(ind_UR);
        v2_DL = eta_2(ind_DL);
        v2_DR = eta_2(ind_DR);
        
        %-- Local Smoothing
        sgn_1 = sign( v1_UL.*v1_UR + v2_UL.*v2_UR );
        v1_UR = sgn_1.*v1_UR;
        v2_UR = sgn_1.*v2_UR;

        sgn_1 = sign( v1_UL.*v1_DL + v2_UL.*v2_DL );
        v1_DL = sgn_1.*v1_DL;
        v2_DL = sgn_1.*v2_DL;

        sgn_1 = sign( v1_UL.*v1_DR + v2_UL.*v2_DR );
        v1_DR = sgn_1.*v1_DR;
        v2_DR = sgn_1.*v2_DR;

        %-- Bilinear interpolation
        % Bilinear interpolation for v1
        c1 = ( xi(ind_UR)-x0 )/dx;
        c2 = ( x0-xi(ind_UL) )/dx;
        c3 = ( yi(ind_DL)-y0 )/dy;
        c4 = ( y0-yi(ind_UL) )/dy;

        v1_0 = c3.*( c1.*v1_UL + c2.*v1_UR ) + c4.*( c1.*v1_DL + c2.*v1_DR );

        % Bilinear interpolation for v2
        c1 = ( xi(ind_UR)-x0 )/dx;
        c2 = ( x0-xi(ind_UL) )/dx;
        c3 = ( yi(ind_DL)-y0 )/dy;
        c4 = ( y0-yi(ind_UL) )/dy;

        v2_0 = c3.*( c1.*v2_UL + c2.*v2_UR ) + c4.*( c1.*v2_DL + c2.*v2_DR );
        
        %-- Normalizing v
        norm_v = sqrt( v1_0.^2+v2_0.^2 );
        v1_0 = v1_0./(norm_v+(norm_v==0));
        v2_0 = v2_0./(norm_v+(norm_v==0));
        
%         if any(isnan(v1_0)) || any(isnan(v1_0))
%             error('... NaN values are detected in %d iteration',nn);
%         end
    end
end

function ind = safe_sub2ind(sz, rr, cc)
      rr(rr < 1) = 1;
      rr(rr > sz(1)) = sz(1);
      cc(cc < 1) = 1;
      cc(cc > sz(2)) = sz(2);
      ind = sub2ind(sz, rr, cc);
    end

