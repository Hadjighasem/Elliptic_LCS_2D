%% References:
% [1]: G. Haller, FJ. Beron-Vera, Coherent Lagrangian vortices: the black
% holes of turbulence. Journal of Fluid Mechanics. (2013)
%%
% function [px0,py0] = OrbitDetection(pxt,pyt,dist_tol)

% Input arguments:
%   pxt: x-coordinates of \lambda lines - size: [#times,#particles]
%   pyt: y-coordinates of \lambda lines - size: [#times,#particles]
%   dist_tol: Maximum permitable distance between the two ends of a potentioal closed orbit

% Output arguments:
%   px0: An structure containing information of outermost closed orbits
%   py0: An structure containing information of all detected closed orbits
%--------------------------------------------------------------------------
% Author: Alireza Hadjighasem  alirezah@ethz.ch
% http://www.zfm.ethz.ch/~hadjighasem/index.html
%--------------------------------------------------------------------------
function [px0,py0] = OrbitDetection(pxt,pyt,dist_tol)

    [Nt,Np] = size(pxt);
    xp = pxt(1,:);       % inital positions on the poincare sections
    yp = pyt(1,:);       % inital positions on the poincare sections

    flag = bsxfun(@minus,pxt(1:end-1,:),xp).*bsxfun(@minus,pxt(2:end,:),xp) < 0;
    %--- Points crossing the plane from both sides are desirable ---%
    % Ncrossing = sum(flag,1);
    tmp  = mat2cell(flag,Nt-1,ones(Np,1));
    ind_cross = cellfun(@(x) find(x),tmp,'UniformOutput',false);
    indt_2nd_cross = ones(Np,1);
    for kk=1:Np
        if  numel(ind_cross{kk})>=2
            indt_2nd_cross(kk) = ind_cross{kk}(2);
%             hold on 
%             plot(pxt(1:indt_2nd_cross(kk)+1,kk),pyt(1:indt_2nd_cross(kk)+1,kk))
        end
    end

    mask1 = indt_2nd_cross ~= 1;   % trajectories returning to the Poincare section
    %- step1: determine the type of spiral
    ind1 = sub2ind([Nt,Np],indt_2nd_cross,(1:Np)');     % before crossing     
    ind2 = sub2ind([Nt,Np],indt_2nd_cross+1,(1:Np)');   % before crossing   

    x1 = pxt(ind1);
    x2 = pxt(ind2);

    y1 = pyt(ind1);
    y2 = pyt(ind2);

    y_cross = mask1.*( y1+(y2-y1).*(xp(:)-x1)./(x2-x1) );
    y_cross(y_cross==0) = NaN;


    % spiral_indicator = mask1.*spiral_indicator;

    y_cross_p = y_cross(1:end/2);
    y_cross_n = y_cross(end/2+1:end);

    dist = y_cross-yp(:);
    dist_p = dist(1:end/2);
    dist_n = dist(end/2+1:end);

    mask2 = (dist_p.*dist_n<0) & (abs(dist_p)<dist_tol) & (abs(dist_n)<dist_tol);

    py0 = y_cross_p+(y_cross_n-y_cross_p).*(0-dist_p)./(dist_n-dist_p);
    px0 = xp(1:end/2)';
    
    px0 = px0(mask2);
    py0 = py0(mask2);

    No = numel(py0);
    disp(sprintf('... %d closed orbits detected',No));

end

