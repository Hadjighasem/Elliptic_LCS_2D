%% References:
% [1]: M. Farazmand, D. Blazevski, G. Haller, Shearless transport barriers
% in unsteady two-dimensional flows and maps. Physica D-Nonlinear
% Phenomena. (2014)
%%
% function [Xs,Ys,SingularityType] =  SingularityDetection(C11,C22,C12,x,y,xi1,Gfilter,R,PickTol)

% Input arguments:
%   C11,C12,C22: Cauchy-Green strain tensor elements, C = [C11,C12;C12,C22]
%   x: x-component of grid points - vector
%   y: y-component of grid points - vector
%   Gfilter: Gaussian filter flag - "true" OR "false"
%   R: Radius of circular neighbourhood around each singularity
%   PickTol: Deviation tol. of peaks from 0 & 1

% Output arguments:
%   Xs: x-coordinates of detected singularities
%   Ys: y-coordinates of detected singularities
%   SingularityType: singularity type indicator - trisector:1 , wedge:-1 , undetermined:0
%--------------------------------------------------------------------------
% Author: Alireza Hadjighasem  alirezah@ethz.ch
% http://www.zfm.ethz.ch/~hadjighasem/index.html
%--------------------------------------------------------------------------
function [Xs,Ys,SingularityType] =  SingularityDetection(C11,C12,C22,x,y,xi1,Gfilter,R,PickTol)

    Ntheta = 2884;   % number of points used to construct a circle around each singularity
                            
    %-- Recall that Cauchy–Green singularities are points, where C = I. We
    %find such points at subgrid-resolution as intersections of the
    %zero-level sets of the functions z1 = C11-C22 and z2=C12=C21, where
    %Cij denote the entries of the Cauchy– Green strain tensor
    [xi,yi] = meshgrid(x,y);
    % Calculate the two surfaces
    z1 = C11-C22;    
    z2 = C12;        
    %%%% Gaussian filtering %%%%
    if Gfilter
        g = fspecial('gaussian', [3,3], 0.5);
        z1 = conv2(z1,g,'same');
        z2 = conv2(z2,g,'same');
    end
    %%%% Visualize the two surfaces %%%%
    % surface(xi, yi, z1, 'FaceColor', [0.5 1.0 0.5], 'EdgeColor', 'none');
    % surface(xi, yi, z2, 'FaceColor', [1.0 0.5 0.0], 'EdgeColor', 'none');
    % surface(xi, yi, z3, 'FaceColor', [0.8 0.9 1.0], 'EdgeColor', 'none');
    % view(3); camlight; axis vis3d

    % Take the difference between the two surface heights and find the contour
    % where that surface is zero.
    zdiff = z1 - z2;
    C = contours(xi, yi, zdiff, [0 0]);

    % Extract the x- and y-locations from the contour matrix C.
    xL = C(1, 2:end);
    yL = C(2, 2:end);

    % Interpolate on the first surface to find z-locations for the intersection
    % line.
    zL = interp2(xi, yi, z1, xL, yL);
    % Visualize the line.
    % line(xL, yL, zL, 'Color', 'r', 'LineWidth', 1);
    %%% finding the intersection with the plane z3=0. %%%
    ind = find( zL(1:end-1).*zL(2:end) <=0 );
    % Xs = zeros(1,numel(ind)); Ys = zeros(1,numel(ind));
    % for k=1:numel(ind)
    %     k
    % Xs(k) = spline(xL(ind(k)-1:ind(k)+1),zL(ind(k)-1:ind(k)+1),xL(ind(k)))
    % Ys(k) = spline(yL(ind(k)-1:ind(k)+1),zL(ind(k)-1:ind(k)+1),yL(ind(k)))
    % end
    Xs = xL(ind)+( xL(ind+1)-xL(ind) ).*( 0-zL(ind) )./( zL(ind+1)-zL(ind) );
    Ys = yL(ind)+( yL(ind+1)-yL(ind) ).*( 0-zL(ind) )./( zL(ind+1)-zL(ind) );
    % Visualize the singularity point.
    % hold on; plot(Xs,Ys,'.r','MarkerSize',18)
    %% --- Finding the type of singularities --- %%
    Ns = numel(Xs);
    theta = linspace(0,360,Ntheta)';

    Xr = bsxfun( @plus,repmat(Xs,Ntheta,1),R*cosd(theta) );
    Yr = bsxfun( @plus,repmat(Ys,Ntheta,1),R*sind(theta) );

    v1r(:,:,1) = interp2(xi,yi,xi1(:,:,1),Xr,Yr,'cubic');
    v1r(:,:,2) = interp2(xi,yi,xi1(:,:,2),Xr,Yr,'cubic');
    vir_norm = sqrt( v1r(:,:,1).^2+v1r(:,:,2).^2 );

    fstrain  = abs( bsxfun(@times,v1r(:,:,1),cosd(theta)) + bsxfun(@times,v1r(:,:,2),sind(theta)) )./vir_norm;

    SingularityType = zeros(Ns,1);   % trisector:1 , wedge:-1 , undetermined:0
    for k = 1:Ns
        [~,locMax] = findpeaks(fstrain(:,k),theta,'MINPEAKHEIGHT',1-PickTol,...
            'MinPeakDistance',35,'Annotate','extents');   % Finding Maxs
%         figure
%         findpeaks(fstrain(:,k),theta,'MINPEAKHEIGHT',1-PickTol,...
%             'MinPeakDistance',35,'Annotate','extents');   % Finding Maxs
        if numel(locMax) ==3
            SingularityType(k) =  1;  % trisector
        elseif numel(locMax) ==2
            SingularityType(k) = -1;  % wedge
        else
            SingularityType(k) = 0;  % undetermined
        end
    end
    Xs = Xs(:);   Ys = Ys(:);
    
    disp(sprintf('... %d singularities are detected with the following types:',Ns));
    disp(sprintf('.......... %d trisectors',sum(SingularityType ==  1)));
    disp(sprintf('.......... %d wedges'    ,sum(SingularityType == -1)));
    disp(sprintf('.......... %d unknowns'  ,sum(SingularityType ==  0)));
    %% plotting:
    trC  = C11+C22;
    detC = C11.*C22-C12.^2;

%     l1 = 0.5*trC-sqrt((0.5*trC).^2-detC);
    l2 = 0.5*trC+sqrt((0.5*trC).^2-detC);

    figure
    imagesc(x,y,log(l2));  hold on; colormap gray;
    plot(Xs(SingularityType ==  1),Ys(SingularityType ==  1),'Sg',...
        'MarkerSize',9,'MarkerFaceColor','g','DisplayName','trisector')   % trisector
    plot(Xs(SingularityType == -1),Ys(SingularityType == -1),'Sr',...
        'MarkerSize',9,'MarkerFaceColor','r','DisplayName','wedge')       % wedge
    plot(Xs(SingularityType ==  0),Ys(SingularityType ==  0),'Sk',...
        'MarkerSize',9,'MarkerFaceColor','y','DisplayName','unknown')     % undetermined
    set(gca,'ydir','normal','fontsize',16); axis equal tight
    for k=1:numel(Xs)
        text(Xs(k),Ys(k),['\fontname{times}','  ',num2str(k)],...
            'FontSize',12,'color','w');
    end
    legend(gca,'show','location','northeastoutside')
    set(gca,'ydir','normal','fontsize',18); axis equal tight; 
    xlabel('Lon[\circ]');   ylabel('Lat[\circ]');
end

