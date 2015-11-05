%% References:
%%
% function [sol_outer,sol] = Cropping_ExtractingOuter(x,y,lambda_2,pxt,pyt,Xe,Ye,L_vec)

% Input arguments:
%   x: x-component of grid points - vector
%   y: y-component of grid points - vector
%   lambda_2: Largest eigenvalue of CG
%   pxt: x-coordinates of \lambda lines - size: [#times,#particles]
%   pyt: y-coordinates of \lambda lines - size: [#times,#particles]
%   Xe: x-coordinates of wedge pairs inside elliptic regions - size: [#ellipticRegions,2]
%   Ye: y-coordinates of wedge pairs inside elliptic regions - size: [#ellipticRegions,2]
%   L_vec: Vector of \lambda-values

% Output arguments:
%   sol_outer: An structure containing information of outermost closed orbits
%   sol: An structure containing information of all detected closed orbits
%--------------------------------------------------------------------------
% Author: Alireza Hadjighasem  alirezah@ethz.ch
% http://www.zfm.ethz.ch/~hadjighasem/index.html
%--------------------------------------------------------------------------
function [sol_outer,sol] = Cropping_ExtractingOuter(x,y,lambda_2,pxt,pyt,Xe,Ye,L_vec)
    Ne = size(Xe,1);
    Npart = numel(L_vec);
    
    xc = mean(Xe,2);     % elliptic region centers
    yc = mean(Ye,2);     % elliptic region centers

    
    
    sol_outer = struct('xt',[],'yt',[],'L',[],'Region_index',[],'length',[]);  % contains the outermost orbits
    
    sol_outer.xt = cell(1,Ne);
    sol_outer.yt = cell(1,Ne);
    sol_outer.L = zeros(Ne,1); 
    sol_outer.length = zeros(Ne,1);
    sol_outer.Region_index = zeros(Ne,1);
    
    
    sol = struct('xt',[],'yt',[],'L',[],'Region_index',[],'length',[]);        % contains all closed orbits
    MaxLength = zeros(Ne,1);
    for jj=1:Npart
        Nc = size(pxt{jj},2);
        if Nc==0;  continue;  end;
        for ii=1:Nc
            % Step 1: Cropping
            xp = pxt{jj}(1,ii);
            yp = pyt{jj}(1,ii);
            sgn = (pxt{jj}(1:end-1,ii)-xp).*(pxt{jj}(2:end,ii)-xp) < 0;
            ind_cross = find(sgn==1,2,'first');
            if numel(ind_cross)<2; continue; end;
            indt_2nd_cross = ind_cross(2);
            % closing the orbit manually
            xt = [pxt{jj}(1:indt_2nd_cross,ii);xp];
            yt = [pyt{jj}(1:indt_2nd_cross,ii);yp];
            length = sum( sqrt( diff(xt).^2 + diff(yt).^2 ) );
            
            sol.xt{end+1} = xt;
            sol.yt{end+1} = yt;
            sol.length(end+1,1) = length;
            sol.L(end+1) = L_vec(jj);
            % Step 2: Finding the corresponding region index.
            for kk=1:Ne
                in_flag = inpolygon(xc(kk),yc(kk),xt,yt);
                if in_flag
                    sol.Region_index(end+1,1) = kk;
                    if length > MaxLength(kk)
                        MaxLength(kk) = max(MaxLength(kk),length);
                        sol_outer.xt{kk} = xt;
                        sol_outer.yt{kk} = yt;
                        sol_outer.length(kk) = length;
                        sol_outer.L(kk) = L_vec(jj);
                        sol_outer.Region_index(kk) = kk;
                    end
                end
            end
        end
    end
    
    NclosedOrbit = numel(sol.xt);
    disp(sprintf('... %d closed orbits in total are detected',NclosedOrbit));
    
    % Create one legend for each detected closed orbit
%     legend_str = cellfun(@(x,y) repmat(sprintf('\\lambda = %1.3f',y),size(x,2),1),...
%     cshrx,num2cell(L_vec),'UniformOutput',false);
%     legend_str = legend_str(~cellfun('isempty',legend_str));
%     legend_str = cellstr(cat(1,legend_str{:}));

    
    %- Fig.1: outer most orbits for each region
    figure; hold on;
    imagesc(x,y,log(lambda_2)); colormap('gray'); hold on; alpha(0.5)
    
    color = parula(Ne);
    for kk=1:Ne
        if isempty(sol_outer.xt{kk});  continue;  end;
        legend_str = sprintf('\\lambda = %1.3f',sol_outer.L(kk));
        hold on; plot(sol_outer.xt{kk},sol_outer.yt{kk},'r','linewidth',3,...
          'color',color(kk,:),'DisplayName',legend_str);
    end
    legend(gca,'show','location','northeastoutside')
    set(gca,'ydir','normal','fontsize',28); axis equal tight; 
    xlabel('Lon[\circ]');   ylabel('Lon[\circ]');
end
