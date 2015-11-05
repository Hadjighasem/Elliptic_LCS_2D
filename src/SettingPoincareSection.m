%% References:
% [1]: D. Karrasch, F. Huhn, G. Haller, Automated detection of coherent
% Lagrangian vortices in two-dimensional unsteady flows. P Roy Soc a-Math
% Phy. (2015)
%%
% function [px,py] = SettingPoincareSection(x,y,Xe,Ye,Np,PsectionRadius)

% Input arguments:
%   x: x-component of grid points - vector
%   y: y-component of grid points - vector
%   Xe: x-coordinates of wedge pairs inside elliptic regions - size: [#ellipticRegions,2]
%   Ye: y-coordinates of wedge pairs inside elliptic regions - size: [#ellipticRegions,2]
%   Np: Number of points on each Poincaré section
%   PsectionRadius: Length of the Poincaré sections

% Output arguments:
%   px: x-coordinates of the Poincaré sections - cell 
%   py: y-coordinates of the Poincaré sections - cell 
%--------------------------------------------------------------------------
% Author: Alireza Hadjighasem  alirezah@ethz.ch
% http://www.zfm.ethz.ch/~hadjighasem/index.html
%--------------------------------------------------------------------------
function [Px,Py] = SettingPoincareSection(x,y,Xe,Ye,Np,PsectionRadius)
Ne = size(Xe,1);   % number of elliptic regions

X1 = mean(Xe,2);
Y1 = mean(Ye,2);

%- Setting a vertical poincare section
X2 = X1;
Y2 = Y1+PsectionRadius;

Px = cell(1,Ne);
Py = cell(1,Ne);
for jj=1:Ne
    tmpx = linspace(X1(jj),X2(jj),Np)';
    tmpy = linspace(Y1(jj),Y2(jj),Np)';
    
    % truncating values located outside the domain
    in = tmpx>min(x) & tmpx<max(x) & tmpy>min(y) & tmpy<max(y);
    Px{jj} = tmpx( in );
    Py{jj} = tmpy( in );
end

%% plotting:
hold on
plot(Px{1},Py{1},'-b','linewidth',2,'DisplayName','Poincare section');
for jj=1:Ne
    h = plot(Px{jj},Py{jj},'-b','linewidth',2);
    hAnnotation = get(h,'Annotation');
    hLegendEntry = get(hAnnotation','LegendInformation');
    set(hLegendEntry,'IconDisplayStyle','off');
end
legend off; legend(gca,'show','location','northeastoutside')
