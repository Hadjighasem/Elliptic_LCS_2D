%% References:
% [1]: D. Karrasch, F. Huhn, G. Haller, Automated detection of coherent
% Lagrangian vortices in two-dimensional unsteady flows. P Roy Soc a-Math
% Phy. (2015)
%%
% function [Xe,Ye] = DetectEllipticRegion(Xs,Ys,SingularityType,MinWedgeDist,MaxWedgeDist,Min2ndDist)

% Input arguments:
%   Xs: x-coordinates of detected singularities
%   Ys: y-coordinates of detected singularities
%   SingularityType: singularity type indicator
%   MinWedgeDist: Minimum permitable distance between a wedge pair
%   MaxWedgeDist: Maximum permitable distance between a wedge pair
%   Min2ndDist: Minimum permitable distance from the 2nd closest wedge

% Output arguments:
%   Xe: x-coordinates of wedge pairs inside elliptic regions - size: [#ellipticRegions,2]
%   Ye: y-coordinates of wedge pairs inside elliptic regions - size: [#ellipticRegions,2]
%--------------------------------------------------------------------------
% Author: Alireza Hadjighasem  alirezah@ethz.ch
% http://www.zfm.ethz.ch/~hadjighasem/index.html
%--------------------------------------------------------------------------
function [Xe,Ye] = DetectEllipticRegion(Xs,Ys,SingularityType,MinWedgeDist,MaxWedgeDist,Min2ndDist)
Ns = numel(Xs);
indwedge = find(SingularityType == -1);
Nw = numel(indwedge);
 
% pairwise distances between singularities
exhaustiveobj = ExhaustiveSearcher([Xs(indwedge),Ys(indwedge)],'Distance','euclidean');
[Id2_KNN,dKNN] = knnsearch(exhaustiveobj,[Xs(indwedge),Ys(indwedge)],'k',3,'IncludeTies',false);
 
id1_KNN = repmat(indwedge,1,2);
id2_KNN = indwedge( Id2_KNN(:,2:end) );
dKNN = dKNN(:,2:end);
 
nn = 0; ind1 = []; ind2 = [];
for kk=1:Nw
    if dKNN(kk,1)<MaxWedgeDist && dKNN(kk,1)>MinWedgeDist && dKNN(kk,2)>Min2ndDist
        nn = nn+1;
        ind1(nn,1) = id1_KNN(kk);
        ind2(nn,1) = id2_KNN(kk);
    end
end
mutual = ismember([ind1,ind2],[ind2,ind1],'rows');
 
ind = sort([ind1(mutual),ind2(mutual)],2);
ind_unique = unique(ind,'rows');
 
Xe = reshape( Xs(ind_unique), size(ind_unique));
Ye = reshape( Ys(ind_unique), size(ind_unique));
 
Ne = size(Xe,1);   % number of elliptic regions
disp(sprintf('... %d elliptic regions are identified',Ne));
 
%% plotting:
hold on
h = plot(mean(Xe,2),mean(Ye,2),'ow');
hAnnotation = get(h,'Annotation');
hLegendEntry = get(hAnnotation','LegendInformation');
set(hLegendEntry,'IconDisplayStyle','off');
