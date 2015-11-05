%% References:
%%
% function [pxt,pyt] = LambdaLine(x0,y0,ArcLength,x,y,l1,l2,xi1,xi2,sgn,L,mode,options)

% Input arguments:
%   x0: x-component
%   y0: x-component of grid points - vector
%   ArcLength: % parameterization of \lambda lines: linspace(0,arclength,NumPointsOnCurve)
%   x: x-component of grid points - vector
%   y: y-component of grid points - vector
%   l1,l2,xi1,xi2: First and second eigenvalues and eigenvectors of Cauchy Green
%   sgn: Sign of the \lambda vector field
%   L: Stretching parameter
%   mode: Processor modes: "serial" OR "parallel"
%   options: ODE solver options

% Output arguments:
%   pxt: x-coordinates of \lambda lines - size: [#times,#particles]
%   pyt: y-coordinates of \lambda lines - size: [#times,#particles]
%--------------------------------------------------------------------------
% Author: Alireza Hadjighasem  alirezah@ethz.ch
% http://www.zfm.ethz.ch/~hadjighasem/index.html
%--------------------------------------------------------------------------
function [pxt,pyt] = LambdaLine(x0,y0,ArcLength,x,y,l1,l2,xi1,xi2,sgn,L,mode,options)
    Np = numel(x0);   % number of trajectories
    [xi,yi] = meshgrid(x,y);

    %% Checks:
    if Np<10;   mode = 'serial';   end;
    x0 = reshape(x0,[],1);
    y0 = reshape(y0,[],1);
    %%
    if strcmp(mode,'serial')
        [pxt,pyt] = eta_tracing(x0,y0,ArcLength,xi,yi,l1,l2,xi1,xi2,sgn,L,options);
    elseif strcmp(mode,'parallel')
%         cpu_num = min(28,Np);
        cpu_num = min(feature('numCores'),Np); 
        id = ceil( linspace(0,Np,cpu_num+1) );
        %% Opening MATLAB Pool
        poolobj = gcp('nocreate'); % If no pool, do not create new one.
        if isempty(poolobj)                                          % if parpool is not open
            parpool('local',cpu_num)
        elseif (~isempty(poolobj)) && (poolobj.NumWorkers~=cpu_num)  % if parpool is not consistent with cpu_num
            delete(gcp)
            parpool('local',cpu_num)
        end
        spmd
            Range = id(labindex)+1:id(labindex+1);
            [pxt,pyt]=eta_tracing(x0(Range),y0(Range),ArcLength,xi,yi,l1,l2,xi1,xi2,sgn,L,options);
        end
        pxt = cat(2,pxt{:});  
        pyt = cat(2,pyt{:});  
    else
        error('The computational mode is not defined. Choose between "serial" OR "parallel".');
    end
    
end

