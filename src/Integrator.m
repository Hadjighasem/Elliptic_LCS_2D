%% References:
%%
% function [xt,yt] = Integrator(x0,y0,tspan,options)

% Input arguments:
%   x0,y0: x- and y-components of the initial positions
%   tspan: time span for advecting particles 
%   options: options structure for ordinary differential equation solvers

% Output arguments:
%   xt: x-component of Lagrangian trajectories - size: [#times,#particles]
%   yt: y-component of Lagrangian trajectories - size: [#times,#particles]
%--------------------------------------------------------------------------
% Author: Alireza Hadjighasem  alirezah@ethz.ch
% http://www.zfm.ethz.ch/~hadjighasem/index.html
%--------------------------------------------------------------------------
function [xt,yt] = Integrator(x0,y0,tspan,mode,options)
Np = numel(x0);               % number of particles
%% Defining the velocity field
%-- SSH data set:
load('Ocean_geostrophic_velocity.mat','lon','lat','UT','VT','time');
u_interp = griddedInterpolant({lon,lat,time},permute(UT,[2,1,3]),'linear','none');
v_interp = griddedInterpolant({lon,lat,time},permute(VT,[2,1,3]),'linear','none');
%% Computing the final positions of the Lagrangian particles:
if strcmp(mode,'serial')
    [~,F] = ode45(@ODEfun,tspan,[x0(:);y0(:)],options,u_interp,v_interp);
    xt = F(end,1:end/2);
    yt = F(end,end/2+1:end);
elseif strcmp(mode,'parallel')
%     cpu_num = min(28,Np);
    cpu_num = min(feature('numCores'),Np);
    id = ceil( linspace(0,Np,cpu_num+1) );
    %- Opening MATLAB Pool
    poolobj = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(poolobj)                                          % if parpool is not open
        parpool('local',cpu_num)
    elseif (~isempty(poolobj)) && (poolobj.NumWorkers~=cpu_num)  % if parpool is not consistent with cpu_num
        delete(gcp)
        parpool('local',cpu_num)
    end
    spmd
        Range = id(labindex)+1:id(labindex+1);
        [~,F] = ode45(@ODEfun,tspan,[x0(Range);y0(Range)],options,u_interp,v_interp);
    end
    xt = cellfun(@(x) x(end,1:end/2),F([1:end]),'UniformOutput',false);
    yt = cellfun(@(x) x(end,end/2+1:end),F([1:end]),'UniformOutput',false);
    
    xt = cat(2,xt{:});
    yt = cat(2,yt{:});
    
else
    error('The computational mode is not defined. Choose between "serial" OR "parallel".');
end  



end

function dy = ODEfun(t,y,u_interp,v_interp)   
    Np = numel(y)/2;
    dy = zeros(2*Np,1);
    dy(1:Np,1)      = u_interp( y(1:Np,1),y(Np+1:2*Np,1),t*ones(Np,1) );
    dy(Np+1:2*Np,1) = v_interp( y(1:Np,1),y(Np+1:2*Np,1),t*ones(Np,1) );
end