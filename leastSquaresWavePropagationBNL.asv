% This function uses time series of vertical displacement (and horizontal velocities) to
% generate a phase-resolved prediction of sea surface elevation (and horizontal velocities) 
% at a specified time & location using an inverse linear model.  
% The solution is constrained by the directional spectra E(f,theta)
%
% INPUT
%   z_in:  vertical displacement time series where M is the length of record and 
%          P is the number of point measurements.
%   u_in,v_in:   measurements of east, north velocities at sea surface (same size as z_in) [m/s]
%           ** the velocities are optional and generally improve prediction accuracy **
%   t_in:  timestamps of input measurements [seconds]
%   x_in,y_in:   positions (easting, northing) of measurement locations [meters]
%   wavespec: data structure containing the following fields
%          Etheta - measured directional wave spectrum (Nautical convention)
%          f - vector of wave frequencies [Hz]
%          theta - vector of wave directions [degrees from North]
%        
%   t_pred:   timestamps for the prediction [seconds]
%   x_pred,y_pred:  position (easting, northing) of target location for prediction [meters]
%
% OUTPUT
%   prediction: predicted z,u,v at location(s) x_pred, y_pred during t_pred
%       if velocities are empty inputs, velocities will not be predicted
%   reconstruction: reconstructed z,u,v motions for the input locations x_in, y_in during t_in
%       if velocities are empty inputs, velocities will not be reconstructed
%   parameters: details of the solution (amplitudes and phases)
%   t: time for computation 
%
% A. Fisher and J. Thomson, 2019-2025

function sol = leastSquaresWavePropagationBNL(params, data, sol, count)
%% Reshape input & check for consistency
%%% Extract data from structures and turn them into column vectors
% wave parameters
k = data.k_i(:);
theta = data.theta_i(:);
kx = k*sind(theta');
ky = k*cosd(theta');
kx = kx(:);
ky = ky(:);
omegaMat = sqrt(9.81*k)*ones(size(theta'));
omega = omegaMat(:);

% time series inputs
x_in_mat = squeeze(data.x_in(:,:,count));
y_in_mat = squeeze(data.y_in(:,:,count));
t_in_mat = squeeze(data.t_in(:,:,count));
z_in_mat = squeeze(data.z_in(:,:,count));
x_pred_in_mat = squeeze(data.x_pred_in(:,:,count));
y_pred_in_mat = squeeze(data.y_pred_in(:,:,count));
z_pred_in_mat = squeeze(data.z_pred_in(:,:,count));
t_pred_in_mat = squeeze(data.t_pred_in(:,:,count));

if(params.use_vel)
    u_in_mat = squeeze(data.u_in(:,:,count));
    v_in_mat = squeeze(data.v_in(:,:,count));

    u_pred_in_mat = squeeze(data.u_pred_in(:,:,count));
    v_pred_in_mat = squeeze(data.v_pred_in(:,:,count));
end

t_in_steps = size(t_in_mat,1);
t_pred_in_steps = size(t_pred_in_mat,1);

x_in = reshape(x_in_mat,[],1);
y_in = reshape(y_in_mat,[],1);
t_in = reshape(t_in_mat,[],1);
z_in = reshape(z_in_mat,[],1);
x_pred_in = reshape(x_pred_in_mat,[],1);
y_pred_in = reshape(y_pred_in_mat,[],1);
z_pred_in = reshape(z_pred_in_mat,[],1);
t_pred_in = reshape(t_pred_in_mat,[],1);
if(params.use_vel)
    u_in = reshape(u_in_mat,[],1);
    v_in = reshape(v_in_mat,[],1);

    u_pred_in = reshape(u_pred_in_mat,[],1);
    v_pred_in = reshape(v_pred_in_mat,[],1);
end

% solution limits
limsVec = data.limsMat(:);
lims = [limsVec; limsVec];

%%% Make sure vectors are the right size
N_input_pts = length(z_in);
if length(x_in) ~= N_input_pts || length(y_in) ~= N_input_pts || length(t_in) ~= N_input_pts
    error('All input vectors must be equal length')
end

N_output_pts = length(t_pred_in);
if length(x_pred_in) ~= N_output_pts || length(y_pred_in) ~= N_output_pts
    error('All output vectors must be equal length')
end

%% Construct Propagator Matrices
phi1 = x_in*kx' + y_in*ky' - t_in*omega';
phi2 = x_pred_in*kx' + y_pred_in*ky' - t_pred_in*omega';

%P1: Used to invert measured wave data (M1 x N)
%P2: Used to predict at target location/time (M2 x N)
%Note: P1 and P2 are consistent formulations, but M1 may be different than M2.

if params.use_vel
    P1 = [[cos(phi1), sin(phi1)];...
        [(kx./sqrt(kx.^2+ky.^2))'.*omega'.*cos(phi1), (kx./sqrt(kx.^2+ky.^2))'.*omega'.*sin(phi1)];...
        [(ky./sqrt(kx.^2+ky.^2))'.*omega'.*cos(phi1), (ky./sqrt(kx.^2+ky.^2))'.*omega'.*sin(phi1)]];
    
    P2 = [[cos(phi2), sin(phi2)];...
        [(kx./sqrt(kx.^2+ky.^2))'.*omega'.*cos(phi2), (kx./sqrt(kx.^2+ky.^2))'.*omega'.*sin(phi2)];...
        [(ky./sqrt(kx.^2+ky.^2))'.*omega'.*cos(phi2), (ky./sqrt(kx.^2+ky.^2))'.*omega'.*sin(phi2)]]; 
    
    inputVec = [z_in; u_in; v_in];
else
    P1 = [cos(phi1), sin(phi1)];

    P2 = [cos(phi2), sin(phi2)];

    inputVec = z_in;
end

% good = find(lims(1:numel(Ei))~=0);
P1(:,lims==0) = [];
P2(:,lims==0) = [];
lims(lims==0) = [];

%% Invert linear model to solve for unknown wave amplitudes
%%% Simple regression
% solve
t_start_bs = tic;
A_bs = P1 \ inputVec;
t_run_bs = toc(t_start_bs);

% reconstruction
reconstruction_bs_vec = P1*A_bs;
reconstruction_bs = reshape(reconstruction_bs_vec,t_in_steps,[]);
z_recon_bs = reconstruction_bs(:,1:params.nbuoys);
if(params.use_vel)
    u_recon_bs = reconstruction_bs(:,params.nbuoys+(1:params.nbuoys));
    v_recon_bs = reconstruction_bs(:,2*params.nbuoys+(1:params.nbuoys));
else
    u_recon_bs = NaN;
    v_recon_bs = NaN;
end

% prediction
prediction_bs_vec = P2*A_bs;
prediction_bs = reshape(prediction_bs_vec,t_pred_in_steps,[]);
z_pred_bs = prediction_bs(:,1:length(data.x_target));
if(params.use_vel)
    u_pred_bs = prediction_bs(:,length(data.x_target)+(1:length(data.x_target)));
    v_pred_bs = prediction_bs(:,2*length(data.x_target)+(1:length(data.x_target)));
else
    u_pred_bs = NaN;
    v_pred_bs = NaN;
end

% error
for jj = 1:params.nbuoys
    z_recon_error_bs(jj) = matrixError(z_in_mat(:,jj),z_recon_bs(:,jj));
    if(params.use_vel)
        u_recon_error_bs(jj) = matrixError(u_in_mat(:,jj),u_recon_bs(:,jj));
        v_recon_error_bs(jj) = matrixError(v_in_mat(:,jj),v_recon_bs(:,jj));
    else
        u_recon_error_bs = NaN;
        v_recon_error_bs = NaN;
    end
    
    z_pred_error_bs(jj) = matrixError(z_pred_in_mat(:,jj),z_pred_bs(:,jj));
    if(params.use_vel)
        u_pred_error_bs(jj) = matrixError(u_pred_in_mat(:,jj),u_pred_bs(:,jj));
        v_pred_error_bs(jj) = matrixError(v_pred_in_mat(:,jj),v_pred_bs(:,jj));
    else
        u_pred_error_bs(jj) = NaN;
        v_pred_error_bs(jj) = NaN;
    end
end

%%% bounded least squares (lsq)
% solve
t_start_lsq = tic;
% options = optimoptions('lsqlin','Algorithm','trust-region-reflective'); % some matlab magic here
A_lsq = lsqlin(P1,inputVec,[],[],[],[],-lims,lims,[]);  
t_run_lsq = toc(t_start_lsq);

% reconstruction
reconstruction_lsq_vec = P1*A_lsq;
reconstruction_lsq = reshape(reconstruction_lsq_vec,t_in_steps,[]);
z_recon_lsq = reconstruction_lsq(:,1:params.nbuoys);
if(params.use_vel)
    u_recon_lsq = reconstruction_lsq(:,params.nbuoys+(1:params.nbuoys));
    v_recon_lsq = reconstruction_lsq(:,2*params.nbuoys+(1:params.nbuoys));
else
    u_recon_lsq = NaN;
    v_recon_lsq = NaN;
end

% prediction
prediction_lsq_vec = P2*A_lsq;
prediction_lsq = reshape(prediction_lsq_vec,t_pred_in_steps,[]);
z_pred_lsq = prediction_lsq(:,1:length(data.x_target));
if(params.use_vel)
    u_pred_lsq = prediction_lsq(:,length(data.x_target)+(1:length(data.x_target)));
    v_pred_lsq = prediction_lsq(:,2*length(data.x_target)+(1:length(data.x_target)));
else
    u_pred_lsq = NaN;
    v_pred_lsq = NaN;
end

% error
for jj = 1:params.nbuoys
    z_recon_error_lsq(jj) = matrixError(z_in_mat(:,jj),z_recon_lsq(:,jj));
    if(params.use_vel)
        u_recon_error_lsq(jj) = matrixError(u_in_mat(:,jj),u_recon_lsq(:,jj));
        v_recon_error_lsq(jj) = matrixError(v_in_mat(:,jj),v_recon_lsq(:,jj));
    else
        u_recon_error_lsq = NaN;
        v_recon_error_lsq = NaN;
    end
    
    z_pred_error_lsq(jj) = matrixError(z_pred_in_mat(:,jj),z_pred_lsq(:,jj));
    if(params.use_vel)
        u_pred_error_lsq(jj) = matrixError(u_pred_in_mat(:,jj),u_pred_lsq(:,jj));
        v_pred_error_lsq(jj) = matrixError(v_pred_in_mat(:,jj),v_pred_lsq(:,jj));
    else
        u_pred_error_lsq(jj) = NaN;
        v_pred_error_lsq(jj) = NaN;
    end
end

%%% lims
% reconstruction
reconstruction_lims_vec = P1*lims;
reconstruction_lims = reshape(reconstruction_lims_vec,t_in_steps,[]);
z_recon_lims = reconstruction_lims(:,1:params.nbuoys);
if(params.use_vel)
    u_recon_lims = reconstruction_lims(:,params.nbuoys+(1:params.nbuoys));
    v_recon_lims = reconstruction_lims(:,2*params.nbuoys+(1:params.nbuoys));
else
    u_recon_lims = NaN;
    v_recon_lims = NaN;
end

% prediction
prediction_lims_vec = P2*lims;
prediction_lims = reshape(prediction_lims_vec,t_pred_in_steps,[]);
z_pred_lims = prediction_lims(:,1:length(data.x_target));
if(params.use_vel)
    u_pred_lims = prediction_lims(:,length(data.x_target)+(1:length(data.x_target)));
    v_pred_lims = prediction_lims(:,2*length(data.x_target)+(1:length(data.x_target)));
else
    u_pred_lims = NaN;
    v_pred_lims = NaN;
end

% error
for jj = 1:params.nbuoys
    z_recon_error_lims(jj) = matrixError(z_in_mat(:,jj),z_recon_lims(:,jj));
    if(params.use_vel)
        u_recon_error_lims(jj) = matrixError(u_in_mat(:,jj),u_recon_lims(:,jj));
        v_recon_error_lims(jj) = matrixError(v_in_mat(:,jj),v_recon_lims(:,jj));
    else
        u_recon_error_lims = NaN;
        v_recon_error_lims = NaN;
    end
    
    z_pred_error_lims(jj) = matrixError(z_pred_in_mat(:,jj),z_pred_lims(:,jj));
    if(params.use_vel)
        u_pred_error_lims(jj) = matrixError(u_pred_in_mat(:,jj),u_pred_lims(:,jj));
        v_pred_error_lims(jj) = matrixError(v_pred_in_mat(:,jj),v_pred_lims(:,jj));
    else
        u_pred_error_lims(jj) = NaN;
        v_pred_error_lims(jj) = NaN;
    end
end

%% Save parameters
% sol.P1(:,:,count) = P1;
% sol.P2(:,:,count) = P2;
% sol.phi1(:,:,count) = phi1;
% sol.phi2(:,:,count) = phi2;

% bs
sol.bs.A(:,count) = A_bs;
sol.bs.t_run(count) = t_run_bs;
sol.bs.z_recon(:,:,count) = z_recon_bs;
sol.bs.u_recon(:,:,count) = u_recon_bs;
sol.bs.v_recon(:,:,count) = v_recon_bs;
sol.bs.z_recon_error(count,:) = z_recon_error_bs;
sol.bs.u_recon_error(count,:) = u_recon_error_bs;
sol.bs.v_recon_error(count,:) = v_recon_error_bs;
sol.bs.z_pred(:,:,count) = z_pred_bs;
sol.bs.u_pred(:,:,count) = u_pred_bs;
sol.bs.v_pred(:,:,count) = v_pred_bs;
sol.bs.z_pred_error(count,:) = z_pred_error_bs;
sol.bs.u_pred_error(count,:) = u_pred_error_bs;
sol.bs.v_pred_error(count,:) = v_pred_error_bs;

% lsq
sol.lsq.A(:,count) = A_lsq;
sol.lsq.lims = lims;
sol.lsq.t_run(count) = t_run_lsq;
sol.lsq.z_recon(:,:,count) = z_recon_lsq;
sol.lsq.u_recon(:,:,count) = u_recon_lsq;
sol.lsq.v_recon(:,:,count) = v_recon_lsq;
sol.lsq.z_recon_error(count,:) = z_recon_error_lsq;
sol.lsq.u_recon_error(count,:) = u_recon_error_lsq;
sol.lsq.v_recon_error(count,:) = v_recon_error_lsq;
sol.lsq.z_pred(:,:,count) = z_pred_lsq;
sol.lsq.u_pred(:,:,count) = u_pred_lsq;
sol.lsq.v_pred(:,:,count) = v_pred_lsq;
sol.lsq.z_pred_error(count,:) = z_pred_error_lsq;
sol.lsq.u_pred_error(count,:) = u_pred_error_lsq;
sol.lsq.v_pred_error(count,:) = v_pred_error_lsq;

% lims
sol.lims.z_recon(:,:,count) = z_recon_lims;
sol.lims.u_recon(:,:,count) = u_recon_lims;
sol.lims.v_recon(:,:,count) = v_recon_lims;
sol.lims.z_recon_error(count,:) = z_recon_error_lims;
sol.lims.u_recon_error(count,:) = u_recon_error_lims;
sol.lims.v_recon_error(count,:) = v_recon_error_lims;
sol.lims.z_pred(:,:,count) = z_pred_lims;
sol.lims.u_pred(:,:,count) = u_pred_lims;
sol.lims.v_pred(:,:,count) = v_pred_lims;
sol.lims.z_pred_error(count,:) = z_pred_error_lims;
sol.lims.u_pred_error(count,:) = u_pred_error_lims;
sol.lims.v_pred_error(count,:) = v_pred_error_lims;

%% Functions
%%% Error functions
% Time series
function eps = timeSeriesError(x1,x2)
    eps = norm(x1-x2,1)/norm(x1,1);
end

% Matrix
function epsMat = matrixError(x1,x2)
    for ii = 1:size(x1,2)
        epsMat(ii) = timeSeriesError(x1(:,ii),x2(:,ii));
    end
end

end
