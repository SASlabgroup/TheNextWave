clear; close all; clc; 

%% Directories
addpath(genpath('..\SWIFT-codes'))
exampleNum = 2;
dataFolder = fullfile(pwd,strcat('ExampleData',num2str(exampleNum)));

%% Initialize structures
data = struct(); % input data, temporal and spectral
params = struct(); % parameters
sol = struct(); % solutions

%% SWIFT time series data
params.latorigin = 41.6878;
params.lonorigin = -9.0545;
params.skipwarmup = 200;
params.burstend = 2700;
params.use_vel = 1;

[params, data] = loadSWIFTTemporalData(dataFolder, params, data);

params.predBuoy = 1:4;
for ii = 1:length(params.predBuoy)
    data.x_target(ii) = mean(data.x(:,params.predBuoy(ii)));
    data.y_target(ii) = mean(data.y(:,params.predBuoy(ii)));
end

%% SWIFT spectra data
params.depth = 95;
spectrumFilename = 'wavespec.mat';
[params, data] = loadSWIFTSpectralData(dataFolder, spectrumFilename, params, data);

%% Run windowed analysis
params.NTe = 10; % number of periods in the analysis window
params.dTe = 1; % slide window by half a wave period
params.predTe = 2; % how much we predict in the future (in periods)

indWindowStart = 1;
indWindowEnd = indWindowStart + round((params.NTe*params.Te*params.fs)) + 1;

predWindowStart = indWindowEnd;
predWindowEnd = predWindowStart + round(params.predTe*params.Te*params.fs) + 1;

% Run loop
count = 1;
while predWindowEnd <  length(data.t)
    %% Set/move window indices and save in data structure
    data.indWindowStart(count) = indWindowStart;
    data.indWindowEnd(count) = indWindowEnd;
    data.predWindowStart(count) = predWindowStart;
    data.predWindowEnd(count) = predWindowEnd;

    %% Input
    % input variables (windowed time series)
    t_in = data.t(indWindowStart:indWindowEnd,:);
    x_in = data.x(indWindowStart:indWindowEnd,:);
    y_in = data.y(indWindowStart:indWindowEnd,:);
    z_in = data.z(indWindowStart:indWindowEnd,:);

    if(params.use_vel)
        u_in = data.u(indWindowStart:indWindowEnd,:);
        v_in = data.v(indWindowStart:indWindowEnd,:);
    else
        u_in = NaN;
        v_in = NaN;
    end

    %% Prediction
    % prediction time
    t_pred_vec = data.t(predWindowStart:predWindowEnd);

    % prediction location
    for ii = 1:length(params.predBuoy)
        x_pred_in(:,ii) = data.x_target(ii) * ones(size(t_pred_vec));
        y_pred_in(:,ii) = data.y_target(ii) * ones(size(t_pred_vec));
        z_pred_in(:,ii) = data.z(predWindowStart:predWindowEnd,params.predBuoy(ii));
        t_pred_in(:,ii) = t_pred_vec;

        if(params.use_vel)
            u_pred_in(:,ii) = data.u(predWindowStart:predWindowEnd,params.predBuoy(ii));
            v_pred_in(:,ii) = data.v(predWindowStart:predWindowEnd,params.predBuoy(ii));
        else
            u_pred_in = NaN;
            v_pred_in = NaN;
        end
    end
    
    %% Save time series inputs into data structure
    % Input
    data.t_in(:,:,count) = t_in;
    data.x_in(:,:,count) = x_in;
    data.y_in(:,:,count) = y_in;
    data.z_in(:,:,count) = z_in;
    data.u_in(:,:,count) = u_in;
    data.v_in(:,:,count) = v_in;

    % Prediction
    data.t_pred_in(:,:,count) = t_pred_in;
    data.x_pred_in(:,:,count) = x_pred_in;
    data.y_pred_in(:,:,count) = y_pred_in;
    data.z_pred_in(:,:,count) = z_pred_in;
    data.u_pred_in(:,:,count) = u_pred_in;
    data.v_pred_in(:,:,count) = v_pred_in;    

    %% Run prediction algorithm
    sol = leastSquaresWavePropagationBNL(params, data, sol, count);

    %% Go to next window
    count = count + 1;
    indWindowStart = indWindowStart + round((params.dTe*params.Te*params.fs)) + 1;
    indWindowEnd = indWindowStart + round((params.NTe*params.Te*params.fs)) + 1;
    predWindowStart = indWindowEnd;
    predWindowEnd = predWindowStart + round(params.predTe*params.Te*params.fs) + 1;
end

%%
save(fullfile(dataFolder,"processedData.mat"),"params","data","sol")


