%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% QUESTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 1. When adding in velocities, accuracy for predicting wave elevation
%%% goes down (underpredicting) when using backslash. Why does this happen
%%% and does lsqlin help?

%%% 2. When using lsqlin, we can lose accuracy over using backslash. Does
%%% this make sense?

%%% 3. What are the constraints on velocity for lsqlin? It looks like the
%%% same amps vec is used with and without velocity included?

%%% 4. What is the overall prediction zone when you have like 4 different
%%% buoys?

%%% 5. It takes much longer to run with the period range is small. I think
%%% it's because the matrix becomes more singular so it takes longer to
%%% solve. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TO ADD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 1. Buoy movement in x. Make a random displacement from initial buoy
%%% location to see how that affects accuracy. (check)
%%% 2. Error metric (check)
%%% 3. Prediction location/error (check)
%%% 4. Maybe make 3D?
%%% 5. Clean
%%% 6. Add noise (check)
%%% 7. Implement lsqlin (check)
%%% 8. Add prediction zone lines to prediction time series plots (check)
%%% 9. Add colormap for u and v
%%% 10. Turn some plots into tiledlayout with big and little titles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear;
colorsBuoys = colormap('lines');
colorsPred = brewermap(8,'Set2');
gray = [1 1 1]* 0.7;
close all; clc; 

%% Simulation parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 5; % number of wave frequencies
numBuoys = 4; % number of buoys
addNoise = 1; % 1 = add gaussian white noise with specified SNR
if(addNoise)
    snr = 2; 
end
moveBuoy = 1; % 1 = buoy moves slightly over time, 0 = buoy is stationary
use_vel = 0; % 1 = use horizontal and vertical velocities in algorithm
numXPred = 3; % Number of locations to predict %%%%%%%%%%%%%%%%%%%%%%%%%%

%% User-defined wave parameters
h = 50; % water depth [m] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g = 9.81;

%%% Wave parameters are chosen as N random numbers between min and max
%%% values that are determined by user
% Wave periods [s]
minPeriod = 3; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxPeriod = 12; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TVec = minPeriod + (maxPeriod-minPeriod).*rand(N,1);

% Wave phases [rad]
minPhase = 0; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxPhase = 2*pi; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phiVec = minPhase + (maxPhase-minPhase).*rand(N,1);

% Wave amplitudes [m]
minAmp = 1; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxAmp = 3; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aVec = minAmp + (maxAmp-minAmp).*rand(N,1);

%% Calculated wave parameters
% Frequencies and wavelengths
omegaVec = 2*pi./TVec; % wave frequencies [rad/s]
const = omegaVec.^2*h/g;
khVec = dispersion(const);
kVec = khVec/h; % wavenumbers [m^-1]
lambdaVec = 2*pi./kVec; % wavelengths [m]

% Wave speeds
cpVec = lambdaVec./TVec; % phase velocities [m/s]
cgVec = 0.5*cpVec.*(ones(size(cpVec))+((2*khVec)./sinh(2*khVec))); % group velocities [m/s]
cgMin = min(cgVec);
cgMax = max(cgVec);

%% Define space-time
%%% Time
% Collection window
T = 2*max(TVec); % duration of collected data %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prediction zone
beta = cgMin/(cgMax - cgMin);
Ti = -beta*T;
Tf = T + beta*T;

% Full time vector
dt = min(TVec)/20; % timestep based on smallest wave period [s]
tVec = [1.1*Ti:dt:1.1*Tf]; % full time range based on Ti, Tf (wave speeds and collection duration) [s]

% Collection time vector
[B,t0ind] = min(abs(tVec)); % index in tVec where collection starts
[B,Tind] = min(abs(tVec - T)); % index in tVec where collection ends
tCollVec = tVec(t0ind:Tind); % collected data time vector
M = length(tCollVec); % number of time steps in collection vector

%%% Space
% Full space vector based on prediction zone and buoy locations [m]
xMin = cgMax*Ti;
xMax = cgMin*Tf;
xVec = (xMin*1.1 + min(lambdaVec):0.1:xMax*1.1 + max(lambdaVec)); 

% Place buoys in random locations between minimum and maximum wavelengths [m]
% Note: this are approximate positions, full ones will be extracted from xVec
buoyLocations = min(lambdaVec)+(max(lambdaVec)-min(lambdaVec)).*rand(numBuoys,1);
if(moveBuoy) % if we are allowing the buoy to move during collection period
    % If buoy moves, only let it move 1/4 of the smallest wavelength
    minDist = -min(lambdaVec)/4;
    maxDist = min(lambdaVec)/4;
else
    minDist = 0;
    maxDist = 0;
end

% Find the exact buoy position using xVec
for ii = 1:numBuoys
    [B,xiInd(ii)] = min(abs(xVec - buoyLocations(ii))); % find index for each buoy position in position vector
    xi(ii) = xVec(xiInd(ii)); % actual buoy position [m]
    xCollVec(:,ii) = xi(ii) + (minDist + (maxDist-minDist).*rand(length(tCollVec),1)); % collected data from buoy [m]
    tCollVecFull(:,ii) = tCollVec; % make tCollVecFull the same size as xCollvec
end

%%% Generate space-time matrices
[xMat,tMat] = meshgrid(xVec,tVec); % genate space and time matrices
etaMat = zeros(size(xMat)); % initialize wave elevation x-t matrix [m]
uMat = zeros(size(xMat)); % itialize horizontal velocity x-t matrix [m/s]
vMat = zeros(size(xMat)); % itialize vertical velocity x-t matrix [m/s]
for ii = 1:length(tVec) % at each time step
    for jj = 1:length(xVec) % and each location
        for kk = 1:length(omegaVec) % add all wave components together using linear wave theory (functions at end)
            etaMat(ii,jj) = etaMat(ii,jj) + etaXT(aVec(kk),kVec(kk),xMat(ii,jj),omegaVec(kk),tMat(ii,jj),phiVec(kk));
            uMat(ii,jj) = uMat(ii,jj) + uXT(aVec(kk),kVec(kk),xMat(ii,jj),omegaVec(kk),tMat(ii,jj),phiVec(kk),h,0);
            vMat(ii,jj) = vMat(ii,jj) + vXT(aVec(kk),kVec(kk),xMat(ii,jj),omegaVec(kk),tMat(ii,jj),phiVec(kk),h,0);
        end
    end
end

%% Collected data
%%% Extract buoy information from full matrices using xi_ind, t0Ind, Tind
% xi_ind is positional index where each buoy is located
% t0Ind is the time index where collection starts
% Tind is the time index where collection ends
for ii = 1:numBuoys
    eta(:,ii) = squeeze(etaMat(t0ind:Tind,xiInd(ii)));
    u(:,ii) = squeeze(uMat(t0ind:Tind,xiInd(ii)));
    v(:,ii) = squeeze(vMat(t0ind:Tind,xiInd(ii)));

    % Add noise with SNR value
    if(addNoise)
        eta(:,ii) = awgn(eta(:,ii),snr);
        u(:,ii) = awgn(u(:,ii),snr);
        v(:,ii) = awgn(v(:,ii),snr);
    end
end

%%% Plot collected data (eta and velocities)
f1 = figure();
for ii = 1:3
    ax2(ii) = subplot(3,1,ii);
    hold on
    grid on
    xlim(ax2(ii),[0 tCollVec(end)])
end

for ii = 1:numBuoys
    plot(ax2(1),tCollVec,eta(:,ii),'Color',colorsBuoys(ii,:))
    plot(ax2(2),tCollVec,u(:,ii),'Color',colorsBuoys(ii,:))
    plot(ax2(3),tCollVec,v(:,ii),'Color',colorsBuoys(ii,:))
end
title(ax2(1),'Input data')
xlabel(ax2(3),'Time [s]')
ylabel(ax2(1),'$\eta$ [m]')
ylabel(ax2(2),'u [m/s]')
ylabel(ax2(3),'v [m/s]')
linkaxes(ax2,'x')

%% Prediction
%%% Time
% Prediction vector starts at T, ends at Tf (can change) %%%%%%%%%%%%%%%%%%
tPredStart = Ti;
tPredEnd = Tf;
[B,tPredIndStart] = min(abs(tVec - tPredStart));
[B,tPredIndEnd] = min(abs(tVec - tPredEnd));
tPredVec = tVec(tPredIndStart:tPredIndEnd);
MPred = length(tPredVec);

%%% Space
% Choose numXPred equidistant locations between the mean buoy location and 
% cg_min*Tf away from the mean buoy location (can change) %%%%%%%%%%%%%%%%%
xPredLocationsApprox = linspace(mean(xi),mean(xi)+(cgMin*Tf),numXPred);

for ii = 1:numXPred
    tPredMat(:,ii) = tPredVec; % Make matrix same size as xPred (will turn to column vector later)

    [B,indLocation(ii)] = min(abs(xVec - xPredLocationsApprox(ii))); % find index in xVec for each approximate location
    xPredVec(ii) = xVec(indLocation(ii)); % actual prediction locations
    xPredMat(:,ii) = xPredVec(ii)*ones(size(tPredVec)); % make these stationary in a time series

    % Temporal prediction range for each prediction location and buoy
    for jj = 1:numBuoys
        if(xPredVec(ii) > xi(jj))
            tPredRangeUpper(jj,ii) = T + ((xPredVec(ii) - xi(jj))/cgMax);
            tPredRangeLower(jj,ii) = (xPredVec(ii) - xi(jj))/cgMin;
        else
            tPredRangeUpper(jj,ii) = T + ((xPredVec(ii) - xi(jj))/cgMin);
            tPredRangeLower(jj,ii) = (xPredVec(ii) - xi(jj))/cgMax;
        end
        if(tPredRangeLower(jj,ii) > tPredRangeUpper(jj,ii))
            tPredRangeUpper(jj,ii) = NaN;
            tPredRangeLower(jj,ii) = NaN;
        end

        [B,tPredRangeUpperInd(jj,ii)] = min(abs(tPredVec - tPredRangeUpper(jj,ii)));
        [B,tPredRangeLowerInd(jj,ii)] = min(abs(tPredVec - tPredRangeLower(jj,ii)));
    end

    % Extract true values in prediction region for eta, u, and v for
    % comparison and error analysis
    etaPred(:,ii) = etaMat(tPredIndStart:tPredIndEnd,indLocation(ii));
    uPred(:,ii) = uMat(tPredIndStart:tPredIndEnd,indLocation(ii));
    vPred(:,ii) = vMat(tPredIndStart:tPredIndEnd,indLocation(ii));
end

%% Plot x-t diagram of wave and prediction window
f2 = figure();
% imagesc(xVec,tVec,etaMat)
hold on
% grid on
% colormap('bone')

% Helpful vectors for plotting prediction zone
t1 = linspace(Ti,0,100); 
t2 = linspace(0,T,100);
t3 = linspace(T,Tf,100);

for ii = 1:numBuoys
    % for Ti < t < 0
    plot((t1-T)*cgMin + xi(ii),t1,'Color',colorsBuoys(ii,:),'linewidth',2)
    plot(t1*cgMax + xi(ii),t1,'Color',colorsBuoys(ii,:),'linewidth',2)
    
    % for 0 < t < T
    plot((t2-T)*cgMin + xi(ii),t2,'Color',colorsBuoys(ii,:),'linewidth',2)
    plot(t2*cgMin + xi(ii),t2,'Color',colorsBuoys(ii,:),'linewidth',2)
    
    % for T < t < Tf
    plot((t3-T)*cgMax + xi(ii),t3,'Color',colorsBuoys(ii,:),'linewidth',2)
    plot(t3*cgMin + xi(ii),t3,'Color',colorsBuoys(ii,:),'linewidth',2)
    
    xline(xi(ii),'--','Color',colorsBuoys(ii,:),'linewidth',1.5)
    xline(min(xCollVec(:,ii)),'.-','Color',colorsBuoys(ii,:))
    xline(max(xCollVec(:,ii)),'.-','Color',colorsBuoys(ii,:))

    xline(Ti*cgMax + xi(ii),'--','Color',colorsBuoys(ii,:))
    xline(Tf*cgMin + xi(ii),'--','Color',colorsBuoys(ii,:))

    for jj = 1:numXPred
        xline(xPredVec(jj),'--','Color',colorsPred(jj,:),'LineWidth',1.5)
        plot([xPredVec(jj) xPredVec(jj)],[tPredRangeLower(ii,jj) tPredRangeUpper(ii,jj)],'Color',colorsBuoys(ii,:),'linewidth',1.5)
    end
end

yline(T,'k--')
yline(0,'k--')
yline(Ti,'k--')
yline(Tf,'k--')
rectangle('position',[min(xPredVec) min(tPredVec) (max(xPredVec)-min(xPredVec)) (max(tPredVec)-min(tPredVec))],"EdgeColor",'c',"linewidth",1.5)

set(gca,'YDir','normal')
xlabel('x [m]')
ylabel('t [s]')

%% Build propagation matrix 
%%% Turn everything to a column vector
% Horizontal to vertical
k = kVec(:);
omega = omegaVec(:);
% Collapse multiple columns from multiple buoys to single column
x_in = xCollVec(:);
t_in = tCollVecFull(:);
z_in = eta(:);
u_in = u(:);
v_in = v(:);
% Collapse multiple columns from multiple prediction locations to single column
x_pred = xPredMat(:);
t_pred = tPredMat(:);

% Construct Propagator Matrices
phi1 = x_in*k' - t_in*omega';
phi2 = x_pred*k' - t_pred*omega';

%P1: Used to invert measured wave data (M1 x N)
%P2: Used to predict at target location/time (M2 x N)
%Note: P1 and P2 are consistent formulations, but M1 may be different than M2.
if use_vel
    P1 = [[cos(phi1),sin(phi1)];... % sines and cosines for eta
        [omega'.*cos(phi1),omega'.*sin(phi1)];... % sines and cosines for u
        [omega'.*cos(phi1),omega'.*sin(phi1)]]; % sines and cosines for v

    P2 = [[cos(phi2),sin(phi2)];...
        [omega'.*cos(phi2),omega'.*sin(phi2)];...
        [omega'.*cos(phi2),omega'.*sin(phi2)]];

    inputVec = [z_in; u_in; v_in];
    lowerEta = -amps;

    amps = aVec.'/1.4142;
    % amps = [aVec.'; aVec.']/1.4142; %%% NEED TO CHECK
else
    P1 = [cos(phi1),sin(phi1)];

    P2 = [cos(phi2),sin(phi2)];

    inputVec = z_in;
    amps = aVec.'/1.4142;
end


%% Solve for coefficients
%%% Backslash
tStart_bs = tic;
A_bs = P1 \ inputVec;
dt_bs = toc(tStart_bs);

%%% LSQ
tStart_lsq = tic;
A_lsq = lsqlin(P1,inputVec,[],[],[],[],-amps,amps,[]);
dt_lsq = toc(tStart_lsq);

%% Reconstruction
% Full reconstructed vector for all buoys and time
full_reconstructed_vec_bs = P1*A_bs;
full_reconstructed_vec_lsq = P1*A_lsq;

% Pull out just eta
eta_reconstructed_vec_bs = full_reconstructed_vec_bs(1:M*numBuoys);
eta_reconstructed_vec_lsq = full_reconstructed_vec_lsq(1:M*numBuoys);

% Separate per buoy
eta_reconstructed_bs = reshape(eta_reconstructed_vec_bs,[],numBuoys);
eta_reconstructed_lsq = reshape(eta_reconstructed_vec_lsq,[],numBuoys);

% Velocities
if(use_vel)
    % Horizontal
    u_reconstructed_vec_bs = full_reconstructed_vec_bs(M*numBuoys+1:2*M*numBuoys);
    u_reconstructed_bs = reshape(u_reconstructed_vec_bs,[],numBuoys);

    u_reconstructed_vec_lsq = full_reconstructed_vec_lsq(M*numBuoys+1:2*M*numBuoys);
    u_reconstructed_lsq = reshape(u_reconstructed_vec_lsq,[],numBuoys);

    % Vertical
    v_reconstructed_vec_bs = full_reconstructed_vec_bs(2*M*numBuoys+1:3*M*numBuoys);
    v_reconstructed_bs = reshape(v_reconstructed_vec_bs,[],numBuoys);

    v_reconstructed_vec_lsq = full_reconstructed_vec_lsq(2*M*numBuoys+1:3*M*numBuoys);
    v_reconstructed_lsq = reshape(v_reconstructed_vec_lsq,[],numBuoys);
end

%%% Plot reconstruction
% Eta
f4 = figure();
for ii = 1:numBuoys
    subplot(numBuoys,1,ii)
    plot(tCollVec,eta(:,ii),'Color',colorsBuoys(ii,:))
    hold on
    grid on
    plot(tCollVec,eta_reconstructed_bs(:,ii),'k')
    plot(tCollVec,eta_reconstructed_lsq(:,ii),'Color',gray)

    xlim([tCollVec(1) tCollVec(end)])
    ylabel('$\eta$ [m]')

    if(ii == 1)
        title('Wave Elevation Reconstruction')
        legend('Measured','Recreated - BS','Recreated - LSQ')
    end
end
xlabel('Time [s]')

% Velocities
if(use_vel)
    f5 = figure();
    for ii = 1:numBuoys
        subplot(numBuoys,1,ii)
        plot(tCollVec,u(:,ii),'Color',colorsBuoys(ii,:))
        hold on
        grid on
        plot(tCollVec,u_reconstructed_bs(:,ii),'k')
        plot(tCollVec,v_reconstructed_lsq(:,ii),'Color',gray)

        xlim([tCollVec(1) tCollVec(end)])
        ylabel('u [m/s]')
    
        if(ii == 1)
            title('Horizontal Velocity Reconstruction')
            legend('Measured','Recreated - BS','Recreated - LSQ')
        end
    end
    xlabel('Time [s]')
    
    f6 = figure();
    for ii = 1:numBuoys
        subplot(numBuoys,1,ii)
        plot(tCollVec,v(:,ii),'Color',colorsBuoys(ii,:))
        hold on
        grid on
        plot(tCollVec,v_reconstructed_bs(:,ii),'k')
        plot(tCollVec,v_reconstructed_lsq(:,ii),'Color',gray)

        xlim([tCollVec(1) tCollVec(end)])
        ylabel('v [m/s]')
    
        if(ii == 1)
            title('Vertical Velocity Reconstruction')
            legend('Measured','Recreated - BS','Recreated - LSQ')
        end
    end
    xlabel('Time [s]')
end

%% Prediction
% Full prediction vector for all buoys and time
full_prediction_vec_bs = P2*A_bs;
full_prediction_vec_lsq = P2*A_lsq;

% Pull out just eta
eta_prediction_vec_bs = full_prediction_vec_bs(1:MPred*numXPred);
eta_prediction_bs = reshape(eta_prediction_vec_bs,[],numXPred);

% Separate per prediction location
eta_prediction_vec_lsq = full_prediction_vec_lsq(1:MPred*numXPred);
eta_prediction_lsq = reshape(eta_prediction_vec_lsq,[],numXPred);

% Velocities
if(use_vel)
    % Horizontal
    u_prediction_vec_bs = full_prediction_vec_bs(MPred*numXPred+1:2*MPred*numXPred);
    u_prediction_bs = reshape(u_prediction_vec_bs,[],numXPred);

    u_prediction_vec_lsq = full_prediction_vec_lsq(MPred*numXPred+1:2*MPred*numXPred);
    u_prediction_lsq = reshape(u_prediction_vec_lsq,[],numXPred);

    % Vertical
    v_prediction_vec_bs = full_prediction_vec_bs(2*MPred*numXPred+1:3*MPred*numXPred);
    v_prediction_bs = reshape(v_prediction_vec_bs,[],numXPred);

    v_prediction_vec_lsq = full_prediction_vec_lsq(2*MPred*numXPred+1:3*MPred*numXPred);
    v_prediction_lsq = reshape(v_prediction_vec_lsq,[],numXPred);
end

% Plotting
if(numXPred <= 10)
    f7 = figure();
    for ii = 1:numXPred
        subplot(numXPred,1,ii)
        plot(tPredVec,etaPred(:,ii),'Color',colorsPred(ii,:))
        hold on
        grid on
        plot(tPredVec,eta_prediction_bs(:,ii),'k')
        plot(tPredVec,eta_prediction_lsq(:,ii),'Color',gray)
        ylabel('$\eta$ [m]')
    
        if(ii == 1)
            title('Wave Elevation Prediction')
            legend('Measured','Recreated - BS','Recreated - LSQ','Autoupdate','off')
        end

        for jj = 1:numBuoys
            xline(tPredRangeLower(jj,ii),'--','Color',colorsBuoys(jj,:),'linewidth',1.5)
            xline(tPredRangeUpper(jj,ii),'--','Color',colorsBuoys(jj,:),'linewidth',1.5)
        end
    end
    xlabel('Time [s]')
    
    if(use_vel)
        f8 = figure();
        for ii = 1:numXPred
            subplot(numXPred,1,ii)
            plot(tPredVec,uPred(:,ii),"Color",colorsPred(ii,:))
            hold on
            grid on
            plot(tPredVec,u_prediction_bs(:,ii),'k')
            plot(tPredVec,u_prediction_lsq(:,ii),'Color',gray)
            ylabel('u [m/s]')
        
            if(ii == 1)
                title('Horizontal Velocity Prediction')
                legend('Measured','Recreated - BS','Recreated - LSQ','Autoupdate','off')
            end

            for jj = 1:numBuoys
                xline(tPredRangeLower(jj,ii),'--','Color',colorsBuoys(jj,:),'linewidth',1.5)
                xline(tPredRangeUpper(jj,ii),'--','Color',colorsBuoys(jj,:),'linewidth',1.5)
            end
        end
        xlabel('Time [s]')
        
        f9 = figure();
        for ii = 1:numXPred
            subplot(numXPred,1,ii)
            plot(tPredVec,vPred(:,ii),'Color',colorsPred(ii,:))
            hold on
            grid on
            plot(tPredVec,v_prediction_bs(:,ii),'k')
            plot(tPredVec,v_prediction_lsq(:,ii),'Color',gray)
            ylabel('v [m/s]')
        
            if(ii == 1)
                title('Vertical Velocity Prediction')
                legend('Measured','Recreated - BS','Recreated - LSQ','Autoupdate','off')
            end

            for jj = 1:numBuoys
                xline(tPredRangeLower(jj,ii),'--','Color',colorsBuoys(jj,:),'linewidth',1.5)
                xline(tPredRangeUpper(jj,ii),'--','Color',colorsBuoys(jj,:),'linewidth',1.5)
            end
        end
        xlabel('Time [s]')
    end
else
    f10 = figure;
    ax10(1) = subplot(1,3,1);
    imagesc(xPredVec,tPredVec,etaPred)
    hold on
    grid on
    colormap(ax10(1),'bone')
    colorbar
    
    for ii = 1:numBuoys
        % for Ti < t < 0
        plot((t1-T)*cgMin + xi(ii),t1,'Color',colorsBuoys(ii,:),'linewidth',2)
        plot(t1*cgMax + xi(ii),t1,'Color',colorsBuoys(ii,:),'linewidth',2)
        
        % for 0 < t < T
        plot((t2-T)*cgMin + xi(ii),t2,'Color',colorsBuoys(ii,:),'linewidth',2)
        plot(t2*cgMin + xi(ii),t2,'Color',colorsBuoys(ii,:),'linewidth',2)
        
        % for T < t < Tf
        plot((t3-T)*cgMax + xi(ii),t3,'Color',colorsBuoys(ii,:),'linewidth',2)
        plot(t3*cgMin + xi(ii),t3,'Color',colorsBuoys(ii,:),'linewidth',2)
        
        xline(xi(ii),'--','Color',colorsBuoys(ii,:),'linewidth',1.5)
        xline(min(xCollVec(:,ii)),'.-','Color',colorsBuoys(ii,:))
        xline(max(xCollVec(:,ii)),'.-','Color',colorsBuoys(ii,:))
    
        xline(Ti*cgMax + xi(ii),'--','Color',colorsBuoys(ii,:))
        xline(Tf*cgMin + xi(ii),'--','Color',colorsBuoys(ii,:))
    end
    
    yline(T,'k--')
    yline(0,'k--')
    yline(Ti,'k--')
    yline(Tf,'k--')
    
    set(gca,'YDir','normal')
    xlabel('x [m]')
    ylabel('t [s]')
    title('Input data')
    xlim([min(xPredVec) max(xPredVec)])
    ylim([min(tPredVec) max(tPredVec)])

    ax10(2) = subplot(1,3,2);
    imagesc(xPredVec,tPredVec,eta_prediction_bs)
    hold on
    grid on
    colormap(ax10(2),'bone')
    colorbar %%% MAKE SAME AS AX10(1)
    
    for ii = 1:numBuoys
        % for Ti < t < 0
        plot((t1-T)*cgMin + xi(ii),t1,'Color',colorsBuoys(ii,:),'linewidth',2)
        plot(t1*cgMax + xi(ii),t1,'Color',colorsBuoys(ii,:),'linewidth',2)
        
        % for 0 < t < T
        plot((t2-T)*cgMin + xi(ii),t2,'Color',colorsBuoys(ii,:),'linewidth',2)
        plot(t2*cgMin + xi(ii),t2,'Color',colorsBuoys(ii,:),'linewidth',2)
        
        % for T < t < Tf
        plot((t3-T)*cgMax + xi(ii),t3,'Color',colorsBuoys(ii,:),'linewidth',2)
        plot(t3*cgMin + xi(ii),t3,'Color',colorsBuoys(ii,:),'linewidth',2)
        
        xline(xi(ii),'--','Color',colorsBuoys(ii,:),'linewidth',1.5)
        xline(min(xCollVec(:,ii)),'.-','Color',colorsBuoys(ii,:))
        xline(max(xCollVec(:,ii)),'.-','Color',colorsBuoys(ii,:))
    
        xline(Ti*cgMax + xi(ii),'--','Color',colorsBuoys(ii,:))
        xline(Tf*cgMin + xi(ii),'--','Color',colorsBuoys(ii,:))
    end
    
    yline(T,'k--')
    yline(0,'k--')
    yline(Ti,'k--')
    yline(Tf,'k--')
    
    set(gca,'YDir','normal')
    xlabel('x [m]')
    ylabel('t [s]')
    title('Predicted')
    xlim([min(xPredVec) max(xPredVec)])
    ylim([min(tPredVec) max(tPredVec)])

    ax10(3) = subplot(1,3,3);
    imagesc(xPredVec,tPredVec,eta_prediction_lsq)
    hold on
    grid on
    colormap(ax10(3),'parula')
    colorbar 

    for ii = 1:numBuoys
        % for Ti < t < 0
        plot((t1-T)*cgMin + xi(ii),t1,'Color',colorsBuoys(ii,:),'linewidth',2)
        plot(t1*cgMax + xi(ii),t1,'Color',colorsBuoys(ii,:),'linewidth',2)
        
        % for 0 < t < T
        plot((t2-T)*cgMin + xi(ii),t2,'Color',colorsBuoys(ii,:),'linewidth',2)
        plot(t2*cgMin + xi(ii),t2,'Color',colorsBuoys(ii,:),'linewidth',2)
        
        % for T < t < Tf
        plot((t3-T)*cgMax + xi(ii),t3,'Color',colorsBuoys(ii,:),'linewidth',2)
        plot(t3*cgMin + xi(ii),t3,'Color',colorsBuoys(ii,:),'linewidth',2)
        
        xline(xi(ii),'--','Color',colorsBuoys(ii,:),'linewidth',1.5)
        xline(min(xCollVec(:,ii)),'.-','Color',colorsBuoys(ii,:))
        xline(max(xCollVec(:,ii)),'.-','Color',colorsBuoys(ii,:))
    
        xline(Ti*cgMax + xi(ii),'--','Color',colorsBuoys(ii,:))
        xline(Tf*cgMin + xi(ii),'--','Color',colorsBuoys(ii,:))
    end
    
    yline(T,'k--')
    yline(0,'k--')
    yline(Ti,'k--')
    yline(Tf,'k--')
    
    set(gca,'YDir','normal')
    xlabel('x [m]')
    ylabel('t [s]')
    title('Difference')
    xlim([min(xPredVec) max(xPredVec)])
    ylim([min(tPredVec) max(tPredVec)])
end

%% Error
%%% Recreation
etaErrorBSRecreation = matrixError(eta,eta_reconstructed_bs);
etaErrorLSQRecreation = matrixError(eta,eta_reconstructed_lsq);

if(use_vel)
    uErrorBSRecreation = matrixError(u,u_reconstructed_bs);
    uErrorLSQRecreation = matrixError(u,u_reconstructed_lsq);
    
    vErrorBSRecreation = matrixError(v,v_reconstructed_bs);
    vErrorLSQRecreation = matrixError(v,v_reconstructed_lsq);
end

%%% Prediction
% Full
etaErrorBSPredFull = matrixError(etaPred,eta_prediction_bs);
etaErrorLSQPredFull = matrixError(etaPred,eta_prediction_lsq);

if(use_vel)
    uErrorBSPredFull = matrixError(uPred,u_prediction_bs);
    uErrorLSQPredFull = matrixError(uPred,u_prediction_lsq);

    vErrorBSPredFull = matrixError(vPred,v_prediction_bs);
    vErrorLSQPredFull = matrixError(vPred,v_prediction_lsq);
end

% In prediction zone
for ii = 1:numBuoys
    for jj = 1:numXPred
        etaErrorBSPredZone(ii,jj) = timeSeriesError(etaPred(tPredRangeLowerInd(ii,jj):tPredRangeUpperInd(ii,jj),jj),eta_prediction_bs(tPredRangeLowerInd(ii,jj):tPredRangeUpperInd(ii,jj),jj));
        etaErrorLSQPredZone(ii,jj) = timeSeriesError(etaPred(tPredRangeLowerInd(ii,jj):tPredRangeUpperInd(ii,jj),jj),eta_prediction_lsq(tPredRangeLowerInd(ii,jj):tPredRangeUpperInd(ii,jj),jj));

        if(use_vel)
            uErrorBSPredZone(ii,jj) = timeSeriesError(uPred(tPredRangeLowerInd(ii,jj):tPredRangeUpperInd(ii,jj),jj),u_prediction_bs(tPredRangeLowerInd(ii,jj):tPredRangeUpperInd(ii,jj),jj));
            uErrorLSQPredZone(ii,jj) = timeSeriesError(uPred(tPredRangeLowerInd(ii,jj):tPredRangeUpperInd(ii,jj),jj),u_prediction_lsq(tPredRangeLowerInd(ii,jj):tPredRangeUpperInd(ii,jj),jj));

            vErrorBSPredZone(ii,jj) = timeSeriesError(vPred(tPredRangeLowerInd(ii,jj):tPredRangeUpperInd(ii,jj),jj),v_prediction_bs(tPredRangeLowerInd(ii,jj):tPredRangeUpperInd(ii,jj),jj));
            vErrorLSQPredZone(ii,jj) = timeSeriesError(vPred(tPredRangeLowerInd(ii,jj):tPredRangeUpperInd(ii,jj),jj),v_prediction_lsq(tPredRangeLowerInd(ii,jj):tPredRangeUpperInd(ii,jj),jj));
        end
    end
end

% Out of prediction zone
for ii = 1:numBuoys
    for jj = 1:numXPred
        % Wave elevation
        etaPredNoZone = etaPred(:,jj);
        etaPredNoZone(tPredRangeLowerInd(ii,jj):tPredRangeUpperInd(ii,jj)) = []; % clear out entries outside the prediction range for this xPred and buoy
        
        eta_prediction_bs_no_zone = eta_prediction_bs(:,jj);
        eta_prediction_bs_no_zone(tPredRangeLowerInd(ii,jj):tPredRangeUpperInd(ii,jj)) = []; % clear out entries outside the prediction range for this xPred and buoy
        
        eta_prediction_lsq_no_zone = eta_prediction_lsq(:,jj);
        eta_prediction_lsq_no_zone(tPredRangeLowerInd(ii,jj):tPredRangeUpperInd(ii,jj)) = []; % clear out entries outside the prediction range for this xPred and buoy

        etaErrorBSPredNoZone(ii,jj) = timeSeriesError(etaPredNoZone,eta_prediction_bs_no_zone); % calc error
        etaErrorLSQPredNoZone(ii,jj) = timeSeriesError(etaPredNoZone,eta_prediction_lsq_no_zone); % calc error
        
        if(use_vel)
            % Horizontal
            uPredNoZone = uPred(:,jj);
            uPredNoZone(tPredRangeLowerInd(ii,jj):tPredRangeUpperInd(ii,jj)) = [];
        
            u_prediction_bs_no_zone = u_prediction_bs(:,jj);
            u_prediction_bs_no_zone(tPredRangeLowerInd(ii,jj):tPredRangeUpperInd(ii,jj)) = [];
        
            u_prediction_lsq_no_zone = u_prediction_lsq(:,jj);
            u_prediction_lsq_no_zone(tPredRangeLowerInd(ii,jj):tPredRangeUpperInd(ii,jj)) = [];

            uErrorBSPredNoZone(ii,jj) = timeSeriesError(uPredNoZone,u_prediction_bs_no_zone);
            uErrorLSQPredNoZone(ii,jj) = timeSeriesError(uPredNoZone,u_prediction_lsq_no_zone);
        
            % Vertical
            vPredNoZone = vPred(:,jj);
            vPredNoZone(tPredRangeLowerInd(ii,jj):tPredRangeUpperInd(ii,jj)) = [];
        
            v_prediction_bs_no_zone = v_prediction_bs(:,jj);
            v_prediction_bs_no_zone(tPredRangeLowerInd(ii,jj):tPredRangeUpperInd(ii,jj)) = [];
        
            v_prediction_lsq_no_zone = v_prediction_lsq(:,jj);
            v_prediction_lsq_no_zone(tPredRangeLowerInd(ii,jj):tPredRangeUpperInd(ii,jj)) = [];

            vErrorBSPredNoZone(ii,jj) = timeSeriesError(vPredNoZone,v_prediction_bs_no_zone);
            vErrorLSQPredNoZone(ii,jj) = timeSeriesError(vPredNoZone,v_prediction_lsq_no_zone);
        end

    end
end

%% Functions
%%% Linear wave theory functions for wave elevation and velocities
% wave elevation, eta
function eta = etaXT(A,k,x,w,t,phi)
    eta = A*cos(k*x - w*t - phi);
end

% horizontal velocity, u
function u = uXT(A,k,x,w,t,phi,h,z)
    u = A*w*(cosh(k*(z+h))/sinh(k*h))*cos(k*x - w*t - phi);
end

% vertical velocity, v
function v = vXT(A,k,x,w,t,phi,h,z)
    v = A*w*(sinh(k*(z+h))/sinh(k*h))*sin(k*x - w*t - phi);
end

%%% Error functions
% Time series
function eps = timeSeriesError(x1,x2)
    eps = norm(x1-x2,2)/norm(x1,2);
end

% Matrix
function epsMat = matrixError(x1,x2)
    for ii = 1:size(x1,2)
        epsMat(ii) = timeSeriesError(x1(:,ii),x2(:,ii));
    end
end