% script to run example data for "NextWave" algorithm
% example data are raw (5 Hz) motion data from four SWIFT buoys
% the data are from a single "burst" (=512 seconds)
%
% the target location (x_out, y_out) for the prediction is arbitrary and
% can be changed by the user, but should be with in the range [0,500] for
% this example
%
% J. Thomson, 1/2026

clear; close all; clc; 

%% set up the example 
examplenum = 1; 
fixedtarget = true;
buoytarget = false;
makevideo = false;


%% prepare a video of the results (fine to skip this)
if makevideo
    if fixedtarget
        vidObj = VideoWriter(['NextWaveExample' num2str(examplenum) '_fixedtarget'] ,'MPEG-4');
    elseif buoytarget
        vidObj = VideoWriter(['NextWaveExample' num2str(examplenum) '_buoytarget'] ,'MPEG-4');
    end
    open(vidObj);
end

%% local geometry (user defined)
latorigin = 41.6878; % origin of local coordinate system (everything needs to be within a few 100 m)
lonorigin = -9.0545; % origin of local coordinate system (everything needs to be within a few 100 m)
rotation = 180;  % rotation of local coordinate system ** THIS MUST BE CONSISTENT WITH USAGE OF VELOCITY COMPONENTS **
% use rotation = 180 if using GenericCoordinateTransform.m from SWIFT codes repo

%% target location(s) for prediction: can be user defined, or a buoy from the array itself (for testing)

if fixedtarget 
    xtarget = 250; % user defined target location for prediction in local coordinate system [meters]
    ytarget = 150; % user defined target location for prediction in local coordinate system [meters]
elseif buoytarget
    targetbuoyindex = 4;  % which buoy to target (best if down-wave of others)
end

%% load example "burst" of raw data from SBG Ellipse sesnor running at 5 Hz on each buoy
skipwarmup = 200; % number of samples to skip at the start of bursts (i.e., skipping AHRS initialization)
burstend = 2700;  % number of samples defining end of burst ... usually 2742, needs to be same for all buoys
nbuoys = 4;  % examples have 4 buoys available.  Usually use 3 and test prediction against 4th one. 

flist = dir(['../ExampleData' num2str(examplenum) '/SWIFT*.mat']); % use '12-Sep-2022 07:00:00', which is burst index 92 from 'SWIFT22_DIGIFLOAT_07Sep2022-04Oct2022_reprocessedSBG.mat'

% preallocate input arrays
zin = NaN( length(skipwarmup:burstend), nbuoys);
uin = NaN( length(skipwarmup:burstend), nbuoys);
vin = NaN( length(skipwarmup:burstend), nbuoys);
tin = NaN( length(skipwarmup:burstend), nbuoys);
xin = NaN( length(skipwarmup:burstend), nbuoys);
yin = NaN( length(skipwarmup:burstend), nbuoys);

% loop through files (one for each buoy) and populate input array
% this has fragile indexing (assuming all data same size)
% this also assumes buoy system clocks are sync'd (and thus relative seconds since start of burst are consistent)
for fi=1:nbuoys

    load(['../ExampleData' num2str(examplenum)  '/' flist(fi).name])
    zin(:,fi) = sbgData.ShipMotion.heave(skipwarmup:burstend)'; % vertical displacement used to invert for wave propagation
    ztime = sbgData.ShipMotion.time_stamp(skipwarmup:burstend)'./1e6; % time since burst started (microseconds --> seconds)
    uin(:,fi) = sbgData.GpsVel.vel_e(skipwarmup:burstend)'; % lateral velocity used to invert for wave propagation
    vin(:,fi) = sbgData.GpsVel.vel_n(skipwarmup:burstend)'; % lateral velocity used to invert for wave propagation
    uvtime = sbgData.GpsVel.time_stamp(skipwarmup:burstend)'./1e6; % time since burst started (microseconds --> seconds)
    lat = sbgData.GpsPos.lat(skipwarmup:burstend)'; % position of measurement
    lon = sbgData.GpsPos.long(skipwarmup:burstend)'; % position of measurement
    lltime = sbgData.GpsPos.time_stamp(skipwarmup:burstend)'./1e6; % time since burst started (microseconds --> seconds)

    % pick a time reference (could improve by interpolating everything to common timestamp first, but probably only a few ms)
    tin(:,fi) = ztime;

    % map everything to a local coordinate system (in meters)
    [ x, y ] = GenericCoordinateTransform(lat, lon, latorigin, lonorigin, rotation);  % function from SWIFT codes repo

    xin(:,fi) = x;
    yin(:,fi) = y;

end

if fixedtarget
    ttarget = tin(:,1);  % create timestamps for the wave predictions with same sampling rate as input (for convenience)
elseif buoytarget
    xtarget = xin(:, targetbuoyindex);
    ytarget = yin(:, targetbuoyindex);
    ztarget = zin(:, targetbuoyindex); % groundtruth measurements from target buoy location (for testing)
    utarget = uin(:, targetbuoyindex); % groundtruth measurements from target buoy location (for testing)
    vtarget = vin(:, targetbuoyindex); % groundtruth measurements from target buoy location (for testing)
    ttarget = tin(:, targetbuoyindex); % groundtruth measurements from target buoy location (for testing)
    zin(:,targetbuoyindex) = []; % remove target buoy from alorithm input
    uin(:,targetbuoyindex) = []; % remove target buoy from alorithm input
    vin(:,targetbuoyindex) = []; % remove target buoy from alorithm input
    tin(:,targetbuoyindex) = []; % remove target buoy from alorithm input
    xin(:,targetbuoyindex) = []; % remove target buoy from alorithm input
    yin(:,targetbuoyindex) = []; % remove target buoy from alorithm input
end

% reverse the direction of the vertical displacement (becasue SBG sensor is mounted upside in SWIFT buoys)
zin = -zin;

fs = 1./mean(mean(diff(tin))); % raw data sampling rate (Hz)

%% Visualize time series data
f1 = figure();
ax1(1) = subplot(5,1,1);
plot(tin,xin)
ylabel('$x_{in}$')

ax1(2) = subplot(5,1,2);
plot(tin,yin)
ylabel('$y_{in}$')

ax1(3) = subplot(5,1,3);
plot(tin,zin)
ylabel('$z_{in}$')

ax1(4) = subplot(5,1,4);
plot(tin,uin)
ylabel('$u_{in}$')

ax1(5) = subplot(5,1,5);
plot(tin,vin)
xlabel('Time [s]')
ylabel('$v_{in}$')

linkaxes(ax1,'x')



%% Visualize spectra
% theta = wavespec.theta;
% Etheta = wavespec.Etheta;
% f = wavespec.f;
% theta_plot = [theta(1:180), 360];
% Etheta_plot = [Etheta(:,1:180), Etheta(:,1)];
% ff = figure();
% polarPcolor(f', theta_plot, log10(Etheta_plot'));
% colormap('parula')
% colorbar


%% load directional spectra for this example (actual processing using SBGwaves.m from "SWIFTcodes" repo)
% can be determined from single buoy or [better] an average of all buoys
% note that this background spectra should be in the nautical direction (i.e., direction FROM which waves are coming, not towards)

%exampleindex = 92; % example 1
%exampleindex = 10; % example 2 '08-Sep-2022 21:00:00': 

%load('/Volumes/Data/DigiFloat/DIGIFLOAT_Portugal/SWIFT22_DIGIFLOAT_fall2022_part1/SWIFT22_DIGIFLOAT_07Sep2022-04Oct2022_reprocessedSBG.mat')
%[Etheta theta E f dir spread spread2 spread2alt ] = SWIFTdirectionalspectra(SWIFT(exampleindex), true, true);
%wavespec.Etheta = Etheta; wavespec.theta = theta; wavespec.f = f; Hs = SWIFT(exampleindex).sigwaveheight; Dp = SWIFT(exampleindex).peakwavedirT; Tp = SWIFT(exampleindex).peakwaveperiod;
%save wavespec wavespec Hs Tp Dp

load(['../ExampleData' num2str(examplenum) '/wavespec.mat']);

Te = sum(wavespec.Etheta(:))./sum(sum(wavespec.Etheta,2) .* wavespec.f); % centroid wave period
fe = 1/Te;
depth = 95;
ce = 9.8 * Te / (2 * 3.14); % phase speed at centroid wave period

ke = wavenumber(fe,depth);
lambdae = 2*pi/ke;

markerVec = 'x*s^';
colors = colormap('lines');
thetaVec = 0:0.01:2*pi;
f2 = figure();
% for ii = 1:length(tin)
    plot(0,0,'ko')
    hold on
    grid on
    for jj = 1:size(xin,2)
        xcirc = mean(xin(:,jj)) + (lambdae/4)*cos(thetaVec);
        ycirc = mean(yin(:,jj)) + (lambdae/4)*sin(thetaVec);
        plot(xin(:,jj)/lambdae,yin(:,jj)/lambdae,'Color',colors(jj,:),'Marker',markerVec(jj),'LineStyle','-')
        plot(xcirc/lambdae,ycirc/lambdae,'r')
    end
    xlim([0 250]/lambdae)
    ylim([0 250]/lambdae)
    axis equal
    xlabel('$\frac{x}{\lambda_e}$')
    xlabel('$\frac{y}{\lambda_e}$')
    % title(strcat("time = ",num2str(tin(ii))," s"))
    pause(1/fs)
    % clf
% end


%% run algorithm for a given target location
% step through temporal windows making predictions
% not the timestamps for prediction are not required to be same sampling rate as input timeseries 
% (though are set up that way for these examples)

NTe = 10; % number of wave periods to use in the input window (usually 10)

for ti = 1:round(fs):length(tin) % for smooth results, increment the windows slowly (every 1 s, thus increment indexing by data sampling rate fs)

    inputwindow = ti + (1:(NTe * Te * fs)); % indices for the samples to be used

    if ~any(inputwindow > length(tin))

        % output times and locations (should be within a few wave periods and wavelengths of input)
        maxtargetdistance = max( max( sqrt( (xin(inputwindow,:)-max(xtarget)).^2 + (yin(inputwindow,:)-max(ytarget)).^2 ) ) ); % how far is target location from buoy array
        leadtime = maxtargetdistance / ce; % how many seconds to predict ih the future (recommend as max distance between target and buoys, divided by phase speed)
        % any more than this and the information from the farthest buoy has already propagated beyond the target
        tlast = max( max(tin(inputwindow,:)) );  % last timestamp of the input window
        tpredindices = find( ttarget >= tlast & ttarget<=(tlast + leadtime ) );
        if fixedtarget && ~isempty(tpredindices)
            tpred = ttarget(tpredindices); % timestamps for the prediction ** these could be user defined instead **
            xpred = xtarget * ones(size(tpred));  % size is [times, locations] just like input
            ypred = ytarget * ones(size(tpred));  % size is [times, locations] just like input
        elseif buoytarget && ~isempty(tpredindices)
            tpred = ttarget( tpredindices ); % timestamps for the prediction ** these could be user defined instead **
            xpred = xtarget( tpredindices )';
            ypred = ytarget( tpredindices )';
        else
            tpred = NaN; % timestamps for the prediction ** these could be user defined instead **
            xpred = Nan(size(tpred));
            ypred = Nan(size(tpred));
        end

        [prediction, reconstruction, params, t] = leastSquaresWavePropagation(zin(inputwindow,:), uin(inputwindow,:), vin(inputwindow,:), ...
            tin(inputwindow,:), xin(inputwindow,:), yin(inputwindow,:), tpred, xpred, ypred, wavespec);

        prediction = reshape(prediction,length(tpred),[]);
        zout = prediction(:,1);
        uout = prediction(:,2);
        vout = prediction(:,3);

        reconstruction = reshape(reconstruction,length(inputwindow),[]);
        zr = reconstruction(:,1:size(zin,2));
        ur = reconstruction(:,size(zin,2)+(1:size(zin,2)));
        vr = reconstruction(:,2*size(zin,2)+(1:size(zin,2)));

        f3 = figure();
        for ii = 1:nbuoys
            subplot(nbuoys,1,ii)
            plot(NaN,NaN,'-','Color',colors(1,:))
            hold on
            grid on
            plot(NaN,NaN,'k-')
            plot(NaN,NaN,'-','Color',[1 1 1]*0.7)
    
            plot(tin(inputwindow,ii),zin(inputwindow,ii),'Color',colors(ii,:))
            hold on
            grid on
            plot(tin(inputwindow,ii),zr(:,ii),'k')

            plot(tpred_true(:,ii),zpred_true(:,ii),'Color',colors(ii,:))
            plot(tpred(:,ii),zout(:,ii),'k')
        
            xlim([40 150])
            ylabel('$z$ [m]')
        
            if(ii == 1)

                title('Wave Elevation Reconstruction and Prediction')
                legend('Measured','Recreated','Prediction')
            end
        end
        xlabel('Time [s]')
    
        % velocities
        if(params.use_vel)
            f4 = figure();
            for ii = 1:nbuoys
                subplot(nbuoys,1,ii)
                plot(NaN,NaN,'-','Color',colors(1,:))
                hold on
                grid on
                plot(NaN,NaN,'k-')
                plot(NaN,NaN,'-','Color',[1 1 1]*0.7)
        
                plot(tin(inputwindow,ii),uin(inputwindow,ii),'Color',colors(ii,:))
                hold on
                grid on
                plot(tin(inputwindow,ii),ur(:,ii),'k')
                 
                xlim([40 150])
                ylabel('$u$ [m/s]')
            
                plot(tpred_true(:,ii),upred_true(:,ii),'Color',colors(ii,:))
                plot(tpred(:,ii),uout(:,ii),'k')
            end
            xlabel('Time [s]')
        
            f5 = figure();
            for ii = 1:nbuoys
                subplot(nbuoys,1,ii)
                plot(NaN,NaN,'-','Color',colors(1,:))
                hold on
                grid on
                plot(NaN,NaN,'k-')
                plot(NaN,NaN,'-','Color',[1 1 1]*0.7)
        
                plot(tin(inputwindow,ii),vin(inputwindow,ii),'Color',colors(ii,:))
                hold on
                grid on
                plot(tin(inputwindow,ii),vr(:,ii),'k')
            
                plot(tpred_true(:,ii),vpred_true(:,ii),'Color',colors(ii,:))
                plot(tpred(:,ii),vout(:,ii),'k')
            
                xlim([40 150])
                ylabel('$v$ [m/s]')
            
            end
            xlabel('Time [s]')
        end

        f6 = figure();
        plot(params.A,'ko')
        hold on
        grid on
        plot(params.amps/sqrt(2),'Color',[1 1 1]*0.7)
        plot(-params.amps/sqrt(2),'Color',[1 1 1]*0.7)
        legend('A','lims')
        xlabel('Index')
        ylabel('Amplitude')
        title('Solution Amplitudes and Limits')

        % option to rerun for a different ouput location using same solution
        %prediction = reprocess_LS_predictions(xpred,ypred,tpred,params)

        % % plot the results
        % figure(1), clf
        % 
        % % incident spectra
        % subplot(2,4,3)
        % polarPcolor(wavespec.f', wavespec.theta, log10(wavespec.Etheta'))
        % title(['Hs = ' num2str(Hs,2) ' m, Dp (FROM) = ' num2str(Dp,3) 'deg'])
        % 
        % % map
        % subplot(2,4,4)
        % plot(xin(inputwindow,:), yin(inputwindow,:),'x','linewidth',2), grid, hold on  % input positions
        %         if buoytarget
        %             plot(xtarget(inputwindow), ytarget(inputwindow),'x','linewidth',2),
        %         end
        % plot(xpred,ypred,'ko','linewidth',2,'markersize',12), hold on  % output (target) positions
        % axis([ (min(xin(:))-200) (max(xin(:))+200) (min(yin(:))-200) (max(yin(:))+200)  ]) , xlabel('x [m]'), ylabel('y [m]'), grid, axis equal
        % quiver(-150,50,-sind(Dp),-cosd(Dp),100,'filled','LineWidth',1,'color',[0 0 0])
        % %legend('buoys')
        % 
        % % input
        % subplot(6,2,1),  plot(tin(inputwindow,:),zin(inputwindow,:)), ylabel('z in [m]'), set(gca,'YLim',round([-1 1]*Hs)), title('input data')
        % subplot(6,2,3),  plot(tin(inputwindow,:),uin(inputwindow,:)), ylabel('u in [m/s]'), set(gca,'YLim',round([-1 1]*Hs/Te*6.28))
        % subplot(6,2,5), plot(tin(inputwindow,:),vin(inputwindow,:)), ylabel('v in [m/s]'), set(gca,'YLim',round([-1 1]*Hs/Te*6.28))
        % set(gca,'YTickLabel',[])
        % 
        % % reconstruction
        % subplot(6,2,7),  plot(tin(inputwindow,:),zr), ylabel('z out [m]'), set(gca,'YLim',round([-1 1]*Hs)), title('reconstructions')
        % subplot(6,2,9),  plot(tin(inputwindow,:),ur), ylabel('u out [m/s]'), set(gca,'YLim',round([-1 1]*Hs/Te*6.28))
        % subplot(6,2,11), plot(tin(inputwindow,:),vr), ylabel('v out [m/s]'), set(gca,'YLim',round([-1 1]*Hs/Te*6.28))
        % xlabel('t [s]')
        % 
        % % predictions
        % subplot(6,2,8), plot(tpred,zout,'k'), hold on, ylabel('z_p [m]'), set(gca,'YLim',round([-1 1]*Hs)), title('predictions')
        % 
        % if fixedtarget
        %     %plot(tpred, zin(tpredindices,4),'k--','linewidth',2),
        % elseif buoytarget
        %     plot(tpred, ztarget(tpredindices),'k--','linewidth',2),
        % end
        % subplot(6,2,10), plot(tpred,uout,'k'), hold on, ylabel('u_p [m/s]'), set(gca,'YLim',round([-1 1]*Hs/Te*6.28))
        % if fixedtarget
        %     %plot(tpred, uin(tpredindices,4),'k--','linewidth',2),
        % elseif buoytarget
        %     plot(tpred, utarget(tpredindices),'k--','linewidth',2),
        % end
        % subplot(6,2,12), plot(tpred,vout,'k'), hold on, ylabel('v_p [m/s]'), set(gca,'YLim',round([-1 1]*Hs/Te*6.28))
        % if fixedtarget
        %     %plot(tpred, vin(tpredindices,4),'k--','linewidth',2),
        % elseif buoytarget
        %     plot(tpred, vtarget(tpredindices),'k--','linewidth',2),
        % end
        % xlabel('t [s]')
        % 
        % if makevideo
        %     currFrame = getframe(gcf);
        %     writeVideo(vidObj,currFrame);
        % end

    end
end

if makevideo
    close(vidObj);
end




