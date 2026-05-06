%%% To Do: Add original time series data to calculate initial variance
%%% before going to the frequency domain

function [params,data] = loadSWIFTTemporalData(swiftDir, params, data)
    % Extract parameters from structure
    skipwarmup = params.skipwarmup;
    burstend = params.burstend;
    latorigin = params.latorigin;
    lonorigin = params.lonorigin;

    % Find all SWIFT data in this folder
    flist = dir(fullfile(swiftDir,"SWIFT*"));
    nbuoys = length(flist); % one data file per SWIFT buoy
    rotation = 180;
    
    % preallocate input arrays
    z = NaN(length(skipwarmup:burstend),nbuoys); % heave position
    u = NaN(length(skipwarmup:burstend),nbuoys); % EW velocity? 
    v = NaN(length(skipwarmup:burstend),nbuoys); % NS velocity?
    t = NaN(length(skipwarmup:burstend),nbuoys); % time
    x = NaN(length(skipwarmup:burstend),nbuoys); % EW position?
    y = NaN(length(skipwarmup:burstend),nbuoys); % NS position?
    lat = NaN(length(skipwarmup:burstend),nbuoys); % NS position?
    lon = NaN(length(skipwarmup:burstend),nbuoys); % NS position?

    % Use sbgData structure to get positions and velocities for each SWIFT
    % buoy
    for fi=1:nbuoys % loop through data file for each buoy
        load(fullfile(flist(fi).folder,flist(fi).name)) % load file
        z_raw(:,fi) = sbgData.ShipMotion.heave(skipwarmup:burstend)'; % vertical displacement used to invert for wave propagation
        ztime = sbgData.ShipMotion.time_stamp(skipwarmup:burstend)'./1e6; % time since burst started (microseconds --> seconds)
        u_raw(:,fi) = sbgData.GpsVel.vel_e(skipwarmup:burstend)'; % lateral velocity used to invert for wave propagation
        v_raw(:,fi) = sbgData.GpsVel.vel_n(skipwarmup:burstend)'; % lateral velocity used to invert for wave propagation
        uvtime = sbgData.GpsVel.time_stamp(skipwarmup:burstend)'./1e6; % time since burst started (microseconds --> seconds)
        lat_raw = sbgData.GpsPos.lat(skipwarmup:burstend)'; % position of measurement
        lon_raw = sbgData.GpsPos.long(skipwarmup:burstend)'; % position of measurement
        lltime = sbgData.GpsPos.time_stamp(skipwarmup:burstend)'./1e6; % time since burst started (microseconds --> seconds)
    
        % pick a time reference (could improve by interpolating everything to common timestamp first, but probably only a few ms)
        t(:,fi) = ztime;
        z(:,fi) = z_raw(:,fi); %interp1(ztime,zin_raw(:,fi),tin(:,fi));
        u(:,fi) = u_raw(:,fi); %interp1(uvtime,uin_raw(:,fi),tin(:,fi));
        v(:,fi) = v_raw(:,fi); %interp1(uvtime,vin_raw(:,fi),tin(:,fi));
        lat(:,fi) = lat_raw; %interp1(lltime,lat_raw,tin(:,fi));
        lon(:,fi) = lon_raw; %interp1(lltime,lon_raw,tin(:,fi));
    
        % map everything to a local coordinate system (in meters)
        [x_cs, y_cs] = GenericCoordinateTransform(lat(:,fi), lon(:,fi), latorigin, lonorigin, rotation);  % function from SWIFT codes repo
        x(:,fi) = x_cs;
        y(:,fi) = y_cs;

        % var(fi) = trapz(t(:,fi),z(:,fi).^2);
    end
    z = -z;
    fs = 1./mean(mean(diff(t))); % raw data sampling rate (Hz)

    %%% Save variables into structures
    data.t = t;
    data.lat = lat;
    data.lon = lon;
    data.x = x;
    data.y = y;
    data.z = z;
    data.u = u;
    data.v = v;
    params.nbuoys = nbuoys;
    params.fs = fs;
end