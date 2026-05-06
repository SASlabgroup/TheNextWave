function [params, data] = loadSWIFTSpectralData(swiftDir, spectrumFilename, params, data)
    %% Load and extract wavespec data
    load(fullfile(swiftDir,spectrumFilename),'wavespec') % wavespec: data structure of processed SWIFT data
    Etheta_orig = wavespec.Etheta; % measured direction wave spectrum (Nautical convention) [m^2/Hz/deg]
    theta_orig = wavespec.theta; % vector of wave directions [deg from North]
    f_orig = wavespec.f; % vector of wave frequencies [Hz]

    % calculate variance
    var_orig = trapz(f_orig',trapz(theta_orig,Etheta_orig'));

    %% Calculate bulk parameters @ centroid
    Te = sum(Etheta_orig(:))./sum(sum(Etheta_orig,2) .* f_orig); % centroid wave period
    fe = 1/Te; % centroid wave frequency
    ke = wavenumber(fe,params.depth);
    lambdae = 2*pi/ke;

    %% Convert wave spectrum to Cartesian coordinates (e.g. direction waves are moving TOWARDS)
    if size(Etheta_orig,1)==length(theta_orig)
        Etheta_orig = Etheta_orig'; % get the matrix in the correct orientation
    end
    [theta_orig,I] = unique(theta_orig,'last'); % find the unique elements of theta, but keep the last entry (rather than the first which is default)
    Etheta_orig = Etheta_orig(:,I); % keep unique indices in Etheta, only have to change the indices for the theta indices, not f

    theta_orig_rotated = theta_orig + 180; % changes the coordinate system (N becomes S? I think this is where theta changes from where they come from to where they are moving towards)
    theta_orig_rotated(theta_orig_rotated>360) = theta_orig_rotated(theta_orig_rotated>360) - 360; % wrap around 
    [~,I] = sort(theta_orig_rotated); % re-sort
    theta_orig_rotated_sorted = theta_orig_rotated(I); %%%%%%% NOT IN ORIGINAL ALGORITHM %%%%%%%%%%%%%%%%%%
    Etheta_orig_rotated_sorted = Etheta_orig(:,I); % sort Etheta based on new coordinate system of theta (don't need to change f indices)
        
    % calculate variance
    var_orig_rotated_sorted = trapz(f_orig',trapz(theta_orig_rotated_sorted,Etheta_orig_rotated_sorted'));

    %% Limit solution space to the important frequency and direction range
    % find frequencies where the energy in that band (?) is above 5\% of the max energy
    df = gradient(f_orig(:)); % basically just the frequency step, but I think gradient helps keep the vector the right size (and diff doesn't)
    Ef_rotated_sorted = trapz(theta_orig_rotated_sorted,Etheta_orig_rotated_sorted'); % E(f), [m^2/Hz], this is integrated over a different theta vector in the original algorithm
    frange = find((df.*Ef_rotated_sorted')./max(df.*Ef_rotated_sorted') >= 0.05); % only keep the frequencies where the wave energy spectrum is greater than 5% of the max energy (creates a cutoff for important frequencies)
    f_i = logspace(log10(f_orig(frange(1))),log10(f_orig(frange(end))),40)'; % I think the logs are mainly for spacing, it is different than the normal linspace
    omega_i = f_i*2*pi; 
    k_i = omega_i.^2./9.81;
    
    % find directions that statisfy DTp-pi/2 < DTp < DTp+pi/2
    [B,maxEIndTheta] = find(Etheta_orig_rotated_sorted == max(Etheta_orig_rotated_sorted(:)),1,'first'); % index of max wave energy
    peakTheta = theta_orig_rotated_sorted(maxEIndTheta); % peak wave direction
    theta_i = linspace(peakTheta - 90, peakTheta + 90, 25); % make a vector of 25 frequencies equally spaced between peak direction +\- pi/2
    theta_i(theta_i>360) = theta_i(theta_i > 360) - 360; % wrap around
    theta_i(theta_i<0) = theta_i(theta_i < 0) + 360; % wrap around
    theta_i = sort(theta_i); % re-sort, [deg]

    %% Interpolate observed spectrum to solution space
    [f_orig_mat,theta_orig_mat] = meshgrid(f_orig,theta_orig_rotated_sorted); % Meshgrid of original frequencies and thetas
    [f_i_mat,theta_i_mat] = meshgrid(f_i,theta_i); % Meshgrid of solution space of frequencies and thetas
    Etheta_i = 10.^griddata(f_orig_mat,theta_orig_mat,log10(Etheta_orig_rotated_sorted'),f_i_mat,theta_i_mat)'; % Interpolate Etheta to new solution space meshgrid (I think the logs are again for spacing?)
    Etheta_i(isnan(Etheta_i)) = 0;
    var_i_init = trapz(f_i',trapz(theta_i,Etheta_i')); % calculate initial variance
    Etheta_i = Etheta_i*var_orig_rotated_sorted/var_i_init; % Re-scaling (rescaling to satisfy Parsevals theorem) bc df and dthetas are different
    var_i = trapz(f_i',trapz(theta_i,Etheta_i')); % calculate final variance to compare to previous spectra/original time series

    %% Calculate amps
    % eta is a Gaussian symmetric about zero, wave heights are Rayleigh
    % distributed, top 1/3 heights of the Rayleigh is 4 * std of eta
    amps_mat = sqrt(Etheta_i.*diff([0 f_i_mat(1,:)]')*mode(diff(theta_i_mat(:,1)))); % m, Ei [m^2/s], bounded integral of (E df dtheta) in the relevant region --> standard deviation of the wave elevation time series
    amps_mat(isnan(amps_mat)) = 0;
    lims_mat = amps_mat / sqrt(2); % sqrt(2) is a free variable based in statistics, standard error on a system with 2 degrees of freedom, most waves would occue within sqrt(2) between Hrms, can mess with

    %% Save variables in structures
    % parameters
    params.Te = Te;
    params.lambdae = lambdae;
    params.ke = ke;
    params.lambdae = lambdae;

    % data
    data.Etheta_orig = Etheta_orig;
    data.theta_orig = theta_orig;
    data.f_orig = f_orig;
    data.var_orig = var_orig;

    data.Etheta_orig_rotated_sorted = Etheta_orig_rotated_sorted;
    data.theta_orig_rotated_sorted = theta_orig_rotated_sorted;
    data.var_orig_rotated_sorted = var_orig_rotated_sorted;

    data.Etheta_i = Etheta_i;
    data.theta_i = theta_i;
    data.f_i = f_i;
    data.omega_i = omega_i;
    data.k_i = k_i;
    data.var_i_init = var_i_init;
    data.var_i = var_i;

    data.f_orig_mat = f_orig_mat;
    data.theta_orig_mat = theta_orig_mat;
    data.f_i_mat = f_i_mat;
    data.theta_i_mat = theta_i_mat;

    data.ampsMat = amps_mat;
    data.limsMat = lims_mat;
end