clear; close all; clc; 

%% Directories
addpath(genpath('..\SWIFT-codes'))
exampleNum = 2;
dataFolder = fullfile(pwd,strcat('ExampleData',num2str(exampleNum)));

%% Load data
load(fullfile(dataFolder,"processedData.mat"))

%% Plotting
markerVec = 'x*s^';
colors = cbrewer('qual','Set1',9);
colors_pred = cbrewer('qual','Pastel1',9);
plot_ind = 1;

%% Plot buoy locations
circleThetaVec = 0:0.01:2*pi;
f1 = figure();
plot(0,0,'ko')
hold on
grid on
for jj = 1:params.nbuoys
    xcirc = mean(data.x(:,jj)) + (params.lambdae/4)*cos(circleThetaVec);
    ycirc = mean(data.y(:,jj)) + (params.lambdae/4)*sin(circleThetaVec);
    plot(data.x(:,jj)/params.lambdae,data.y(:,jj)/params.lambdae,'Color',colors(jj,:),'Marker',markerVec(jj))
    plot(xcirc/params.lambdae,ycirc/params.lambdae,'Color',colors(jj,:))
end
plot(data.x_target/params.lambdae,data.y_target/params.lambdae,'v','Color',colors(params.nbuoys+1,:))
axis equal
xlabel('$\frac{x}{\lambda_e}$')
ylabel('$\frac{y}{\lambda_e}$')
title(strcat("$\lambda_{e}$ = ",num2str(round(params.lambdae,2))," m"))

%% Plot time series input
if(params.use_vel)
    numSubPlots = 5;
else
    numSubPlots = 3;
end

f2 = figure();
f2.Position = [488,50.2,560,835.6];
ax2(1) = subplot(numSubPlots,1,1);
hold on
grid on
xline(data.t(data.indWindowStart(plot_ind),1),'Color',colors(params.nbuoys+1,:),'linewidth',1.5)
xline(data.t(data.indWindowEnd(plot_ind),1),'Color',colors(params.nbuoys+1,:),'linewidth',1.5)
xline(data.t(data.predWindowStart(plot_ind),1),'Color',colors(params.nbuoys+2,:),'linewidth',1.5)
xline(data.t(data.predWindowEnd(plot_ind),1),'Color',colors(params.nbuoys+2,:),'linewidth',1.5)
ylabel('$x_{in}$')
title('Time Series Input')

ax2(2) = subplot(numSubPlots,1,2);
hold on
grid on
xline(data.t(data.indWindowStart(plot_ind),1),'Color',colors(params.nbuoys+1,:),'linewidth',1.5)
xline(data.t(data.indWindowEnd(plot_ind),1),'Color',colors(params.nbuoys+1,:),'linewidth',1.5)
xline(data.t(data.predWindowStart(plot_ind),1),'Color',colors(params.nbuoys+2,:),'linewidth',1.5)
xline(data.t(data.predWindowEnd(plot_ind),1),'Color',colors(params.nbuoys+2,:),'linewidth',1.5)
ylabel('$y_{in}$')

ax2(3) = subplot(numSubPlots,1,3);
hold on
grid on
xline(data.t(data.indWindowStart(plot_ind),1),'Color',colors(params.nbuoys+1,:),'linewidth',1.5)
xline(data.t(data.indWindowEnd(plot_ind),1),'Color',colors(params.nbuoys+1,:),'linewidth',1.5)
xline(data.t(data.predWindowStart(plot_ind),1),'Color',colors(params.nbuoys+2,:),'linewidth',1.5)
xline(data.t(data.predWindowEnd(plot_ind),1),'Color',colors(params.nbuoys+2,:),'linewidth',1.5)
ylabel('$z_{in}$')

for ii = 1:params.nbuoys
    plot(ax2(1),data.t,data.x(:,ii),'Color',colors(ii,:))
    plot(ax2(2),data.t,data.y(:,ii),'Color',colors(ii,:))
    plot(ax2(3),data.t,data.z(:,ii),'Color',colors(ii,:))
end

if(params.use_vel)
    ax2(4) = subplot(numSubPlots,1,4);
    hold on
    grid on
    xline(data.t(data.indWindowStart(plot_ind),1),'Color',colors(params.nbuoys+1,:),'linewidth',1.5)
    xline(data.t(data.indWindowEnd(plot_ind),1),'Color',colors(params.nbuoys+1,:),'linewidth',1.5)
    xline(data.t(data.predWindowStart(plot_ind),1),'Color',colors(params.nbuoys+2,:),'linewidth',1.5)
    xline(data.t(data.predWindowEnd(plot_ind),1),'Color',colors(params.nbuoys+2,:),'linewidth',1.5)
    ylabel('$u_{in}$')

    ax2(5) = subplot(5,1,5);
    hold on
    grid on
    xline(data.t(data.indWindowStart(plot_ind),1),'Color',colors(params.nbuoys+1,:),'linewidth',1.5)
    xline(data.t(data.indWindowEnd(plot_ind),1),'Color',colors(params.nbuoys+1,:),'linewidth',1.5)
    xline(data.t(data.predWindowStart(plot_ind),1),'Color',colors(params.nbuoys+2,:),'linewidth',1.5)
    xline(data.t(data.predWindowEnd(plot_ind),1),'Color',colors(params.nbuoys+2,:),'linewidth',1.5)
    xlabel('Time [s]')
    ylabel('$v_{in}$')

    for ii = 1:params.nbuoys
        plot(ax2(4),data.t,data.u(:,ii),'Color',colors(ii,:))
        plot(ax2(5),data.t,data.v(:,ii),'Color',colors(ii,:))
    end
else
    xlabel('Time [s]')
end

linkaxes(ax2,'x')

%%
f3 = figure();
f3.Position = [312.2,99,735.8,658.8];
subplot(2,2,1)
polarPcolor(data.f_orig', data.theta_orig, log10(data.Etheta_orig'));
title('Original Spectra from wavespec.mat')
subtitle(strcat("var = ",num2str(round(data.var_orig,2))," m$^2$"))

subplot(2,2,2)
polarPcolor(data.f_orig', data.theta_orig_rotated_sorted, log10(data.Etheta_orig_rotated_sorted'));
title('Rotated Original Spectra')
subtitle(strcat("var = ",num2str(round(data.var_orig_rotated_sorted,2))," m$^2$"))

subplot(2,2,3)
polarPcolor(data.f_i', data.theta_i, log10(data.Etheta_i'));
title('Interpolated Spectra')
subtitle(strcat("var = ",num2str(round(data.var_i,2))," m$^2$"))


%% Plot reconstruction and prediction
%%% wave elevation
f4 = figure();
f4.Position = [1,49,1536,843.2];
mm4 = tiledlayout(params.nbuoys,3);
title(mm4,'Wave Elevation Reconstruction and Prediction','interpreter','latex')
for ii = 1:params.nbuoys
    % BS
    nexttile
    plot(NaN,NaN,'-','Color',colors(1,:))
    hold on
    grid on
    plot(NaN,NaN,'k-')
    plot(NaN,NaN,'-','Color',[1 1 1]*0.7)

    plot(data.t_in(:,ii),data.z_in(:,ii),'Color',colors(ii,:))
    hold on
    grid on
    plot(data.t_in(:,ii),sol.bs.z_recon(:,ii),'k')

    for jj = 1:length(params.predBuoy)
        if(ii == params.predBuoy(jj))
            plot(data.t_pred_in(:,jj),data.z_pred_in(:,jj),'Color',colors_pred(ii,:))
            plot(data.t_pred_in(:,jj),sol.bs.z_pred(:,jj),'Color',[1 1 1]*0.7)
        end
    end

    xlim([data.t_in(1,1) data.t_pred_in(end,1)])
    ylim([-1 1]*2.2)
    ylabel('$z$ [m]')

    if(ii == 1)
        title('BS')
        legend('Measured','Recreated','Prediction')
    end

    % lims
    nexttile
    plot(NaN,NaN,'-','Color',colors(1,:))
    hold on
    grid on
    plot(NaN,NaN,'k-')
    plot(NaN,NaN,'-','Color',[1 1 1]*0.7)

    plot(data.t_in(:,ii),data.z_in(:,ii),'Color',colors(ii,:))
    hold on
    grid on
    plot(data.t_in(:,ii),sol.lims.z_recon(:,ii),'k')

    for jj = 1:length(params.predBuoy)
        if(ii == params.predBuoy(jj))
            plot(data.t_pred_in(:,jj),data.z_pred_in(:,jj),'Color',colors_pred(ii,:))
            plot(data.t_pred_in(:,jj),sol.lims.z_pred(:,jj),'Color',[1 1 1]*0.7)
        end
    end

    xlim([data.t_in(1,1) data.t_pred_in(end,1)])
    ylabel('$z$ [m]')

    if(ii == 1)
        title('lims')
        legend('Measured','Recreated','Prediction')
    end

    % LSQ
    nexttile
    plot(NaN,NaN,'-','Color',colors(1,:))
    hold on
    grid on
    plot(NaN,NaN,'k-')
    plot(NaN,NaN,'-','Color',[1 1 1]*0.7)

    plot(data.t_in(:,ii),data.z_in(:,ii),'Color',colors(ii,:))
    hold on
    grid on
    plot(data.t_in(:,ii),sol.lsq.z_recon(:,ii),'k')

    for jj = 1:length(params.predBuoy)
        if(ii == params.predBuoy(jj))
            plot(data.t_pred_in(:,jj),data.z_pred_in(:,jj),'Color',colors_pred(ii,:))
            plot(data.t_pred_in(:,jj),sol.lsq.z_pred(:,jj),'Color',[1 1 1]*0.7)
        end
    end

    xlim([data.t_in(1,1) data.t_pred_in(end,1)])
    ylabel('$z$ [m]')

    if(ii == 1)
        title('LSQ')
        legend('Measured','Recreated','Prediction')
    end
end
xlabel('Time [s]')

%%% velocities
if(params.use_vel)
    f5 = figure();
    f5.Position = [1,49,1536,843.2];
    mm5 = tiledlayout(params.nbuoys,3);
    title(mm5,'Horizontal Velocity ($u$) Reconstruction and Prediction','interpreter','latex')
    for ii = 1:params.nbuoys
        % BS
        nexttile
        plot(NaN,NaN,'-','Color',colors(1,:))
        hold on
        grid on
        plot(NaN,NaN,'k-')
        plot(NaN,NaN,'-','Color',[1 1 1]*0.7)
    
        plot(data.t_in(:,ii),data.u_in(:,ii),'Color',colors(ii,:))
        hold on
        grid on
        plot(data.t_in(:,ii),sol.bs.u_recon(:,ii),'k')
    
        for jj = 1:length(params.predBuoy)
            if(ii == params.predBuoy(jj))
                plot(data.t_pred_in(:,jj),data.u_pred_in(:,jj),'Color',colors_pred(ii,:))
                plot(data.t_pred_in(:,jj),sol.bs.u_pred(:,jj),'Color',[1 1 1]*0.7)
            end
        end
    
        xlim([data.t_in(1,1) data.t_pred_in(end,1)])
        ylim([-1 1]*2.2)
        ylabel('$u$ [m/s]')
    
        if(ii == 1)
            title('BS')
            legend('Measured','Recreated','Prediction')
        end
    
        % lims
        nexttile
        plot(NaN,NaN,'-','Color',colors(1,:))
        hold on
        grid on
        plot(NaN,NaN,'k-')
        plot(NaN,NaN,'-','Color',[1 1 1]*0.7)
    
        plot(data.t_in(:,ii),data.u_in(:,ii),'Color',colors(ii,:))
        hold on
        grid on
        plot(data.t_in(:,ii),sol.lims.u_recon(:,ii),'k')
    
        for jj = 1:length(params.predBuoy)
            if(ii == params.predBuoy(jj))
                plot(data.t_pred_in(:,jj),data.u_pred_in(:,jj),'Color',colors_pred(ii,:))
                plot(data.t_pred_in(:,jj),sol.lims.u_pred(:,jj),'Color',[1 1 1]*0.7)
            end
        end
    
        xlim([data.t_in(1,1) data.t_pred_in(end,1)])
        ylabel('$u$ [m/s]')
    
        if(ii == 1)
            title('lims')
            legend('Measured','Recreated','Prediction')
        end
    
        % LSQ
        nexttile
        plot(NaN,NaN,'-','Color',colors(1,:))
        hold on
        grid on
        plot(NaN,NaN,'k-')
        plot(NaN,NaN,'-','Color',[1 1 1]*0.7)
    
        plot(data.t_in(:,ii),data.u_in(:,ii),'Color',colors(ii,:))
        hold on
        grid on
        plot(data.t_in(:,ii),sol.lsq.u_recon(:,ii),'k')
    
        for jj = 1:length(params.predBuoy)
            if(ii == params.predBuoy(jj))
                plot(data.t_pred_in(:,jj),data.u_pred_in(:,jj),'Color',colors_pred(ii,:))
                plot(data.t_pred_in(:,jj),sol.lsq.u_pred(:,jj),'Color',[1 1 1]*0.7)
            end
        end
    
        xlim([data.t_in(1,1) data.t_pred_in(end,1)])
        ylabel('$u$ [m/s]')
    
        if(ii == 1)
            title('LSQ')
            legend('Measured','Recreated','Prediction')
        end
    end
    xlabel('Time [s]')
    
    f6 = figure();
    f6.Position = [1,49,1536,843.2];
    mm6 = tiledlayout(params.nbuoys,3);
    title(mm6,'Horizontal Velocity ($v$) Reconstruction and Prediction','interpreter','latex')
    for ii = 1:params.nbuoys
        % BS
        nexttile
        plot(NaN,NaN,'-','Color',colors(1,:))
        hold on
        grid on
        plot(NaN,NaN,'k-')
        plot(NaN,NaN,'-','Color',[1 1 1]*0.7)
    
        plot(data.t_in(:,ii),data.v_in(:,ii),'Color',colors(ii,:))
        hold on
        grid on
        plot(data.t_in(:,ii),sol.bs.v_recon(:,ii),'k')
    
        for jj = 1:length(params.predBuoy)
            if(ii == params.predBuoy(jj))
                plot(data.t_pred_in(:,jj),data.v_pred_in(:,jj),'Color',colors_pred(ii,:))
                plot(data.t_pred_in(:,jj),sol.bs.v_pred(:,jj),'Color',[1 1 1]*0.7)
            end
        end
    
        xlim([data.t_in(1,1) data.t_pred_in(end,1)])
        ylim([-1 1]*2.2)
        ylabel('$v$ [m/s]')
    
        if(ii == 1)
            title('BS')
            legend('Measured','Recreated','Prediction')
        end
    
        % lims
        nexttile
        plot(NaN,NaN,'-','Color',colors(1,:))
        hold on
        grid on
        plot(NaN,NaN,'k-')
        plot(NaN,NaN,'-','Color',[1 1 1]*0.7)
    
        plot(data.t_in(:,ii),data.v_in(:,ii),'Color',colors(ii,:))
        hold on
        grid on
        plot(data.t_in(:,ii),sol.lims.v_recon(:,ii),'k')
    
        for jj = 1:length(params.predBuoy)
            if(ii == params.predBuoy(jj))
                plot(data.t_pred_in(:,jj),data.v_pred_in(:,jj),'Color',colors_pred(ii,:))
                plot(data.t_pred_in(:,jj),sol.lims.v_pred(:,jj),'Color',[1 1 1]*0.7)
            end
        end
    
        xlim([data.t_in(1,1) data.t_pred_in(end,1)])
        ylabel('$v$ [m/s]')
    
        if(ii == 1)
            title('lims')
            legend('Measured','Recreated','Prediction')
        end
    
        % LSQ
        nexttile
        plot(NaN,NaN,'-','Color',colors(1,:))
        hold on
        grid on
        plot(NaN,NaN,'k-')
        plot(NaN,NaN,'-','Color',[1 1 1]*0.7)
    
        plot(data.t_in(:,ii),data.v_in(:,ii),'Color',colors(ii,:))
        hold on
        grid on
        plot(data.t_in(:,ii),sol.lsq.v_recon(:,ii),'k')
    
        for jj = 1:length(params.predBuoy)
            if(ii == params.predBuoy(jj))
                plot(data.t_pred_in(:,jj),data.v_pred_in(:,jj),'Color',colors_pred(ii,:))
                plot(data.t_pred_in(:,jj),sol.lsq.v_pred(:,jj),'Color',[1 1 1]*0.7)
            end
        end
    
        xlim([data.t_in(1,1) data.t_pred_in(end,1)])
        ylabel('$v$ [m/s]')
    
        if(ii == 1)
            title('LSQ')
            legend('Measured','Recreated','Prediction')
        end
    end
    xlabel('Time [s]')
end

%% amplitudes
f7 = figure();
f7.Position = [132.2,308.6,1222,420];
subplot(1,2,1)
plot(NaN,NaN,'ko')
hold on
grid on
plot(NaN,NaN,'Color',[1 1 1]*0.7)
plot(sol.lsq.A,'ko')
plot(sol.lsq.lims,'Color',[1 1 1]*0.7)
plot(-sol.lsq.lims,'Color',[1 1 1]*0.7)
legend('$A_{lsq}$','lims')
xlabel('Index')
ylabel('Amplitude')
title('Solution Amplitudes and Limits - LSQ')

subplot(1,2,2)
plot(sol.bs.A,'ko')
hold on
grid on
legend('$A_{bs}$')
xlabel('Index')
ylabel('Amplitude')
title('Solution Amplitudes and Limits - BS')

%% error
f8 = figure();
f8.Position = [1,49,1536,843.2];
mm8 = tiledlayout(2,3);
title(mm8,strcat("Error at index ",num2str(plot_ind)),'interpreter','latex')
if(params.use_vel)
    nexttile
    b1 = bar([sol.bs.z_recon_error(plot_ind,:); sol.bs.u_recon_error(plot_ind,:); sol.bs.v_recon_error(plot_ind,:)]);
    for ii = 1:params.nbuoys
        b1(ii).FaceColor = colors(ii,:);
    end
    hold on
    grid on
    legend('buoy 1','buoy 2','buoy 3','buoy 4','location','northwest')
    xticklabels({'z','u','v'})
    ylabel('Normalized error')
    title('Reconstruction - BS')

    nexttile
    b1 = bar([sol.lims.z_recon_error(plot_ind,:); sol.lims.u_recon_error(plot_ind,:); sol.lims.v_recon_error(plot_ind,:)]);
    for ii = 1:params.nbuoys
        b1(ii).FaceColor = colors(ii,:);
    end
    hold on
    grid on
    legend('buoy 1','buoy 2','buoy 3','buoy 4','location','northwest')
    xticklabels({'z','u','v'})
    title('Reconstruction - lims')

    nexttile
    b1 = bar([sol.lsq.z_recon_error(plot_ind,:); sol.lsq.u_recon_error(plot_ind,:); sol.lsq.v_recon_error(plot_ind,:)]);
    for ii = 1:params.nbuoys
        b1(ii).FaceColor = colors(ii,:);
    end
    hold on
    grid on
    legend('buoy 1','buoy 2','buoy 3','buoy 4','location','northwest')
    xticklabels({'z','u','v'})
    title('Reconstruction - LSQ')

    nexttile
    b2 = bar([sol.bs.z_pred_error(plot_ind,:); sol.bs.u_pred_error(plot_ind,:); sol.bs.v_pred_error(plot_ind,:)]);
    for ii = 1:length(params.predBuoy)
        b2(ii).FaceColor = colors_pred(params.predBuoy(ii),:);
    end
    hold on
    grid on
    xticklabels({'z','u','v'})
    ylabel('Normalized error')
    title('Prediction - BS')

    nexttile
    b2 = bar([sol.lims.z_pred_error(plot_ind,:); sol.lims.u_pred_error(plot_ind,:); sol.lims.v_pred_error(plot_ind,:)]);
    for ii = 1:length(params.predBuoy)
        b2(ii).FaceColor = colors_pred(params.predBuoy(ii),:);
    end
    hold on
    grid on
    xticklabels({'z','u','v'})
    title('Prediction - lims')

    nexttile
    b2 = bar([sol.lsq.z_pred_error(plot_ind,:); sol.lsq.u_pred_error(plot_ind,:); sol.lsq.v_pred_error(plot_ind,:)]);
    for ii = 1:length(params.predBuoy)
        b2(ii).FaceColor = colors_pred(params.predBuoy(ii),:);
    end
    hold on
    grid on
    xticklabels({'z','u','v'})
    title('Prediction')
else
    nexttile
    bar([sol.bs.z_recon_error(plot_ind,:)])
    hold on
    grid on
    legend('buoy 1','buoy 2','buoy 3','buoy 4','location','northwest')
    xticklabels({'z'})
    ylabel('Normalized error')
    title('Reconstruction - BS')

    nexttile
    bar([sol.lims.z_recon_error(plot_ind,:)])
    hold on
    grid on
    legend('buoy 1','buoy 2','buoy 3','buoy 4','location','northwest')
    xticklabels({'z'})
    title('Reconstruction - lims')

    nexttile
    bar([sol.lsq.z_recon_error(plot_ind,:)])
    hold on
    grid on
    legend('buoy 1','buoy 2','buoy 3','buoy 4','location','northwest')
    xticklabels({'z'})
    title('Reconstruction - LSQ')

    nexttile
    bar([sol.bs.z_pred_error(plot_ind,:)])
    hold on
    grid on
    xticklabels({'z'})
    ylabel('Normalized error')
    title('Prediction - BS')

    nexttile
    bar([sol.lims.z_pred_error(plot_ind,:)])
    hold on
    grid on
    xticklabels({'z'})
    title('Prediction - lims')

    nexttile
    bar([sol.lsq.z_pred_error(plot_ind,:)])
    hold on
    grid on
    xticklabels({'z'})
    title('Prediction - LSQ')
end

%% Error time series
t_plot_recon = data.t(round((data.indWindowEnd-data.indWindowStart)/2)+data.indWindowStart,1);
t_plot_pred = data.t(round((data.predWindowEnd-data.predWindowStart)/2)+data.predWindowStart,1);

data.t_plot_recon = t_plot_recon;
data.t_plot_pred = t_plot_pred;

% z
f9 = figure;
f9.Position = [1,49,1536,843.2];
mm9 = tiledlayout(2,3);
title(mm9,'Error time series - wave elevation','interpreter','latex')

nexttile
for ii = 1:params.nbuoys
    plot(t_plot_recon,sol.bs.z_recon_error(:,ii),'Color',colors(ii,:))
    hold on
    grid on
end
xlim([data.t(1,1) data.t(end,1)])
ylim([0 1])
ylabel('Normalized error')
title('Reconstruction - BS')

nexttile
for ii = 1:params.nbuoys
    plot(t_plot_recon,sol.lims.z_recon_error(:,ii),'Color',colors(ii,:))
    hold on
    grid on
end
xlim([data.t(1,1) data.t(end,1)])
% ylim([0 1])
ylabel('Normalized error')
title('Reconstruction - lims')

nexttile
for ii = 1:params.nbuoys
    plot(t_plot_recon,sol.lsq.z_recon_error(:,ii),'Color',colors(ii,:))
    hold on
    grid on
end
xlim([data.t(1,1) data.t(end,1)])
ylim([0 1])
title('Reconstruction - LSQ')

nexttile
for ii = 1:params.nbuoys
    for jj = 1:length(params.predBuoy)
        if(ii == params.predBuoy(jj))
            plot(t_plot_pred,sol.bs.z_pred_error(:,jj),'Color',colors_pred(ii,:))
            hold on
            grid on
        end
    end
end
xlim([data.t(1,1) data.t(end,1)])
xlabel('Time [s]')
ylabel('Normalized error')
title('Prediction - BS')

nexttile
for ii = 1:params.nbuoys
    for jj = 1:length(params.predBuoy)
        if(ii == params.predBuoy(jj))
            plot(t_plot_pred,sol.lims.z_pred_error(:,jj),'Color',colors_pred(ii,:))
            hold on
            grid on
        end
    end
end
xlim([data.t(1,1) data.t(end,1)])
xlabel('Time [s]')
title('Prediction - lims')

nexttile
for ii = 1:params.nbuoys
    for jj = 1:length(params.predBuoy)
        if(ii == params.predBuoy(jj))
            plot(t_plot_pred,sol.lsq.z_pred_error(:,jj),'Color',colors_pred(ii,:))
            hold on
            grid on
        end
    end
end
xlim([data.t(1,1) data.t(end,1)])
xlabel('Time [s]')
title('Prediction - LSQ')

%%% velocities
if(params.use_vel)
    % u
    f10 = figure;
    f10.Position = [1,49,1536,843.2];
    mm10 = tiledlayout(2,3);
    title(mm10,'Error time series - horizontal velocity, $u$','interpreter','latex')
    
    nexttile
    for ii = 1:params.nbuoys
        plot(t_plot_recon,sol.bs.u_recon_error(:,ii),'Color',colors(ii,:))
        hold on
        grid on
    end
    xlim([data.t(1,1) data.t(end,1)])
    ylim([0 1])
    ylabel('Normalized error')
    title('Reconstruction - BS')
    
    nexttile
    for ii = 1:params.nbuoys
        plot(t_plot_recon,sol.lims.u_recon_error(:,ii),'Color',colors(ii,:))
        hold on
        grid on
    end
    xlim([data.t(1,1) data.t(end,1)])
    % ylim([0 1])
    ylabel('Normalized error')
    title('Reconstruction - lims')
    
    nexttile
    for ii = 1:params.nbuoys
        plot(t_plot_recon,sol.lsq.u_recon_error(:,ii),'Color',colors(ii,:))
        hold on
        grid on
    end
    xlim([data.t(1,1) data.t(end,1)])
    ylim([0 1])
    title('Reconstruction - LSQ')
    
    nexttile
    for ii = 1:params.nbuoys
        for jj = 1:length(params.predBuoy)
            if(ii == params.predBuoy(jj))
                plot(t_plot_pred,sol.bs.u_pred_error(:,jj),'Color',colors_pred(ii,:))
                hold on
                grid on
            end
        end
    end
    xlim([data.t(1,1) data.t(end,1)])
    xlabel('Time [s]')
    ylabel('Normalized error')
    title('Prediction - BS')
    
    nexttile
    for ii = 1:params.nbuoys
        for jj = 1:length(params.predBuoy)
            if(ii == params.predBuoy(jj))
                plot(t_plot_pred,sol.lims.u_pred_error(:,jj),'Color',colors_pred(ii,:))
                hold on
                grid on
            end
        end
    end
    xlim([data.t(1,1) data.t(end,1)])
    xlabel('Time [s]')
    title('Prediction - lims')
    
    nexttile
    for ii = 1:params.nbuoys
        for jj = 1:length(params.predBuoy)
            if(ii == params.predBuoy(jj))
                plot(t_plot_pred,sol.lsq.u_pred_error(:,jj),'Color',colors_pred(ii,:))
                hold on
                grid on
            end
        end
    end
    xlim([data.t(1,1) data.t(end,1)])
    xlabel('Time [s]')
    title('Prediction - LSQ')
    
    % v
    f9 = figure;
    f9.Position = [1,49,1536,843.2];
    mm9 = tiledlayout(2,3);
    title(mm9,'Error time series - horizontal velocity, $v$','interpreter','latex')
    
    nexttile
    for ii = 1:params.nbuoys
        plot(t_plot_recon,sol.bs.v_recon_error(:,ii),'Color',colors(ii,:))
        hold on
        grid on
    end
    xlim([data.t(1,1) data.t(end,1)])
    ylim([0 1])
    ylabel('Normalized error')
    title('Reconstruction - BS')
    
    nexttile
    for ii = 1:params.nbuoys
        plot(t_plot_recon,sol.lims.v_recon_error(:,ii),'Color',colors(ii,:))
        hold on
        grid on
    end
    xlim([data.t(1,1) data.t(end,1)])
    % ylim([0 1])
    ylabel('Normalized error')
    title('Reconstruction - lims')
    
    nexttile
    for ii = 1:params.nbuoys
        plot(t_plot_recon,sol.lsq.v_recon_error(:,ii),'Color',colors(ii,:))
        hold on
        grid on
    end
    xlim([data.t(1,1) data.t(end,1)])
    ylim([0 1])
    title('Reconstruction - LSQ')
    
    nexttile
    for ii = 1:params.nbuoys
        for jj = 1:length(params.predBuoy)
            if(ii == params.predBuoy(jj))
                plot(t_plot_pred,sol.bs.v_pred_error(:,jj),'Color',colors_pred(ii,:))
                hold on
                grid on
            end
        end
    end
    xlim([data.t(1,1) data.t(end,1)])
    xlabel('Time [s]')
    ylabel('Normalized error')
    title('Prediction - BS')
    
    nexttile
    for ii = 1:params.nbuoys
        for jj = 1:length(params.predBuoy)
            if(ii == params.predBuoy(jj))
                plot(t_plot_pred,sol.lims.v_pred_error(:,jj),'Color',colors_pred(ii,:))
                hold on
                grid on
            end
        end
    end
    xlim([data.t(1,1) data.t(end,1)])
    xlabel('Time [s]')
    title('Prediction - lims')
    
    nexttile
    for ii = 1:params.nbuoys
        for jj = 1:length(params.predBuoy)
            if(ii == params.predBuoy(jj))
                plot(t_plot_pred,sol.lsq.v_pred_error(:,jj),'Color',colors_pred(ii,:))
                hold on
                grid on
            end
        end
    end
    xlim([data.t(1,1) data.t(end,1)])
    xlabel('Time [s]')
    title('Prediction - LSQ')
end
