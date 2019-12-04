% This calculates the Column Integrated Moist Static Energy from the
% statictics file and the imbalance between surface precipitation and evaporation - indicating convection. 

clear
clc
data = sim_data;
c = atm.load_constants('SAM');
vars = ["QV","QN","TABS","p","x","y","z"];
load_var = {'qv','qn','T','x','y','z'};
time = (51:1:100);
time_th = [(1201:24:2400) 2401];
for j = 1:length(data.SSTs)
    for i = 1:length(vars)
        SST = data.SSTs(j);
        exper = ['WTG_SST' num2str(SST) '_96x96x74-1000m-10s'];
        stat_file = [data.output_dir exper '/OUT_STAT/' exper '.nc'];
        P = ncread(stat_file,'PREC');
        E = ncread(stat_file, 'QPEVP');
        LHF = ncread(stat_file, 'LHF');
        P_mean(:,j) = mean(P(1200:2400));
        LHF_mean(:,j) = mean(LHF(1200:2400));
        for k = 1:length(time)
            load_var{i,k} = load_3D(exper,vars{i},time(k));
        end
    end
    rho{j} = ncread(stat_file, 'RHO');
    h1{j} = ncread(stat_file, 'MSE').*c.cp;
    for k = 1:length(time)
        T = load_var{3,k};
        rt = load_var{1,k};
        rt = rt./1000;
        p = load_var{4,k};
        p = p.*100;
        z = load_var{7,k};
        dz = [37;diff(z)];
        p_mat = repmat(p,[1 size(rt,1) size(rt,2)]);
        p_mat = permute(p_mat,[2 3 1]);
        z_mat = repmat(z,[1 size(rt,1) size(rt,2)]);
        z_mat = permute(z_mat,[2 3 1]);
        
        h{j,k} = atm.calculate_MSE(T,p_mat,z_mat,rt);
        h_lat = squeeze(mean(h{j,k},1));
        h_lon = squeeze(mean(h_lat,1));
        h_d{j}(:,k) = h_lon';
        
    end
    mse_in = rho{j}.*h1{j}.*dz;
    for k = 1:length(time_th)-1
        rho_mean{j}(:,k) = mean(rho{j}(:,time_th(k):time_th(k+1)-1),2);
    end
    mse_in_th = rho_mean{j}.*h_d{j}.*dz;
    I = find(z>=1000 & z<=15000);
    mse_ci(j,:) = sum(mse_in(I,:),1);
    for k = 1:length(time_th)-1
        mse_ci_mean(j,k) = mean(mse_ci(j,time_th(k):time_th(k+1)-1),2);
    end
    mse_ci_th(j,:) = sum(mse_in_th(I,:),1);
end
%% Plotting Column Integrated MSE
figure(1)
subplot(2,3,1)
for i = 1:5
    hold on
    plot (time,mse_ci_th(i,:),'Linewidth',2)
end
xlim([min(time) max(time)])
legend('WTG-297K','WTG-299K','WTG-301K','WTG-303K','WTG-305K')
xlabel('Days')
ylabel('MSE')
%view ([90 -90])
title('Column Integrated MSE - Theory (1 K warming)')
hold off
subplot(2,3,2)
for i = 6:10
    %i = j+length(SST1);
    hold on
    plot (time,mse_ci_th(i,:),'Linewidth',2)
end
xlim([min(time) max(time)])
legend('WTG-298K','WTG-300K','WTG-302K','WTG-304K','WTG-306K')
xlabel('Days')
ylabel('MSE')
%view ([90 -90])
title('Column Integrated MSE - Theory (2 K warming)')
hold off
subplot(2,3,3)
for i = 11:15
    %i = j+length(SST1);
    hold on
    plot (time,mse_ci_th(i,:),'Linewidth',2)
end
xlim([min(time) max(time)])
legend('WTG-297.5K','WTG-299.5K','WTG-301.5K','WTG-303.5K','WTG-305.5K')
xlabel('Days')
ylabel('MSE')
%view ([90 -90])
title('Column Integrated MSE - Theory (1.5K warming)')
hold off
subplot(2,3,4)
for i = 1:5
    hold on
    plot(time,mse_ci_mean(i,:),'Linewidth',2)
end
xlim([min(time) max(time)])
legend('WTG-297K','WTG-299K','WTG-301K','WTG-303K','WTG-305K')
xlabel('time (day)')
ylabel('MSE')
%view ([90 -90])
title('Column Integrated MSE (1K warming)')
hold off

subplot(2,3,5)
for i = 6:10
    hold on
    plot(time,mse_ci_mean(i,:),'Linewidth',2)
end
xlim([min(time) max(time)])
legend('WTG-298K','WTG-300K','WTG-302K','WTG-304K','WTG-306K')
xlabel('time (day)')
ylabel('MSE')
%view ([90 -90])
title('Column Integrated MSE (2K warming)')
hold off


subplot(2,3,6)
for i = 6:10
    hold on
    plot(time,mse_ci_mean(i,:),'Linewidth',2)
end
xlim([min(time) max(time)])
legend('WTG-297.5K','WTG-299.5K','WTG-301.5K','WTG-303.5K','WTG-305.5K')
xlabel('time (day)')
ylabel('MSE')
%view ([90 -90])
title('Column Integrated MSE (1.5K warming)')
hold off

saveas(figure(1),[data.output_dir '/Plots/Column_Integrated_MSE.fig'])


%% Plotting Surface Precipitation and Evaporation rates
figure(2)

hold on

plot (data.SSTs(1:5),P_mean(1:5),'-*','Linewidth',2)
plot (data.SSTs(1:5),LHF_mean(1:5)./c.Lv0.*86400,'-o','Linewidth',2)

plot (data.SSTs(6:10),P_mean(6:10),'-*','Linewidth',2)
plot (data.SSTs(6:10),LHF_mean(6:10)./c.Lv0.*86400,'-o','Linewidth',2)

plot (data.SSTs(11:15),P_mean(11:15),'-*','Linewidth',2)
plot (data.SSTs(11:15),LHF_mean(11:15)./c.Lv0.*86400,'-o','Linewidth',2)


%xlim([297 305])
legend('1K warming-P','1K warming-E','2K warming-P','2K warming-E','1.5K warming-P', '1.5K warming-E')
xlabel('Background SST + warming')
ylabel('Surface Precipitation and evaporation')
%view ([90 -90])
title('Surface Precipitation and Evaporation vs Background SST')
hold off
saveas(figure(2), [data.output_dir '/Plots/precip_and_evap.fig'])