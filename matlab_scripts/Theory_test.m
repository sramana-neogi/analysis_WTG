%%%%%%% Theory Check %%%%%%
clear
clc
data = sim_data;
c = atm.load_constants('SAM');
vars = ["QV","QN","TABS","p","x","y","z"];
time = (49:1:100);
time_th = [(1153:24:2400) 2401];
t = input('Theory 1: Wwtg = F \nTheory 2: Calc P vs Model P\nTheory 3: P-E budget\n Enter choice of theory: ');
switch t
    case 1 % Theory 1
        n = input('Enter Choice of MSE (1 = theory, 2 = stat_file): ');
        % Wwtg.GMS = F
        
        % Calculate forcing (SHF + LHF)
        
        for j = 1:length(data.SSTs)
            SST = data.SSTs(j);
            exper = ['WTG_SST' num2str(SST) '_96x96x74-1000m-10s'];
            stat_file = [data.output_dir exper '/OUT_STAT/' exper '.nc'];
            LHF(:,j) = ncread(stat_file,'LHF');
            SHF(:,j) = ncread(stat_file,'SHF');
            z(:,j) = ncread(stat_file, 'z');
            dz = [75;diff(z(:,j))];
            rho{j} = ncread(stat_file, 'RHO'); % Density
            Qr(:,j) = ncread(stat_file,'SWNT')-ncread(stat_file,'LWNT')+ncread(stat_file,'LWNS')-ncread(stat_file,'SWNS');
            %Qr_int(:,j) = Qr{j}(1,:)-Qr{j}(74,:);
            F = LHF + SHF + Qr;
            mean_F(j) = mean(F(1200:2400,j),1);
%             for k = 1:length(time_th)-1
%                 mean_F(k,j) = mean(F(time_th(k):time_th(k+1)-1,j),1);
%             end
        end
        
        % Calculate GMS (Neelin and Held '87)
        
        for j = 1:length(data.SSTs)
            SST = data.SSTs(j);
            exper = ['WTG_SST' num2str(SST) '_96x96x74-1000m-10s'];
            stat_file = [data.output_dir exper '/OUT_STAT/' exper '.nc'];
            for i = 1:length(vars)
                for k = 1:length(time)
                    load_var{i,k} = load_3D(exper,vars{i},time(k));
                end
            end
            rho{j} = ncread(stat_file, 'RHO'); % Density
%             mu{j} = ncread(stat_file, 'MCUP'); % Upward mass flux
%             mds{j} = ncread(stat_file, 'MCDNS'); % Downward saturated mass flux
%             mdu{j} = ncread(stat_file, 'MCDNU'); % Downward unsaturated mass flux
%             mass_flux{j} = mu{j}+mds{j}+mdu{j}; % Total mass flux
              mass_flux{j} = ncread(stat_file, 'MFCLD');
            switch n
                case 1 % Calculate MSE
                    for k = 1:length(time)
                        T = load_var{3,k};
                        rt = load_var{1,k};
                        rt = rt./1000;
                        p = load_var{4,k};
                        p = p.*100;
                        z = load_var{7,k};
                        dz = [37;diff(z)];
                        x = load_var{5,k};
                        y = load_var{6,k};
                        p_mat = repmat(p,[1 size(rt,1) size(rt,2)]);
                        p_mat = permute(p_mat,[2 3 1]);
                        z_mat = repmat(z,[1 size(rt,1) size(rt,2)]);
                        z_mat = permute(z_mat,[2 3 1]);
                        
                        h{j,k} = atm.calculate_MSE(T,p_mat,z_mat,rt);
%                         dhdx(1,:,:) = (h{j,k}(2,:,:) - h{j,k}(1,:,:))./1000;
%                         dhdy(:,1,:) = (h{j,k}(:,2,:) - h{j,k}(:,1,:))./1000;
%                         for i = 3:length(x)
%                             dhdx(i-1,:,:) =  (h{j,k}(i,:,:) - h{j,k}(i-2,:,:))./1000;
%                             dhdy(:,i-1,:) =  (h{j,k}(:,i,:) - h{j,k}(:,i-2,:))./1000;
%                         end
%                         dhdx(length(x),:,:) = (h{j,k}(end,:,:) - h{j,k}(end-1,:,:))./1000;
%                         dhdy(:,length(x),:) = (h{j,k}(:,end,:) - h{j,k}(:,end-1,:))./1000;
%                         div_h{j,k} = dhdx + dhdy;
                        h_lat = squeeze(mean(h{j,k},1));
                        h_lon = squeeze(mean(h_lat,1));
                        h_mean{j}(:,k) = h_lon';
                        for l = 1:length(z)-1
                            if (mean(w{j}(l,1200:2400),2)<=0)
                                rdz = 1./(z(l+1)-z(l));
                                l1 = l+1;
                                l2 = l;
                            else
                                rdz = 1./(z(l)-z(l-1));
                                l1 = l;
                                l2 = l-1;
                            end
                            dhdz(l,:) = (h_mean{j}(l1,:)-h_mean{j}(l2,:)).*rdz;
                            %dsdz(l,:) = (s{j}(l1,:)-s{j}(l2,:)).*rdz;
                            dz(l,j) = z(l1)-z(l2);
                        end
                        l1 = length(z);
                        dhdz(l1,:) = (h_mean{j}(l1,:)-h_mean{j}(l1-1,:))./(z(l1)-z(l1-1));          
                        %dsdz(l1,:) = (s{j}(l1,:)-s{j}(l1-1,:))./(z(l1)-z(l1-1));
                        dz(l1,j) = z(l1)-z(l1-1);
                    end
                    for k = 1:length(time_th)-1
                        rho_mean{j}(:,k) = mean(rho{j}(:,time_th(k):time_th(k+1)-1),2);
                        mean_mu{j}(:,k) = mean(mass_flux{j}(:,time_th(k):time_th(k+1)-1),2);
                    end
                    mse_in_mean = rho_mean{j}.*dhdz(:,j).*dz(:,j);
                    mu_in_mean = rho_mean{j}.*mean_mu{j}.*dz;
                    I = find(z>=1000 & z<=15000);
                    
                case 2 % MSE from stat file
                    h{j} = ncread(stat_file, 'MSE').*c.cp;
                    dhdz{j}(1,:) = (h{j}(2,:) - h{j}(1,:))./(z(2,j)-z(1,j));
                    dz = [37;diff(z(:,j))];
                    for i = 3:length(z)
                        dhdz{j}(i-1,:) =  (h{j}(i,:) - h{j}(i-2,:))./(z(i,j)-z(i-2,j));
                    end
                    dhdz{j}(length(z),:) = (h{j}(end,:) - h{j}(end-1,:))./(z(end,j)-z(end-1,j));
                    mse_in = rho{j}.*(-dhdz{j}).*dz;
                    mu_in = rho{j}.*mass_flux{j}.*dz;
                    for k = 1:length(time_th)-1
                        mse_in_mean(:,k) = mean(mse_in(:,time_th(k):time_th(k+1)-1),2);
                        mu_in_mean(:,k) = mean(mu_in(:,time_th(k):time_th(k+1)-1),2);
                    end
                    I = find(z(:,j)>=1000 & z(:,j)<=15000);
            end
            
            num(:,j) = sum(mse_in_mean(I,:),1);
            den(:,j) = sum(mu_in_mean(I,:),1);
            
            GMS(:,j) = num(:,j)./den(:,j);
        end
        
        %% Plotting GMS and Forcing
        
        figure(1)
        subplot(2,3,1)
        for i = 1:5
            hold on
            plot (time,GMS(:,i),'Linewidth',2)
        end
        xlim([min(time) max(time)])
        legend('WTG-297K','WTG-299K','WTG-301K','WTG-303K','WTG-305K')
        xlabel('Days')
        ylabel('GMS')
        %view ([90 -90])
        title('Gross Moist Stability (1 K warming)')
        hold off
        subplot(2,3,2)
        for i = 6:10
            %i = j+length(SST1);
            hold on
            plot (time,GMS(:,i),'Linewidth',2)
        end
        xlim([min(time) max(time)])
        legend('WTG-298K','WTG-300K','WTG-302K','WTG-304K','WTG-306K')
        xlabel('Days')
        ylabel('GMS')
        %view ([90 -90])
        title('Gross Moist Stability (2 K warming)')
        hold off
        subplot(2,3,3)
        for i = 11:15
            %i = j+length(SST1);
            hold on
            plot (time,GMS(:,i),'Linewidth',2)
        end
        xlim([min(time) max(time)])
        legend('WTG-297.5K','WTG-299.5K','WTG-301.5K','WTG-303.5K','WTG-305.5K')
        xlabel('Days')
        ylabel('GMS')
        %view ([90 -90])
        title('Gross Moist Stability (1.5K warming)')
        hold off
        subplot(2,3,4)
        for i = 1:5
            hold on
            plot(data.SSTs(i),mean_F(:,i),'*')
        end
        %xlim([min(time) max(time)])
        %legend('WTG-297K','WTG-299K','WTG-301K','WTG-303K','WTG-305K')
        xlabel('time (day)')
        ylabel('Forcing')
        %view ([90 -90])
        title('Background warming as forcing (1K warming)')
        hold off
        
        subplot(2,3,5)
        for i = 6:10
            hold on
            plot(data.SSTs(i),mean_F(:,i),'*')
        end
        %xlim([min(time) max(time)])
        %legend('WTG-298K','WTG-300K','WTG-302K','WTG-304K','WTG-306K')
        xlabel('time (day)')
        ylabel('Forcing')
        %view ([90 -90])
        title('Background warming as forcing (2K warming)')
        hold off
        
        
        subplot(2,3,6)
        for i = 6:10
            hold on
            plot(data.SSTs(i),mean_F(:,i),'*')
        end
        %xlim([min(time) max(time)])
        %legend('WTG-297.5K','WTG-299.5K','WTG-301.5K','WTG-303.5K','WTG-305.5K')
        xlabel('time (day)')
        ylabel('Forcing')
        %view ([90 -90])
        title('Background warming as forcing (1.5K warming)')
        hold off
        
%         if (n == 1)
%             saveas(figure(1),[data.output_dir '/Plots/Theory_test_calculated_MSE.fig'])
%         elseif (n==2)
%             saveas(figure(1),[data.output_dir '/Plots/Theory_test_stat_MSE.fig'])
%         end
        
    case 2
        LHF = zeros(2400,15);
        SHF = zeros(2400,15);
        Qr = cell(1,15);
        z = zeros(74,15);
        rho = cell(1,15);
       % Qr_int = zeros(2400,15);
        h = cell(1,15);
        w = cell(1,15);
        for j = 1:length(data.SSTs)
            SST = data.SSTs(j);
            exper = ['WTG_SST' num2str(SST) '_96x96x74-1000m-10s'];
            stat_file = [data.output_dir exper '/OUT_STAT/' exper '.nc'];
            LHF(:,j) = ncread(stat_file,'LHF');
            SHF(:,j) = ncread(stat_file,'SHF');
            Qr{j} = ncread(stat_file,'SWNT')-ncread(stat_file,'LWNT')+ncread(stat_file,'LWNS')-ncread(stat_file,'SWNS');
            z(:,j) = ncread(stat_file, 'z');
            rho{j} = ncread(stat_file, 'RHO'); % Density
           
            h{j} = ncread(stat_file, 'MSE').*c.cp;
            s{j} = ncread(stat_file, 'DSE').*c.cp;
            w{j} = ncread(stat_file, 'WOBS');
            P(:,j) = ncread(stat_file,'PREC');
            P = c.Lv0.*P./86400;
            P_mean(:,j) = mean(P(1200:2400,j),1);
            for k = 1:length(z(:,j))-1
                if (mean(w{j}(k,1200:2400),2)<=0)
                    rdz = 1./(z(k+1,j)-z(k,j));
                    k1 = k+1;
                    k2 = k;
                else
                    rdz = 1./(z(k,j)-z(k-1,j));
                    k1 = k;
                    k2 = k-1;
                end
                dhdz(k,:) = (h{j}(k1,:)-h{j}(k2,:)).*rdz;
                dsdz(k,:) = (s{j}(k1,:)-s{j}(k2,:)).*rdz;
                dz(k,j) = z(k1,j)-z(k2,j);
            end
            k1 = length(z(:,j));
            dhdz(k1,:) = (h{j}(k1,:)-h{j}(k1-1,:))./(z(k1,j)-z(k1-1,j));          
            dsdz(k1,:) = (s{j}(k1,:)-s{j}(k1-1,:))./(z(k1,j)-z(k1-1,j));
            dz(k1,j) = z(k1,j)-z(k1-1,j);
            num = sum(mean(rho{j}(:,1200:2400),2).*mean(w{j}(:,1200:2400).*dhdz(:,1200:2400),2).*dz(:,j),1);
            den = sum(mean(rho{j}(:,1200:2400),2).*mean(w{j}(:,1200:2400).*dsdz(:,1200:2400),2).*dz(:,j),1);
            ngms(j) = num/den;
             % Vertical integral of Qr
            Qr_int(j) = mean(Qr{j}(1200:2400),1);
            F(j) = mean(LHF(1200:2400,j)) + mean(SHF(1200:2400,j)) + Qr_int(j);
            P_cal(:,j) = F(j)./ngms(j) - Qr_int(j) - mean(SHF(1200:2400,j));
            budget_1(:,j) = (mean(LHF(1200:2400,j))+mean(SHF(1200:2400,j))+Qr_int(j))./(mean(SHF(1200:2400,j))+mean(P(1200:2400,j))+Qr_int(j));
            budget_2(:,j) = (P_mean(j)-mean(LHF(1200:2400,j),1))./(den - num);
        end
end