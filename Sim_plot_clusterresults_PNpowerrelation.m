clc;
clear;
% close all;
%%
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis')
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis\function')
filedir=('e:\EXP RoDrtest\Exp\sharpened\');
dp=0.002;%particle diameter
Foldersdir_atominfo = (['f:\RoDrtest\RoDr_alphastart0-1\model']);
[num,txt,raw]=xlsread(['C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\parametric_study.xlsx'],1);
keypara=num(4:end,:);
omega=keypara(:,5);
d_p=keypara(:,7);
liquid_content=keypara(:,6);
%% plot cluster size vs variables
% model_var="omega";
% model_var="mu_s";
% model_var="mu_r";
model_var="gamma";
% model_var="liquid_content";
% model_var="dry_mu_r";
% model_var="all";
dbscan="on";
if model_var=="omega"
modelnum=[86:89]';% vary rotation speed
elseif model_var == "mu_s"
modelnum=[90,87,91]';% mu_s
elseif model_var == "mu_r"
modelnum=[92,87,93,94]';% mu_r
elseif model_var == "gamma"
% modelnum=[77,95,96,87,97]';% vary surface tension gamma
% modelnum=[139,124,119]';% omega90
% modelnum=[119,120,134,121:123]';% omega90
% modelnum=[124,125,151,126:128,]';% omega120
modelnum=[139,140,152,141:143]';% omega150
elseif model_var == "liquid_content"
modelnum=[77,87,99:102]';% vary liquid content
elseif model_var == "dry_mu_r"
modelnum=[77,81,82]';
elseif model_var == "all"
% modelnum=[121,126]';
end
%% marker
if modelnum(1,1)==119
    marker="s";
elseif modelnum(1,1)==124
    marker="o";
elseif modelnum(1,1)==139
    marker="^";
else
    marker=".";
end
%% custom colormap
color=colormap(turbo(8));
%%
for modeli=1:size(modelnum,1)
    %% model parameters
    omega=num(modelnum(modeli,1)+3,5);%degree/seconds
    dp=num(modelnum(modeli,1)+3,7);% meter
    mur(modeli,1)=keypara(modelnum(modeli,1),10);
    mus=keypara(modelnum(modeli,1),9);
    surften(modeli,1)=keypara(modelnum(modeli,1),11);
    liquidcontent=keypara(modelnum(modeli,1),6);
    rho=2460;%kg/m3
    WE(modeli,1)=2460*dp/2*(omega/180*pi*0.1)^2/surften(modeli,1);
    %% if not dbscan 
    if dbscan=="off"
%         load([Foldersdir_atominfo num2str(modelnum(modeli,1)) '\model' num2str(modelnum(modeli,1)) '_clustermethod_dotvel_coe1.mat']);
        load([Foldersdir_atominfo num2str(modelnum(modeli,1)) '\model' num2str(modelnum(modeli,1)) '_clustermethod_angle_Angle60_vcorr_conn_pldist_saminter5.mat']);
%         manuel_edges=[1:9,10:10:90,100:100:900,1000:1000:9000];
        manuel_edges=[1:1e4];
        log_edge=log10(manuel_edges);
    
%% extract sizerecord each step and do the fitting for each step
        for stepi=1:size(cluster_binrecord,1)
            cluster_binrecord_estp=cluster_binrecord{stepi,1}';
            cluster_sizerecord_estp=cellfun(@length,cluster_binrecord_estp,'UniformOutput',false);
            cluster_sizerecord_estp=cell2mat(cluster_sizerecord_estp);
            
            %% ------start add single particle in cluster calculation------
            id_clustered_particles=[];
            cluster_count=0;
            for clui=1:size(cluster_binrecord_estp,1)
                id_clustered_particles=[id_clustered_particles;cluster_binrecord_estp{clui,1}'];
                if size(cluster_binrecord_estp{clui,1}',1)>0 && size(cluster_binrecord_estp{clui,1}',1)<10000
                cluster_count=cluster_count+1;
                cluster_id_size(cluster_count,1:2)=[clui size(cluster_binrecord_estp{clui,1}',1)];
                else
                    %continue
                end
            end
            cluster_sizerecord=[cluster_sizerecord;cluster_id_size];

            %% -------end 
            [p_estp,edges_estp]=histcounts(cluster_sizerecord_estp,manuel_edges,"Normalization","Probability");
            [alpha, xmin, L]=plfit(cluster_sizerecord_estp,'xmin',3);
            alpha_estp(stepi,1)=alpha;
            xmin_estp(stepi,1)=xmin;

        end
        [alpha_overall, xmin_overall, L_overall]=plfit(cluster_sizerecord(:,2),'xmin',2);
        alpha_overall_record(modeli,1)=alpha_overall;
        xmin_overall_record(modeli,1)=xmin_overall;
        % % ----------------figure 3--------------------
%         figure(3);hold on;
        power(modeli,1)=-round(mean(alpha_estp),2);
        power_std(modeli,1)=round(std(alpha_estp),2);
%         fig3=plot(effect_edge,p,'.','Color',color(modeli*3,:));%No Shift
%         fig3=plot(effect_edge,effect_edge.^power(modeli,1),'-','Color',color(modeli*3,:));

        %% all flowzone particle
%         flow_zone_single_particles_est=[];
%         flow_zone_single_particles=[];
%         for stepj=1:size(id_centroids_radius_vxyz_record,1)
%             all_particle_tmp=id_centroids_radius_vxyz_record{stepj,1};
%             all_particle_tmp(all_particle_tmp(:,6)>0,:)=[];
%             all_particle_tmp(all_particle_tmp(:,8)>0,:)=[];
%             clustered_particles_id=sort(cell2mat(cluster_binrecord{stepj,1}')');
%             flow_zone_single_particles_est{stepj,1}=all_particle_tmp(~ismember(all_particle_tmp(:,1),clustered_particles_id));
%             flow_zone_single_particles=[flow_zone_single_particles;ones(size(all_particle_tmp,1),1)];
%         end

%         [p,edges]=histcounts([cluster_sizerecord(:,2);flow_zone_single_particles],manuel_edges,"Normalization","Probability");
        [p,edges]=histcounts(cluster_sizerecord(:,2),manuel_edges,"Normalization","Probability");
    else
    %% load strain rate boundary
        boundary_info=load(['f:\RoDrtest\RoDr_alphastart0-1\model' num2str(modelnum(modeli,1)) ...
            '\model' num2str(modelnum(modeli,1)) '_strrbound_strrcoe7_surface.mat']);
        full_bound_filter=[boundary_info.strr_botbound; ...
                           boundary_info.strr_botbound(end,1)+5*dp,boundary_info.strr_botbound(end,2); ...
                           boundary_info.strr_botbound(end,1)+5*dp,boundary_info.strr_botbound(end,2)+10*dp; ...
                           boundary_info.strr_botbound(1,1)-5*dp,boundary_info.strr_botbound(end,2)+10*dp; ...
                           boundary_info.strr_botbound(1,1)-5*dp,boundary_info.strr_botbound(1,2); ...
                           boundary_info.strr_botbound(1,1),boundary_info.strr_botbound(1,2)
                           ];
%         figure(2)
%         plot(boundary_info.strr_botbound(:,1),boundary_info.strr_botbound(:,2),'k-')
%         hold on
%         plot(boundary_info.strr_topbound(:,1),boundary_info.strr_topbound(:,2),'o-')
%         plot(boundary_info.freesurfacepoints(:,1),boundary_info.freesurfacepoints(:,2),'r-')
%         plot(full_bound_filter(:,1),full_bound_filter(:,2),'*-');
%         load([Foldersdir_atominfo num2str(modelnum(modeli,1)) '\model' num2str(modelnum(modeli,1)) '_clustermethod_dotvel_coe1.mat']);
%         load([Foldersdir_atominfo num2str(modelnum(modeli,1)) '\model' num2str(modelnum(modeli,1)) '_clustermethod_angle_Angle60_dbs_minpts2_strrfilt_avefnfilt_pldistcoe0_00000001.mat']);
%         load([Foldersdir_atominfo num2str(modelnum(modeli,1)) '\model' num2str(modelnum(modeli,1)) '_clustermethod_angle_Angle60_dbs_minpts2_strrfilt_avefnfilt.mat']);
load([Foldersdir_atominfo num2str(modelnum(modeli,1)) '\model' num2str(modelnum(modeli,1)) '_clustermethod_angle_Angle60_dbs_minpts2_strrfilt_5surftendpfilt_.mat']);
%         manuel_edges=[1:9,10:10:90,100:100:900,1000:1000:9000];
        manuel_edges=[1:1e4];
        log_edge=log10(manuel_edges(1:end-1));
    %% extract sizerecord each step and do the fitting for each step
        dbs_cluster_size_all=[];
        for stepi=1:size(dbs_cluster_size_record,1)
        %% filter flow zone particles
            id_centroids_radius_vxyz=id_centroids_radius_vxyz_record{stepi,1};
            dbs_cluster_size=dbs_cluster_size_record{stepi,1};
            dbs_cluster_id=dbs_cluster_id_record{stepi,1};
    
            in=inpolygon(id_centroids_radius_vxyz(:,2),id_centroids_radius_vxyz(:,3), ...
                         full_bound_filter(:,1),full_bound_filter(:,2));
            dbs_cluster_size_flowzone=dbs_cluster_size(in,:);
            dbs_cluster_id_flowzone=dbs_cluster_id(in,:);
            id_centroids_flowzone=id_centroids_radius_vxyz(in,1:4);

        %%  unique flowzone cluster size
            [~,idu]=unique(dbs_cluster_id_flowzone);
            dbs_cluster_size_flowzone_unique=dbs_cluster_size_flowzone(idu,:);
            dbs_cluster_size_flowzone_unique=[dbs_cluster_size_flowzone_unique;dbs_cluster_id_flowzone(dbs_cluster_id_flowzone==1)];
    %         dbs_cluster_size_flowzone_unique=[dbs_cluster_size_flowzone_unique];
            dbs_cluster_size_all=[dbs_cluster_size_all;dbs_cluster_size_flowzone_unique];
        %% plot clusters in 3D
%             figure(3)
%             scatter3(id_centroids_radius_vxyz(dbs_cluster_size>1,2), ...
%                 id_centroids_radius_vxyz(dbs_cluster_size>1,3), ...
%                 id_centroids_radius_vxyz(dbs_cluster_size>1,4), ...
%                 [],dbs_cluster_size(dbs_cluster_size>1));
%             % scatter3(centroids_est(:,1), ...
%             %     centroids_est(:,2), ...
%             %     centroids_est(:,3), ...
%             %     [],dbs_cluster_size);
%             set(gca, 'ColorScale', 'log')
%             axis equal;
%             figure(4);hold on;
%             [~,idu]=unique(dbs_cluster_id);
%             dbs_cluster_size=[dbs_cluster_size(idu,:);dbs_cluster_id(dbs_cluster_id==1,:) ];
%             histogram(dbs_cluster_size,[1:1:2000],'normalization','probability');
%             set(gca, 'XScale', 'log')
%             set(gca, 'YScale', 'log')
%             % stop=1;
        %% fitting cluster distribution
            [alpha, xmin, L]=plfit(dbs_cluster_size_flowzone_unique,'range',[1.001:0.001:10.001]);
            alpha_estp(stepi,1)=alpha;
            xmin_estp(stepi,1)=xmin;
        end
    %% all flowzone particle
        [p,edges]=histcounts(dbs_cluster_size_all,manuel_edges,"Normalization","Probability");

    end
    %% organising power fitting parameters
    [alpha_overall, xmin_overall, L_overall]=plfit(dbs_cluster_size_all,'range',[1.001:0.001:10.001]);
    alpha_overall_record(modeli,1)=alpha_overall;
    xmin_overall_record(modeli,1)=xmin_overall;
    power(modeli,1)=-round(mean(alpha_estp),2);
    power_std(modeli,1)=round(std(alpha_estp),2);
    %% mean cluster size
    Nmean(modeli,1)=mean(dbs_cluster_size_all(dbs_cluster_size_all>1));
        %% figure 1 cluster vs variables
        figure(1);hold on;
%         edges(1,1)=2;
        effect_edge=edges(1:end-1);
%         fig1=plot(edges(1:end-1),p*10^(modeli-1),'*','Color',color(modeli*3,:));%Shift
        fig1=plot(effect_edge,p,marker,'Color',color(modeli+1,:));%No Shift
%         fig1=plot(effect_edge,effect_edge.^-mean(alpha_overall_record(modeli,1)),'-','Color',color(modeli*3,:));
        %% figure 2 fitting figure 1
%         log_p=log10(p);
%         figure(2);hold on;
%         fig2=plot(log_edge,log_p,'s','Color',color(modeli+1,:));%No Shift
        % %------------- Set up fittype and options------------
%         ft = fittype( 'a*x+b', 'independent', 'x', 'dependent', 'y' );
%         opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%         opts.Display = 'Off';
%         opts.StartPoint = [0.395894484205228 0.936450101590675];
%         fit_log_p=log_p(log_edge<=3)';%settle range
%         fit_log_p=fit_log_p(~isinf(fit_log_p));
%         fit_log_edge=log_edge(log_edge<=3)';
%         fit_log_edge=fit_log_edge(~isinf(fit_log_p));
%         [fitresult, gof] = fit(fit_log_edge,fit_log_p,ft,opts);
%         fita(modeli,1)=round(fitresult.a,2);
%         fitb(modeli,1)=round(fitresult.b,2);
%         fitgof(modeli,1)=gof;

end

%% figure 1 plot settings
figure(1);
xlabel('Cluster size (\itN)');
ylabel('\itP(N)');
set(gca,'FontSize',18);
set(gca,'FontName','Times New Roman');
if model_var=="omega"
legend(['\omega=60^{\circ}/s; \alpha=' num2str(power(1))], ...
    ['\omega=90^{\circ}/s; \alpha=' num2str(power(2))], ...
    ['\omega=120^{\circ}/s; \alpha=' num2str(power(3))], ...
    ['\omega=150^{\circ}/s; \alpha=' num2str(power(4))],'Location','northeast');
elseif model_var=="mu_s"
legend(['\mu_s=0.1; \alpha=' num2str(power(1))], ...
    ['\mu_s=0.5; \alpha=' num2str(power(2))], ...
    ['\mu_s=0.8; \alpha=' num2str(power(3))], ...
    'Location','northeast');
elseif model_var=="mu_r"
legend(['\mu_r=0.001; \alpha=' num2str(power(1))], ...
    ['\mu_r=0.01; \alpha=' num2str(power(2))], ...
    ['\mu_r=0.1; \alpha=' num2str(power(3))], ...
    ['\mu_r=1; \alpha=' num2str(power(4))], ...
    'Location','northeast', ...
    'FontSize',12);
elseif model_var=="gamma"
    legend( ['\alpha = ' num2str( -power(1)) '\pm' num2str(power_std(1))], ...
    ['\alpha = ' num2str( -power(2)) '\pm' num2str(power_std(2))], ...
    ['\alpha = ' num2str( -power(3)) '\pm' num2str(power_std(3))], ...
    ['\alpha = ' num2str( -power(4)) '\pm' num2str(power_std(4))], ...
    ['\alpha = ' num2str( -power(5)) '\pm' num2str(power_std(5))], ...
    ['\alpha = ' num2str( -power(6)) '\pm' num2str(power_std(6))], ...
    'FontSize',15, ...
    'Location','northeast');
% legend( ['{\itWe} = ' num2str(round(WE(1,1),0)) '; \alpha = ' num2str( -power(1)) '\pm' num2str(power_std(1))], ...
%     ['{\itWe} = ' num2str(round(WE(2,1),2)) '; \alpha = ' num2str( -power(2)) '\pm' num2str(power_std(2))], ...
%     ['{\itWe} = ' num2str(round(WE(3,1),2)) '; \alpha = ' num2str( -power(3)) '\pm' num2str(power_std(3))], ...
%     ['{\itWe} = ' num2str(round(WE(4,1),2)) '; \alpha = ' num2str( -power(4)) '\pm' num2str(power_std(4))], ...
%     ['{\itWe} = ' num2str(round(WE(5,1),2)) '; \alpha = ' num2str( -power(5)) '\pm' num2str(power_std(5))], ...
%     ['{\itWe} = ' num2str(round(WE(6,1),2)) '; \alpha = ' num2str( -power(6)) '\pm' num2str(power_std(6))], ...
%     'FontSize',15, ...
%     'Location','northeast');
% legend( ['dry; \alpha = ' num2str( -power(1)) '\pm' num2str(power_std(1))], ...
%     ['\gamma = 0.0073{\itN/m}; \alpha = ' num2str( -power(2)) '\pm' num2str(power_std(2))], ...
%     ['\gamma = 0.0128{\itN/m}; \alpha = ' num2str( -power(3)) '\pm' num2str(power_std(3))], ...
%     ['\gamma = 0.0365{\itN/m}; \alpha = ' num2str( -power(4)) '\pm' num2str(power_std(4))], ...
%     ['\gamma = 0.073{\itN/m}; \alpha = ' num2str( -power(5)) '\pm' num2str(power_std(5))], ...
%     ['\gamma = 0.146{\itN/m}; \alpha = ' num2str( -power(6)) '\pm' num2str(power_std(6))], ...
%     'FontSize',15, ...
%     'Location','northeast');
% legend( ['dry; \alpha=' num2str( power(1)) '\pm' num2str(power_std(1))], ...
%     '', ...
%     ['\gamma=0.0365N/m; \alpha=' num2str( power(2)) '\pm' num2str(power_std(2))], ...
%     '', ...
%     ['\gamma=0.073N/m; \alpha=' num2str( power(3)) '\pm' num2str(power_std(3))], ...
%     '', ...
%     ['\gamma=0.146N/m; \alpha=' num2str( power(4)) '\pm' num2str(power_std(4))], ...
%     '', ...
%     'FontSize',12, ...
%     'Location','northeast');
elseif model_var=="liquid_content"
legend(['dry; \alpha=' num2str(power(1))], ...
    ['m=0.02; \alpha=' num2str(power(2))], ...
    ['m=0.04; \alpha=' num2str(power(3))], ...
    ['m=0.06; \alpha=' num2str(power(4))], ...
    ['m=0.08; \alpha=' num2str(power(5))], ...
    ['m=0.1; \alpha=' num2str(power(6))], ...
    'Location','northeast');
elseif model_var=="dry_mu_r"
legend(['\mu_r=0.01; \alpha=' num2str(power(1))], ...
    ['\mu_r=0.5; \alpha=' num2str(power(2))], ...
    ['\mu_r=1; \alpha=' num2str(power(3))], ...
    'Location','northeast');
end
box on
xlim([0.8 1e3]);
ylim([1e-7 2]);
% title(["Cluster method: vel fluct vec angle, 5^{\circ}"]);
set(gca,'LineWidth',2);
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
set(gca, 'XTick',[1,10,100,1000,10000]);
set(gca, 'YTick',[0.0000001,0.000001,0.00001,0.0001,0.001,0.01,0.1,1]);
set(gca, 'FontName',"Times New Roman");

%% figure 2 plot settings
% figure(2);
% xlabel('log(N)');
% ylabel('log(P)');
% set(gca,'FontSize',18);
% set(gca,'FontName','Times New Roman');
% if model_var=="omega"
% legend(['\omega=60^{\circ}/s; \alpha=' num2str(power(1))], ...
%     ['\omega=90^{\circ}/s; \alpha=' num2str(power(2))], ...
%     ['\omega=120^{\circ}/s; \alpha=' num2str(power(3))], ...
%     ['\omega=150^{\circ}/s; \alpha=' num2str(power(4))],'Location','northeast');
% elseif model_var=="mu_s"
% legend(['\mu_s=0.1; \alpha=' num2str(power(1))], ...
%     ['\mu_s=0.5; \alpha=' num2str(power(2))], ...
%     ['\mu_s=0.8; \alpha=' num2str(power(3))], ...
%     'Location','northeast');
% elseif model_var=="mu_r"
% legend(['\mu_r=0.001; \alpha=' num2str(power(1))], ...
%     ['\mu_r=0.01; \alpha=' num2str(power(2))], ...
%     ['\mu_r=0.1; \alpha=' num2str(power(3))], ...
%     ['\mu_r=1; \alpha=' num2str(power(4))], ...
%     'Location','northeast');
% elseif model_var=="gamma"
% legend(['dry; \beta=' num2str(fitb(1))], ...
%     ['\gamma=0.0365N/m; \beta=' num2str(fitb(2))], ...
%     ['\gamma=0.073N/m; \beta=' num2str(fitb(3))], ...
%     ['\gamma=0.146N/m; \beta=' num2str(fitb(4))], ...
%     'Location','northeast');
% elseif model_var=="liquid_content"
% legend(['dry; \alpha=' num2str(power(1))], ...
%     ['m=0.02; \alpha=' num2str(power(2))], ...
%     ['m=0.04; \alpha=' num2str(power(3))], ...
%     ['m=0.06; \alpha=' num2str(power(4))], ...
%     ['m=0.08; \alpha=' num2str(power(5))], ...
%     ['m=0.1; \alpha=' num2str(power(6))], ...
%     'Location','northeast');
% elseif model_var=="dry_mu_r"
% legend(['\mu_r=0.01; \alpha=' num2str(power(1))], ...
%     ['\mu_r=0.5; \alpha=' num2str(power(2))], ...
%     ['\mu_r=1; \alpha=' num2str(power(3))], ...
%     'Location','northeast');
% end
% set(gca,'LineWidth',1);
% % set(gca, 'YScale', 'log');
% % set(gca, 'XScale', 'log');
% box on
% % xlim([0.8 1e4]);
% % ylim([1e-6 1e5]);
% title(["Cluster method: vel fluct vec angle, 5^{\circ}"]);
% % title(["Cluster method: velocity correlation","dot vel coe 1"]);

%% figure 2 plot settings
% figure(3);
% xlabel('Cluster size (N particles)');
% ylabel('Probablity');
% set(gca,'FontSize',18);
% set(gca,'FontName','Times New Roman');
% if model_var=="omega"
% legend('\omega=60^{\circ}/s', '\omega=90^{\circ}/s','\omega=120^{\circ}/s','\omega=150^{\circ}/s','Location','northeast');
% elseif model_var=="mu_s"
% legend('\mu_s=0.1', '\mu_s=0.5','\mu_s=0.8','Location','northeast');
% elseif model_var=="mu_r"
% legend('\mu_r=0.001', '\mu_r=0.01','\mu_r=0.1','\mu_r=1','Location','northeast');
% elseif model_var=="gamma"
% legend('dry','\gamma=0.0073N/m', '\gamma=0.0365N/m','\gamma=0.073N/m','\gamma=0.146N/m','Location','northeast');
% elseif model_var=="liquid_content"
% legend('dry','m=0.02', 'm=0.04','m=0.06','m=0.08','m=0.1','Location','northeast');
% end
% set(gca,'LineWidth',1);
% set(gca, 'YScale', 'log');
% set(gca, 'XScale', 'log');
% box on
% xlim([0.8 1e4]);
% ylim([1e-6 1e5]);
% title(["Cluster method: vel fluct vec angle, 5^{\circ}"]);
% % title(["Cluster method: velocity correlation","dot vel coe 1"]);

%% legend for plotting dry and comparing with botamy2002
% set(gca, 'XScale', 'log')
% set(gca, 'YScale', 'log')
% plot([10 10]',[100 110]','k--');
% plot([10 10]',[100 110]','k-.');
% legend( ['\omega = 25rpm; \alpha = ' num2str( -power(1)) '\pm' num2str(power_std(1))], ...
% ['\omega = 20rpm; \alpha = ' num2str( -power(2)) '\pm' num2str(power_std(2))], ...
% ['\omega = 15rpm; \alpha = ' num2str( -power(3)) '\pm' num2str(power_std(3))], ...
% ['\omega \leq 8rpm, quasi-2D; \alpha = ' num2str( -2.9) '\pm' num2str(0.1)], ...
% ['\omega \leq 8rpm, 2D; \alpha = ' num2str(-2) '\pm' num2str(0.1)], ...
% 'FontSize',14, ...
% 'Location','northeast');
% xlim([0.8 10000])
% ylim([1e-7 2]);
% set(gca, 'YTick',[0.0000001,0.000001,0.00001,0.0001,0.001,0.01,0.1,1]);
% set(gca, 'XTick',[1,10,100,1000,10000]);
% box on;