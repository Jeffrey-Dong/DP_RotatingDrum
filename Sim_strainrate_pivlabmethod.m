clc;
clear;
%%
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\dataanalysis\');
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis\function');
[num,txt,raw]=xlsread(['C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\parametric_study.xlsx'],1);
keypara=num(4:end,:);
density_threshold=0.6;
%% select model and parameters
% model_var="omega";
% model_var="mu_s";
% model_var="mu_r";
model_var='gamma';
% model_var='liquid_content'
% modrel_var='all';
if model_var=="omega"
modelnum=[86:89]';% vary rotation speed
elseif model_var == "mu_s"
modelnum=[90,87,91]';% mu_s
elseif model_var == "mu_r"
modelnum=[92,87,93,94]';% mu_r
elseif model_var == "gamma"
% modelnum=[103,95,96,87,97]';% m_r=0.01
% modelnum=[112:114,92]';% m_r=0.001
modelnum=[119,120,134,121:123]';% w90
% modelnum=[124,125,151,126:128]';% w120
% modelnum=[139,140,152,141:143]';% w150
% modelnum=[120,134,125,140,151,152];
if modelnum(1,1)==119
    marker="s";
elseif modelnum(1,1)==124
    marker="o";
elseif modelnum(1,1)==139
    marker="^";
else
    marker=".";
end
elseif model_var == "liquid_content"
modelnum=[103,87,99:102]';% vary liquid content
else
% modelnum=[112,103,109:111]';%dry mu_r
modelnum=[110:114]';%dry mu_r
end
% modelnum=[86]';
% omega_sample=60;
strainrate_coef=7;
%% custom colormap
color=colormap(turbo(8));
%% running the algorithm
for modeli=1:size(modelnum)

    datafile=(['f:\RoDrtest\RoDr_alphastart0-1\model' num2str(modelnum(modeli,1))]);
    load([datafile '\model' num2str(modelnum(modeli,1)) '_velmatrix_frames.mat']);
 
%% model parameters
omega=num(modelnum(modeli,1)+3,5);%degree/seconds
dp=num(modelnum(modeli,1)+3,7);% meter
mur(modeli,1)=keypara(modelnum(modeli,1),10);
mus=keypara(modelnum(modeli,1),9);
surften(modeli,1)=keypara(modelnum(modeli,1),11);
liquidcontent=keypara(modelnum(modeli,1),6);
rho=2460;%kg/m3
WE(modeli,1)=rho*dp/2*(omega/180*pi*0.1)^2/surften(modeli,1);
vave_max=sqrt(max(meanvelxincell,[],"all")^2+max(meanvelzincell,[],"all")^2);
fluid_viscousity=0.00089;%Pas
Ca(modeli,1)=vave_max*fluid_viscousity/surften(modeli,1);
%% ========= plot strain rate contour: accumulate u&v then obtain strain rate ===========    
    strrave2=strain(Xpf',Zpf',-meanvelxincell,-meanvelzincell);
    strrave2_tmp=strrave2;
% % ------------strain rate contour plot----------------
%     figure(1)
%     for tmpstri=1:size(strrave2_tmp,1)
%         for tmpstrj=1:size(strrave2_tmp,2)
%             if strrave2_tmp(tmpstri,tmpstrj)<=-10 || strrave2_tmp(tmpstri,tmpstrj)==0
%                 strrave2_tmp(tmpstri,tmpstrj)=NaN;
%             else
%             end
%         end
%     end
%     [c2,h2]=contourf(Xpf,Zpf,strrave2_tmp,200);
%     set(h2, 'edgecolor','none');
%     axis equal;
%     cbar=colorbar;
%     clim([0 40]);
%     title(["model" num2str(modelnum(modeli,1))]);
% %     cbar_title=('$\dot{char(949)}$','Interpreter','latex');
%     title(cbar,'$\dot{\varepsilon}$','Interpreter','latex', ...
%         'FontSize',18, ...
%         'FontName','Times New Roman');
%     cbar.Position=[0.83,0.11,0.038,0.75];

%%
strrave2_tmp1=strrave2_tmp;
strrave2_tmp1(isnan(strrave2_tmp1))=0;
% figure(5);
% [c2,h2]=contourf(Xpf,Zpf,strrave2_tmp1,200);
%     set(h2, 'edgecolor','none');
%     axis equal;
%     colorbar;

strr_thresh=omega/180*pi*strainrate_coef;
[strr_bound_tmp]=contours(Xpf,Zpf,strrave2_tmp1,[strr_thresh strr_thresh]);
strr_bound_tmp=strr_bound_tmp';
% strr_bound_tmp(strr_bound_tmp(:,1)>=1,:)=[];%remove index array
circle_index=find(strr_bound_tmp(:,2)==max(strr_bound_tmp(:,2)));
strr_bound_tmp1=strr_bound_tmp(circle_index+1:circle_index+max(strr_bound_tmp(:,2)),:);
% figure(4)
% hold on;
% plot(strr_bound_tmp1(:,1),strr_bound_tmp1(:,2),'s','Color',color(modeli+1,:));

%% method 2 extract shear region lower boundary: clusterdata function
% dp=0.003;%particle size
xws=2*dp;%x window size
steps=1*dp;%moving step size
[strr_botbound,strr_topbound]=strainrate_bound(strr_bound_tmp1,xws,steps);
figure(1);hold on;
plot(strr_botbound(:,1),strr_botbound(:,2),'k.');
plot(strr_topbound(:,1),strr_topbound(:,2),'k*');
% datafilename=[char(expnum) 'strinrate_boundary.mat'];
% save(datafilename,'boundary_in_pixel');

%% fitting and find the strainrate depth
datafilename=[datafile '\model' num2str(modelnum(modeli,1)) '_up_surface_points.mat'];
load(datafilename)
% load([datafile '\model' num2str(modelnum(modeli,1)) '_velmatrix_frames.mat']);
accu_up_surf_points=[];
for i=1:size(up_surf_points,1)
    accu_up_surf_points=[accu_up_surf_points;up_surf_points{i,1}];
end
%% ------------------- surface std shade ----------------
[xfit_range2,inBetween,xfit_range_extracted,surf_mid]=errorshade(accu_up_surf_points);% xfit_range_extracted is the surfacex and surf_mid is the surfacey
freesurfacepoints=[xfit_range_extracted,surf_mid];
marker_color_type='ko';
xwsfit=10*dp;%x window size
stepsfit=1*dp;%moving step size
[strr_boundary_depth,flow_depth,flow_depth_normbydp,maxflowdepth_norm,swfit,fitangle_a,maxflowdepth_std]=strr_depth_fun(strr_botbound,strr_topbound,xfit_range_extracted,surf_mid,dp,xwsfit,stepsfit);
maxflowdepth_norm_record(modeli,1)=round(maxflowdepth_norm,1);
maxflowdepth_std_record(modeli,1)=round(maxflowdepth_std/dp,1);
%% ------------------- strain area ----------------
poly=polyshape([strr_botbound;flip(strr_topbound)]);
strr_area(modeli,1)=area(poly);
%% save strain rate data
% cd(datafile);
% datafilename=(['model' num2str(modelnum(modeli,1)) '_strainrate_depth.mat']);
% save(datafilename,'flow_depth_normbydp');
% save(datafilename,'maxflowdepth_norm','-append');
% save(datafilename,'swfit','-append');

% % ----------------save boundary data to mat file-------------   
% atominfo_datafilename=['E:\RoDrtest\RoDr_alphastart0-1\model' num2str(modelnum(modeli,1)) '\model' num2str(modelnum(modeli,1)) '_strrbound_strrcoe' num2str(strainrate_coef) '_surface.mat'];   
% save(atominfo_datafilename,'strr_botbound','-v7.3');
% save(atominfo_datafilename,'strr_topbound','-append');
% save(atominfo_datafilename,'freesurfacepoints','-append');

%% -----------figure (1) plot strain rate bound---------------------
if model_var=="all"
%     color=colormap(turbo(36));
    figure(1);hold on;
    plot(xfit_range_extracted,surf_mid, ...
        "r--","LineWidth",2);
    plot(strr_botbound(:,1),strr_botbound(:,2), ...
        "r--","LineWidth",2);
    title(model_var);


%     cd('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Manuscript figures\Sim_strainrate_boundary');
%     saveas(gca,['model' num2str(modelnum(modeli)) '.tif']);
%     close figure 1;
else
    figure(1);hold on;
    plot([xfit_range_extracted(1:end-1);strr_botbound(:,1)],[surf_mid(1:end-1);strr_botbound(:,2)], ...
        ".","LineWidth",2, ...
        'Color',color(modeli+1,:));
    title(model_var);
%     cd('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Manuscript figures\Sim_strainrate_boundary');
%     saveas(gca,['model' num2str(modelnum(modeli)) '.tif']);
%     close figure 1;
end

%% figure (2) WE vs max flow depth
figure(2); hold on;
if model_var=="omega"
   fig2=plot(omega,maxflowdepth_norm,'ko');
%     xlim([-0.02 1.1]);
    xlabel("\omega","fontsize",20);
    ylabel("Normalized flow depth","fontsize",15);
elseif model_var == "mu_s"
   fig2=plot(mus,maxflowdepth_norm,'ko');
%       xlim([-0.02 1.1]);
    xlabel("\mu_s","fontsize",20);
    ylabel("Normalized flow depth","fontsize",15);
elseif model_var == "mu_r"
   fig2=plot(mur(modeli,1),maxflowdepth_norm,'ko');
   xlim([-0.02 1.1]);
    xlabel("\mu_r","fontsize",20);
    ylabel("Normalized flow depth","fontsize",15);
elseif model_var == "gamma"
%    fig2=plot(surften(modeli,1),maxflowdepth_norm,'k*');
% if WE(modeli,1) == inf
%     WE(modeli,1)=100;
% else
% end
   fig2=plot(WE(modeli,1),maxflowdepth_norm,marker,'MarkerFaceColor','k','MarkerEdgeColor','none');
%    xlim([-0.02 0.16]);
    xlabel("\gamma","fontsize",20);
    ylabel("Normalized flow depth","fontsize",15);
elseif model_var == "liquid_content"
   fig2=plot(liquidcontent,maxflowdepth_norm,'ko');
else
   fig2=plot(1/WE(modeli,1),maxflowdepth_norm,'ko');
end
box on;
set(fig2,"MarkerSize",10);

figure(2)
xlabel("\itWe");
ylabel("{\itD}_{max}/\itd");
set(gca, 'XScale', 'log');
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',15);
xlim([0.1 200]);
ylim([11 20]);

%% figure settings
figure(2);
set(gca, 'XScale', 'log')
xlim([0.05 300])
ylim([11 20])
set(gca,'XTick',[0.1,1,10,100]);
legend('\omega=15rpm', ...
'', ...
'', ...
'', ...
'', ...
'', ...
'\omega=20rpm', ...
'', ...
'', ...
'', ...
'', ...
'', ...
'\omega=25rpm', ...
'', ...
'', ...
'', ...
'', ...
'', ...
'Location','northwest', ...
'fontsize',12);

%% -----------figure (3) plot flow depth along x-axis---------------------
if model_var=="all"
    figure(3);hold on;
    plot(swfit(fitangle_a>0),flow_depth(fitangle_a>0)/dp, ...
        "r--","LineWidth",2, ...
        'Color',color(modeli+1,:));
    title(model_var);

%     cd('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Manuscript figures\Sim_strainrate_boundary');
%     saveas(gca,['model' num2str(modelnum(modeli)) '.tif']);
%     close figure 1;
else
    figure(3);hold on;
    plot(swfit(find(fitangle_a>0,1)+1:end,1)/dp,flow_depth(find(fitangle_a>0,1)+1:end,1)/dp, ...
        [char(marker) char("-")],"LineWidth",1, ...
        'Color',color(modeli+1,:), ...
        'MarkerIndices',[1:3:size(swfit(fitangle_a>0),1)], ...
        'MarkerSize',6);
    title(model_var);
%     cd('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Manuscript figures\Sim_strainrate_boundary');
%     saveas(gca,['model' num2str(modelnum(modeli)) '.tif']);
%     close figure 1;
end
%%
figure(3)
ylabel("\ith/d_p","fontsize",15);
xlabel("\itx/d_p")
title '';
xlim([-0.06/dp 0.1/dp]);
ylim([2 18]);
set(gca,'XTick',[-0.06 -0.04 -0.02 0 0.02 0.04 0.06 0.08 0.1]/dp);
set(gca,'YTick',[2:2:18]);
set(gca, 'FontSize',18);
set(gca, 'FontName','Times New Roman');
set(gca,"LineWidth",2);
box on;
%% modify lines widthes for insets
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 3;
  lines(i).MarkerSize = 14.0;
end
set(gca,'FontSize',30);
% set(gca,'XTick',[15,20,25,30])
% set(gca,'YTick',[11,12,13,14])
% xlabel('')
% ylabel('')

%% -----------figure (4) plot strain rate depth along x-axis---------------------
if model_var=="all"
    figure(4);hold on;
    plot(swfit(fitangle_a>0),strr_boundary_depth(fitangle_a>0)/dp, ...
        "r--","LineWidth",2, ...
        'Color',color(modeli+1,:));
    title(model_var);

%     cd('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Manuscript figures\Sim_strainrate_boundary');
%     saveas(gca,['model' num2str(modelnum(modeli)) '.tif']);
%     close figure 1;
else
    figure(4);hold on;
    plot(swfit(fitangle_a>0)/dp,strr_boundary_depth(fitangle_a>0)/dp, ...
        [char(marker) char("-")],"LineWidth",1, ...
        'Color',color(modeli+1,:), ...
        'MarkerIndices',[1:3:size(swfit(fitangle_a>0),1)]);
    title(model_var);
%     cd('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Manuscript figures\Sim_strainrate_boundary');
%     saveas(gca,['model' num2str(modelnum(modeli)) '.tif']);
%     close figure 1;
end
%%
figure(4)
ylabel("{\itD}_{strr}/\itd","fontsize",15);
xlabel("\itx/d")
title '';
xlim([-0.06/dp 0.1/dp]);
ylim([6 18]);
set(gca,'XTick',[-0.06 -0.04 -0.02 0 0.02 0.04 0.06 0.08 0.1]/dp);
set(gca, 'FontSize',15);
set(gca, 'FontName','Times New Roman');
box on;
%% -----------figure (5) plot plug depth along x-axis---------------------
if model_var=="all"
    figure(5);hold on;
    plot(swfit(fitangle_a>0),(flow_depth(fitangle_a>0)-strr_boundary_depth(fitangle_a>0))/dp, ...
        "r--","LineWidth",2, ...
        'Color',color(modeli+1,:));
    title(model_var);

%     cd('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Manuscript figures\Sim_strainrate_boundary');
%     saveas(gca,['model' num2str(modelnum(modeli)) '.tif']);
%     close figure 1;
else
    figure(5);hold on;
    plot(swfit(fitangle_a>0)/dp,(flow_depth(fitangle_a>0)-strr_boundary_depth(fitangle_a>0))/dp, ...
        [char(marker) char("-")],"LineWidth",1, ...
        'Color',color(modeli+1,:), ...
        'MarkerIndices',[1:3:size(swfit(fitangle_a>0),1)]);
    title(model_var);
%     cd('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Manuscript figures\Sim_strainrate_boundary');
%     saveas(gca,['model' num2str(modelnum(modeli)) '.tif']);
%     close figure 1;
end
%%
figure(5)
ylabel("{\itD}_{plug}/\itd","fontsize",15);
xlabel("\itx/d")
title '';
xlim([-0.06/dp 0.1/dp]);
ylim([6 18]);
set(gca,'XTick',[-0.06 -0.04 -0.02 0 0.02 0.04 0.06 0.08 0.1]/dp);
set(gca, 'FontSize',15);
set(gca, 'FontName','Times New Roman');
box on;
%% figure (6) WE vs max flow depth
figure(6); hold on;
if model_var=="omega"
   fig2=plot(omega,maxflowdepth_norm,'ko');
%     xlim([-0.02 1.1]);
    xlabel("\omega","fontsize",20);
    ylabel("Normalized flow depth","fontsize",15);
elseif model_var == "mu_s"
   fig2=plot(mus,maxflowdepth_norm,'ko');
%       xlim([-0.02 1.1]);
    xlabel("\mu_s","fontsize",20);
    ylabel("Normalized flow depth","fontsize",15);
elseif model_var == "mu_r"
   fig2=plot(mur(modeli,1),maxflowdepth_norm,'ko');
   xlim([-0.02 1.1]);
    xlabel("\mu_r","fontsize",20);
    ylabel("Normalized flow depth","fontsize",15);
elseif model_var == "gamma"
%    fig2=plot(surften(modeli,1),maxflowdepth_norm,'k*');
% if WE(modeli,1) == inf
%     WE(modeli,1)=100;
% else
% end
   fig2=plot(WE(modeli,1),strr_area(modeli,1)*0.03/(4/3*pi*(dp/2)^3),marker,'MarkerFaceColor','k','MarkerEdgeColor','none');
%    xlim([-0.02 0.16]);
elseif model_var == "liquid_content"
   fig2=plot(liquidcontent,maxflowdepth_norm,'ko');
else
   fig2=plot(1/WE(modeli,1),maxflowdepth_norm,'ko');
end
box on;
set(fig2,"MarkerSize",10);

figure(6)
xlabel("\itWe");
ylabel("{\itV}_{strr}/\itV_p");
set(gca, 'XScale', 'log');
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',15);
xlim([0.05 300]);
% ylim([0.002 0.005])

%% figure settings
figure(6);
set(gca, 'XScale', 'log')

set(gca,'XTick',[0.1,1,10,100]);
legend('\omega=15rpm', ...
'', ...
'', ...
'', ...
'\omega=20rpm', ...
'', ...
'', ...
'', ...
'\omega=25rpm', ...
'', ...
'', ...
'', ...
'Location','northwest', ...
'fontsize',12);

end
%% legend
figure(1)
if model_var=="omega"
    legend(['\omega=60^{\circ}/s; \delta=' num2str(maxflowdepth_norm_record(1))], ...
        ['\omega=90^{\circ}/s; \delta=' num2str(maxflowdepth_norm_record(2))], ...
        ['\omega=120^{\circ}/s; \delta=' num2str(maxflowdepth_norm_record(3))], ...
        ['\omega=150^{\circ}/s; \delta=' num2str(maxflowdepth_norm_record(4))], ...
        'Location','southeast', ...
        'fontsize',12);
    elseif model_var=="mu_s"
    legend(['\mu_s=0.1; \delta=' num2str(maxflowdepth_norm_record(1))], ...
        ['\mu_s=0.5; \delta=' num2str(maxflowdepth_norm_record(2))], ...
        ['\mu_s=0.8; \delta=' num2str(maxflowdepth_norm_record(3))], ...
        'Location','southeast', ...
        'fontsize',12);
    elseif model_var=="mu_r"
    ylim([-0.1 0.1]);
    legend(['\mu_r=0.001; \delta=' num2str(maxflowdepth_norm_record(1))], ...
        ['\mu_r=0.01; \delta=' num2str(maxflowdepth_norm_record(2))], ...
        ['\mu_r=0.1; \delta=' num2str(maxflowdepth_norm_record(3))], ...
        ['\mu_r=1; \delta=' num2str(maxflowdepth_norm_record(4))], ...
        'Location','southeast', ...
        'fontsize',12);
    elseif model_var=="gamma"
    legend(['dry; \delta=' num2str(maxflowdepth_norm_record(1))], ...
        ['\gamma=0.0365N/m; \delta=' num2str(maxflowdepth_norm_record(2))], ...
        ['\gamma=0.073N/m; \delta=' num2str(maxflowdepth_norm_record(3))], ...
        ['\gamma=0.146N/m; \delta=' num2str(maxflowdepth_norm_record(4))], ...
        'Location','southeast', ...
        'fontsize',12);
    elseif model_var=="liquid_content"
    legend(['dry; \delta=' num2str(maxflowdepth_norm_record(1))], ...
        ['m=0.02; \delta=' num2str(maxflowdepth_norm_record(2))], ...
        ['m=0.04; \delta=' num2str(maxflowdepth_norm_record(3))], ...
        ['m=0.06; \delta=' num2str(maxflowdepth_norm_record(4))], ...
        ['m=0.08; \delta=' num2str(maxflowdepth_norm_record(5))], ...
        ['m=0.1; \delta=' num2str(maxflowdepth_norm_record(6))], ...
        'Location','southeast', ...
        'fontsize',12);
    end
    set(gca,'FontName','Times New Roman');



%% =======================Functions===========================================================================
% =======================Functions===========================================================================
function out=strain(x,y,vx,vy)
hx = x(1,:);
hy = y(:,1);
[pvxx, pvxy] = gradient(vx, hx, hy);
[pvyx, pvyy] = gradient(vy, hx, hy); %#ok<*ASGLU>
% out = (-pvxy+pvyx);
out = (pvxx-pvyy);
end

function out=shear(x,y,u,v)
hx = x(1,:);
hy = y(:,1);
[junk, py] = gradient(u, hx, hy);
[qx, junk] = gradient(v, hx, hy);
out= 0.5*(qx+py);
end

function [strr_bound,strr_top]=strainrate_bound(strr_bound_points,xws,steps)
xrange=[min(strr_bound_points(:,1)) max(strr_bound_points(:,1))];%fron extreme left point to right point
sw=[xrange(1):steps:xrange(2)]';%sample window, set steps to xws to extract boundary points (otherwise overlap points)
strr_bound=[];
strr_top=[];
    for wi=1:size(sw,1)
        window=[sw(wi)-xws/2 sw(wi)+xws/2];
        points_index1=find(strr_bound_points(:,1)>window(1));
        points_index2=find(strr_bound_points(:,1)<window(2));
        points_index=points_index1(ismember(points_index1,points_index2));
        points_tmp=strr_bound_points(points_index,:);
        points_tmp_cell{wi,1}=points_tmp;
        cluster_index=clusterdata(points_tmp,'Maxclust',3);
        ycluster1=mean(points_tmp(cluster_index==1,2));
        ycluster2=mean(points_tmp(cluster_index==2,2));
        ycluster3=mean(points_tmp(cluster_index==3,2));
        ycluster4=mean(points_tmp(cluster_index==4,2));
        ycluster=[ycluster1;ycluster2;ycluster3;ycluster4];
        ycluster_minindex=find(ycluster==min(ycluster));
        ycluster_maxindex=find(ycluster==max(ycluster));
    %         if ycluster1>ycluster2
                strr_bound=[strr_bound;points_tmp(cluster_index==ycluster_minindex,:)];
                strr_top=[strr_top;points_tmp(cluster_index==ycluster_maxindex,:)];
    %         else
    %             strr_bound=[strr_bound;points_tmp(cluster_index==1,:)];
    %         end
    
    end

end

function [strr_boundary_depth,flow_depth,flow_depth_normbydp,maxflowdepth_norm,swfit,fitangle_a,maxflowdepth_std]=strr_depth_fun(strr_botbound,strr_topbound,surface_x,surface_y,dp,xwsfit,stepsfit)
xrangefit=[min(strr_botbound(:,1)) max(strr_botbound(:,1))];%fron extreme left point to right point
swfit=[xrangefit(1):stepsfit:xrangefit(2)]';%sample window, set steps to xws to extract boundary points (otherwise overlap points)
fitangle_a=[];
fitangle_b=[];
gof_fitangle=[];
fitdatamidx=[];
fitdatamidy=[];
% %--------------angle fitting settings----------------
ft = fittype( 'a*x+b', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.395894484205228 0.936450101590675];
% %--------------end angle fitting settings----------------

    for fiti=1:size(swfit,1)
    %     :size(swfit,1)
        % Set up fittype and options.
        windowfit=[swfit(fiti)-xwsfit/2 swfit(fiti)+xwsfit/2];
        points_index1=find(strr_botbound(:,1)>windowfit(1));
        points_index2=find(strr_botbound(:,1)<windowfit(2));
        points_index=points_index1(ismember(points_index1,points_index2));
        points_tmp=strr_botbound(points_index,:);
            [fitresult, gof] = fit(points_tmp(:,1),points_tmp(:,2),ft,opts);
            fitangle_a(fiti,1)=atan(fitresult.a);%angle in radian
            fitangle_b(fiti,1)=fitresult.b;%intersection between y axis
            gof_fitangle(fiti,1)=gof.rsquare;
    %         DAOR(samplei,1)=max(fitangle_a);     
    %         DAOR_gof(samplei,1)=gof_fitangle(fitangle_a==max(fitangle_a));
        %% mid point of bottom
        midbot=[swfit(fiti) mean(points_tmp(:,2))];
        %% left and right point based on bottom mid point
        botleft_coor=[midbot(1)-xwsfit/2,midbot(2)-xwsfit/2*fitangle_a(fiti,1)];
        botright_coor=[midbot(1)+xwsfit/2,midbot(2)+xwsfit/2*fitangle_a(fiti,1)];
        %% left and right line and intersection point y axis
        lr_tan=tan(fitangle_a(fiti,1)+pi/2);
        left_b=botleft_coor(2)-botleft_coor(1)*lr_tan;
        right_b=botright_coor(2)-botright_coor(1)*lr_tan;
        left_2nd_coor=[0 left_b];
        right_2nd_coor=[0 right_b];
        %% dist to lines of surface points & extract between lines
        dist_linetoline=sqrt(sum((botleft_coor-botright_coor).^2));%this distance is between two lines perpendicular to the tangent of strain rate boudnary slope
        dist_surf_left_right(:,1)=point_to_line_distance([surface_x surface_y],botleft_coor,left_2nd_coor);
        dist_surf_left_right(:,2)=point_to_line_distance([surface_x surface_y],botright_coor,right_2nd_coor);
        betweenlines_indexleft=find(dist_surf_left_right(:,1)<dist_linetoline);
        betweenlines_indexright=find(dist_surf_left_right(:,2)<dist_linetoline);
        betweenlines_index=betweenlines_indexleft(ismember(betweenlines_indexleft,betweenlines_indexright));
        betweenlines_pointscoor=[surface_x(betweenlines_index) surface_y(betweenlines_index)];
        midtop=[mean(betweenlines_pointscoor(:,1)) mean(betweenlines_pointscoor(:,2))];
        %% dist to lines of top strain rate boundary points & extract between lines
        dist_topstrr_left_right(:,1)=point_to_line_distance([strr_topbound(:,1) strr_topbound(:,2)],botleft_coor,left_2nd_coor);
        dist_topstrr_left_right(:,2)=point_to_line_distance([strr_topbound(:,1) strr_topbound(:,2)],botright_coor,right_2nd_coor);
        betweenlines_indexleft=find(dist_topstrr_left_right(:,1)<dist_linetoline);
        betweenlines_indexright=find(dist_topstrr_left_right(:,2)<dist_linetoline);
        betweenlines_index=betweenlines_indexleft(ismember(betweenlines_indexleft,betweenlines_indexright));
        betweenlines_pointscoor=strr_topbound(betweenlines_index,:);
        midstrrtop=[mean(betweenlines_pointscoor(:,1)) mean(betweenlines_pointscoor(:,2))];
        %% flow region depth; dist between midbot and midtop
        flow_depth(fiti,1)=sqrt(sum((midtop-midbot).^2));
        strr_boundary_depth(fiti,1)=sqrt(sum((midstrrtop-midbot).^2));
        flow_depth_normbydp=flow_depth/dp;
    end
%% smooth the flow depth by average within a window
flow_depth_tmp=flow_depth;
for smi=1:size(flow_depth,1)
    left_bound=swfit(smi,1)-5*dp;
    right_bound=swfit(smi,1)+5*dp;
    flow_depth(smi,1)= mean(flow_depth_tmp(swfit>=left_bound & swfit<=right_bound));
    flow_depth_std(smi,1)=std(flow_depth_tmp(swfit>=left_bound & swfit<=right_bound));
end
maxflowdepth_norm=max(flow_depth)/dp;
flow_depth_std_tmp=flow_depth_std(flow_depth/dp==maxflowdepth_norm);
maxflowdepth_std=flow_depth_std_tmp(1);
% % max by manuel window 
% leftbound=find(swfit/dp>=10,1);
% rightbound=find(swfit/dp>=20,1);
% maxflowdepth_norm=mean(flow_depth_normbydp(leftbound:rightbound));
% %-----------max by window before and after the max value
% flow_depth_normbydp=flow_depth_normbydp(fitangle_a>0);
% max_flowdepth=max(flow_depth_normbydp);
% maxflowdepth_norm=mean(flow_depth_normbydp(find(flow_depth_normbydp==max_flowdepth)-10:find(flow_depth_normbydp==max_flowdepth)+10));

% figure(1);hold on;plot(swfit,fitangle_a/pi*180);
% figure(3);hold on;plot(swfit,flow_depth_normbydp);
end

% % function ------------------- surface std shade ----------------
function [xfit_range2,inBetween,xfit_range,surf_mid]=errorshade(accu_up_surf_points)
surf_std=[];
surf_mid=[];
xfit_range=[-0.1:0.001:0.1]';
for fiti=1:size(xfit_range,1)
index_left=find(accu_up_surf_points(:,1)>=(xfit_range(fiti,1)-0.002));
index_right=find(accu_up_surf_points(:,1)<(xfit_range(fiti,1)+0.002));
index_all=index_right(ismember(index_right,index_left));
    if size(index_all,1)==0
        surf_mid(fiti,1)=NaN;
        surf_std(fiti,1)=NaN;
    else
        surf_mid(fiti,1)=mean(accu_up_surf_points(index_all,2));
        surf_std(fiti,1)=std(accu_up_surf_points(index_all,2));
    end
end
% ----------extract points in R<0.085-------------
dist_to00=sqrt(xfit_range.^2+surf_mid.^2);
xfit_range=xfit_range(dist_to00<0.085);
surf_mid=surf_mid(dist_to00<0.085);
surf_std=surf_std(dist_to00<0.085);
% ---------plot shade std -----------
%         y = rand(1,10); % your mean vector;
%         x = 1:numel(y);
%         std_dev = 1;
curve1 = surf_mid + surf_std;
curve2 = surf_mid - surf_std;
xfit_range2 = [xfit_range(~isnan(curve1)); flipud(xfit_range(~isnan(curve1)))];
inBetween = [curve1(~isnan(curve1),:); flipud(curve2(~isnan(curve2),:))];
end
% figure(5);hold on;
% load('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis\C0083_strainrate_depth.mat')
% plot(swfit,flow_depth_normbydp,'c.');
% load('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis\C0082_strainrate_depth.mat')
% plot(swfit,flow_depth_normbydp,'g.');
% load('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis\C0081_strainrate_depth.mat')
% plot(swfit,flow_depth_normbydp,'r.');
% load('E:\RoDrtest\RoDr_alphastart0-1\model45\model45_strainrate_depth.mat');
% plot(swfit,flow_depth_normbydp,'co');
% load('E:\RoDrtest\RoDr_alphastart0-1\model46\model46_strainrate_depth.mat');
% plot(swfit,flow_depth_normbydp,'go');
% load('E:\RoDrtest\RoDr_alphastart0-1\model47\model47_strainrate_depth.mat');
% plot(swfit,flow_depth_normbydp,'ro');
% set(gca,'FontSize',12)
% xlabel('x');
% ylabel('Normalised maximum depth (depth/d_p)');
% ylabel('Normalised depth (depth/d_p)');