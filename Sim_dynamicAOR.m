clc;
clearvars -except We;
%%
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\dataanalysis\');
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis\function');
[num,txt,raw]=xlsread(['C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\parametric_study.xlsx'],1);
keypara=num(4:end,:);
density_threshold=0.6;

% modelnum=[31,33,35,37,39,41]';
% modelnum=[33,53,54,59:69]';
%% plot cluster size vs variables
% model_var="omega";
% model_var="mu_s";
% model_var="mu_r";
model_var='gamma';
% model_var='liquid_content'
% model_var='all';

if model_var=="omega"
modelnum=[86:89]';% vary rotation speed
elseif model_var == "mu_s"
modelnum=[90,87,91]';% mu_s
elseif model_var == "mu_r"
modelnum=[92,87,93,94]';% mu_r
elseif model_var == "gamma"
% modelnum=[103,95,96,87,97]';% vary surface tension gamma
modelnum=[119,120,134,121:123]';% w90
% modelnum=[124,125,151,126:128]';% w120
% modelnum=[139,140,152,141:143]';% w150
elseif model_var == "liquid_content"
modelnum=[103,87,99:102]';% vary liquid content
elseif model_var == "all"
% modelnum=[95,96,87,97]';
% modelnum=[86:89,95:97]';
% modelnum=[81,82]';
% modelnum=[49:52,70]';
% modelnum=[112,103,109:111]';
% modelnum=[112:114,92]';%dry mu_r
% modelnum=[129,121:123]';%dry mu_r
% modelnum=[124,126]';%compare with exp w=120
modelnum=[119,121]';%compare with exp w=90
end
% modelnum=[77];
%% results record
DAORt_den=[];
DAORt_ras=[];
goft_den=[];
goft_raw=[];
revolution=[];
DAOR=[];
DAOR_gof=[];

%% custom colormap
color=colormap(turbo(8));

%%
for modeli=1:size(modelnum,1)

%% marker & shift
if modelnum(1,1)==119
    marker="-";
    shift=-10;
    column=1;
elseif modelnum(1,1)==124
    marker="--";
    shift=60;
    column=2;
elseif modelnum(1,1)==139
    marker=":";
    shift=90;
    column=3;
end


datafile=(['f:\RoDrtest\RoDr_alphastart0-1\model' num2str(modelnum(modeli,1))]);
load([datafile '\model' num2str(modelnum(modeli,1)) '_up_surface_points.mat']);
% load([datafile '\model' num2str(modelnum(modeli,1)) '_velmatrix_frames.mat']);
accu_up_surf_points=[];
for i=1:size(up_surf_points,1)
    accu_up_surf_points=[accu_up_surf_points;up_surf_points{i,1}];
end
% plot(accu_up_surf_points(:,1),accu_up_surf_points(:,2),'.');
% [f,ep]=ksdensity([bxtmp bytmp],xi); % remove the outputs to see a 3D plot of the distribution
DEMtimestep=10^-7;
DEMstepinterval=50000;

omega=keypara(modelnum(modeli,1),5);%revolution/s
dp=keypara(modelnum(modeli,1),7);
V_p=4/3*pi*(dp/2)^3;
Rp=dp/2;
Rcon=0.1;%container radius
Centercir=[0,0];
container_d=0.03;
phi_sample_meshsize=1*dp;
mur=keypara(modelnum(modeli,1),10);
mus=keypara(modelnum(modeli,1),9);
mur=keypara(modelnum(modeli,1),10);
mus=keypara(modelnum(modeli,1),9);
surften=keypara(modelnum(modeli,1),11);
liquidcontent=keypara(modelnum(modeli,1),6);
rho=2460;%kg/m3
% vave_max=sqrt(max(meanvelxincell,[],"all")^2+max(meanvelzincell,[],"all")^2);
We(modeli,column)=rho*dp/2*(omega/180*pi*0.1)^2/surften;
% We(modeli,1)=rho*dp/2*(vave_max)^2/surften;
% We(modeli,2)=modelnum(modeli,1);


        %% ------------------- surface std shade ----------------
        [xfit_range2,inBetween,xfit_range_extracted,surf_mid]=errorshade(accu_up_surf_points);

        %% ---------------fitting DAOR-------------------
        % Set up fittype and options.
        ft = fittype( 'a*x+b', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.StartPoint = [0.395894484205228 0.936450101590675];
        fitangle_a=[];
        fitangle_b=[];
        gof_fitangle=[];
        fitdatamidx=[];
        fitdatamidy=[];
        num_points_to_fit=20;
        xfit_range=[-0.1:0.001:0.1]';
        for fiti=1:size(xfit_range,1)
        index_left=find(accu_up_surf_points(:,1)>=(xfit_range(fiti,1)-0.01));
        index_right=find(accu_up_surf_points(:,1)<(xfit_range(fiti,1)+0.01));
        index_all=index_right(ismember(index_right,index_left));
            if size(index_all,1)<2
                fitangle_a(fiti,1)=NaN;%dynamic angle of repose in radius
                fitangle_b(fiti,1)=NaN;%dynamic angle of repose in degree
                gof_fitangle(fiti,1)=NaN;
            else
                [fitresult, gof] = fit(accu_up_surf_points(index_all,1),accu_up_surf_points(index_all,2),ft,opts);
                fitangle_a(fiti,1)=round(atan(fitresult.a)/pi*180,1);%dynamic angle of repose in radius
                fitangle_b(fiti,1)=fitresult.b;%dynamic angle of repose in degree
                gof_fitangle(fiti,1)=gof.rsquare;
            end
        end
%         DAOR(modeli,1)=max(fitangle_a);
          DAOR(modeli,1)=round(mean(fitangle_a(find(xfit_range>0,1):find(xfit_range>0.03,1))),1);
          DAOR_tmp=fitangle_a;
            for smi=1:size(DAOR_tmp,1)
                left_bound=xfit_range(smi,1)-5*dp;
                right_bound=xfit_range(smi,1)+5*dp;
                DAOR_mean(smi,1)= mean(DAOR_tmp(xfit_range>=left_bound & xfit_range<=right_bound));
                DAOR_std(smi,1)= std(DAOR_tmp(xfit_range>=left_bound & xfit_range<=right_bound));
            end
          DAOR_mean_max(modeli,1)=max(DAOR_mean);
          DAOR_std_tmp=DAOR_std(DAOR_mean==DAOR_mean_max(modeli,1));
          DAOR_std_max(modeli,1)=DAOR_std_tmp(1);
%         DAOR_gof(modeli,1)=gof_fitangle(fitangle_a==max(fitangle_a));

        %  %-------plot---------
        
        figure(1);hold on;
        min_xdistance=find(xfit_range>min(accu_up_surf_points(:,1)),1);
        max_xdistance=find(xfit_range>max(accu_up_surf_points(:,1)),1);
        if model_var=="all"
        color=colormap(turbo(36));
        figure(3);hold on;
        fig3=plot(1/We(modeli,1),DAOR(modeli,1), ...
            'k.', ...
            'MarkerSize',10);
%%
        figure(4);hold on;%plot surface profile with shade std area
        fig4=plot(xfit_range_extracted/dp,surf_mid/dp, ...
            '--', ...
            'LineWidth',1, ...
            'Color',color(modeli+1,:));
        fill(xfit_range2/dp, inBetween/dp, ...
            color(modeli+1,:), ...
            'FaceAlpha',0.3,'LineStyle','none');
        hold on;
        figure(5);hold on;%plot slope transition from left to right
        fig5=plot(xfit_range(min_xdistance:max_xdistance,1),fitangle_a(min_xdistance:max_xdistance,1), ...
            '.', ...
            'MarkerSize',10, ...
            'Color',color(modeli+1,:));
        else
        figure(1);hold on;%plot surface profile with shade std area
        fig1=plot(xfit_range_extracted/dp+shift,surf_mid/dp, ...
            marker, ...
            'LineWidth',1.5, ...
            'Color',color(modeli+1,:));
        fill(xfit_range2/dp+shift, inBetween/dp, ...
            color(modeli+1,:), ...
            'FaceAlpha',0.3,'LineStyle','none');
        hold on;
        figure(2);hold on;%plot slope transition from left to right
        fig2=plot(xfit_range(min_xdistance:max_xdistance,1),fitangle_a(min_xdistance:max_xdistance,1), ...
            '.', ...
            'MarkerSize',10, ...
            'Color',color(modeli+1,:));
        end
        

        % % -----plot----------
%         figure(2);
% %         hold on;
% %         plot(boundline(up_surf_points(10:end-4),1),boundline(up_surf_points(10:end-4),2),'ko');
% %         xlim([0 0.12]);
% %         ylim([-0.1 0.05]);
%         hold on
%         xft=[-0.1:0.01:0.1];
%         yft=tan(max(fitangle_a))*xft+fitangle_b(fitangle_a==max(fitangle_a));
%         plot(xft,yft,'r-','LineWidth',2);
%         ylim([-0.14 0.05]);
%         xlim([-0.1 0.1])
% 
% %         figure(3)
% %         hold on
% %         plot(boundline(up_surf_points([1:end-num_points_to_fit+1]+2),1),fitangle_a);
        %% --------plot max angle point-------------
%         coord_maxanglex=fitdatamidx(fitangle_a==max(fitangle_a));
%         coord_maxangley=fitdatamidy(fitangle_a==max(fitangle_a));
%         plot(coord_maxanglex,coord_maxangley,'go','MarkerSize',10,'LineWidth',1);
        
        %% ----------plot perpendicular line of the fitting line
%         currslope=tan(max(fitangle_a));
%         perpslope=-1/currslope;
%         perp_b=fitdatamidy(fitangle_a==max(fitangle_a))-(perpslope*fitdatamidx(fitangle_a==max(fitangle_a)));
% %         figure(2); hold on;
% %         yft_perp=perpslope*xft+perp_b;
% %         plot(xft,yft_perp,'g-','LineWidth',1);
        %% ----------plot intersection between perpendicular line and circle
% %         figure(2);hold on;
%         [interx,intery] = linecirc(perpslope,perp_b,Centercir(1),Centercir(2),Rcon);
% %         plot(interx(1),intery(1),'y*','MarkerSize',10,'LineWidth',1);
%         dist_OI=sqrt((coord_maxanglex-interx(1)).^2+(coord_maxangley-intery(1)).^2);
        %% -----plot rotated data according to AOR point---------------
% %         coordtmp=[reshape(xextacu,[],1) (reshape(yextacu,[],1)-1080)*-1];% (data-1080)*-1 means a reverse in y direction
%         rota_datacoord= rotation([tmpx tmpz],[coord_maxanglex,coord_maxangley],-DAOR(samplei,1)/pi*180);
% %         after_rotx=reshape(rota_datacoord(:,1),size(xextacu,1),size(xextacu,2));
% %         after_roty=(reshape(rota_datacoord(:,2),size(xextacu,1),size(xextacu,2)));
%         figure(7);
% %         quiver(after_rotx,after_roty,uoriginal_extacu,voriginal_extacu,3);
%         scatter(rota_datacoord(:,1),rota_datacoord(:,2));
%         hold on; plot(coord_maxanglex,coord_maxangley,'go','MarkerSize',10,'LineWidth',1);
%         axis equal;
%         ylim([-0.14 0.05]);
%         xlim([-0.1 0.1]);
%     %     set(gca, 'ydir', 'reverse');
        %% ---------plot sampling region---------
%         figure(7);hold on;
%         atomrange=d_p*5;%the assumed pixel for one particle
%         S1(1,1)=coord_maxanglex+atomrange;S1(1,2)=coord_maxangley;
%         S2(1,1)=coord_maxanglex-atomrange;S2(1,2)=coord_maxangley;
%         S3(1,1)=coord_maxanglex+atomrange;S3(1,2)=coord_maxangley-dist_OI;
%         S4(1,1)=coord_maxanglex-atomrange;S4(1,2)=coord_maxangley-dist_OI;
%         S_coord=[S1;S2;S4;S3];
%         plot(S_coord(:,1),S_coord(:,2),'o','MarkerSize',10,'LineWidth',1,'MarkerFaceColor',[0.9290 0.6940 0.1250]);
%         plot(polyshape(S_coord(:,1),S_coord(:,2)));
        %% -----extract velocity profile--------------
%         
%         vel_dircorr=ones(size(tmpvel,1),1);%velocity direction correction
%         vel_dircorr(tmpvel(:,1)>0)=-1;
%         vel_totacu=sqrt(tmpvel(:,1).^2+tmpvel(:,3).^2)/sqrt(9.81*d_p).*vel_dircorr;
%         vel_isrx_index=intersect(find(rota_datacoord(:,1)<S1(1)),find(rota_datacoord(:,1)>S2(1)));%extract data in sampling region
%         vel_isry_index=intersect(find(rota_datacoord(:,2)<S1(2)),find(rota_datacoord(:,2)>S3(2)));%extract data in sampling region
%         vel_isr_index=intersect(vel_isrx_index,vel_isry_index);
%         vel_isr_coord=rota_datacoord(vel_isr_index,:);
%         vel_isr=vel_totacu(vel_isr_index,:);
%     %     figure(8);quiver(after_rotx,after_roty,uoriginal_extacu,voriginal_extacu,3);
%         figure(9)
%         subplot(2,1,1);
%         scatter((vel_isr_coord(:,2)-S1(2))/d_p,vel_isr_coord(:,1)/d_p,5,vel_isr);
%         xlim([-30 0]);
%         axis equal;
%         wsv=1*d_p;%window size vertical; unit = d_p particle diameter
%         wsh=S1(1)-S2(1);
%         movingsteps=ceil((S1(2)-S3(2))/wsv);
%         vel_el=zeros(movingsteps,1);%velocity each layer
%         middleline=zeros(movingsteps,1);
%         for velisri=1:movingsteps
%             upperbound=S3(2)+velisri*wsv;
%             bottombound=S3(2)+(velisri-1)*wsv;
%             middleline(velisri,1)=(upperbound+bottombound)*0.5;
%             veltmp_extract=[];
%             veltmp_extract=vel_isr(intersect(find(vel_isr_coord(:,2)>bottombound),find(vel_isr_coord(:,2)<upperbound)));
%             veltmp_extract(isnan(veltmp_extract))=[];
%             vel_el(velisri,1)=mean(veltmp_extract);
% 
%         end
%         subplot(2,1,2);
%         plot((middleline-S1(2))/d_p,vel_el);
%         xlim([-30 0]);
%         

    %% ===========save results=============
%    cd(datafile);
%    datafilename=['model' num2str(modelnum(modeli,1)) '_DAOR_gof.mat'];
%    save(datafilename,'DAOR');
%    save(datafilename,'DAOR_gof','-append');
% %    save(datafilename,'DAORt_fitraw','-append');
% %    save(datafilename,'goft_raw','-append');
%    save(datafilename,'revolution','-append');
%    clear DAOR DAOR_gof revolution
   
end
%% figure settings
figure(1)
xlabel('\itx/d');
ylabel('\ity/d');
xlim([-50 130]);
ylim([-60 60]);
box on;
set(gca,'LineWidth',2);
set(gca,'FontSize',18);
set(gca,'XTick',[-40,-20,0,20,40,60,80,100,120]);
set(gca,'FontName','times new roman')
if model_var=="omega"
legend(['\omega=60^{\circ}/s; DAOR=' num2str(DAOR(1))], ...
    ['\omega=90^{\circ}/s; DAOR=' num2str(DAOR(2))], ...
    ['\omega=120^{\circ}/s; DAOR=' num2str(DAOR(3))], ...
    ['\omega=150^{\circ}/s; DAOR=' num2str(DAOR(4))], ...
    'Location','south', ...
    'fontsize',12);
elseif model_var=="mu_s"
legend(['\mu_s=0.1; DAOR=' num2str(DAOR(1))], ...
    ['\mu_s=0.5; DAOR=' num2str(DAOR(2))], ...
    ['\mu_s=0.8; DAOR=' num2str(DAOR(3))], ...
    'Location','south', ...
    'fontsize',12);
elseif model_var=="mu_r"
legend(['\mu_r=0.001; DAOR=' num2str(DAOR(1))], ...
    '', ...
    ['\mu_r=0.01; DAOR=' num2str(DAOR(2))], ...
    '', ...
    ['\mu_r=0.1; DAOR=' num2str(DAOR(3))], ...
    '', ...
    ['\mu_r=1; DAOR=' num2str(DAOR(4))], ...
    'Location','northwest', ...
    'fontsize',12);
elseif model_var=="gamma"
legend(['\gamma = 0'], ...
    '', ...
    ['\gamma = 0.0073 \itN/m'], ...
    '', ...
    ['\gamma = 0.0128 \itN/m'], ...
    '', ...
    ['\gamma = 0.0365 \itN/m'], ...
    '', ...
    ['\gamma = 0.073 \itN/m'], ...
    '', ...
    ['\gamma = 0.146 \itN/m'], ...
    '', ...
    'Location','northwest', ...
    'fontsize',12, ...
    'numColumns',2);
%     ['\itWe = ' num2str(round(We(1,1),3,'significant'))], ...
%     '', ...
%     ['\itWe = ' num2str(round(We(2,1),3,'significant'))], ...
%     '', ...
%     ['\itWe = ' num2str(round(We(3,1),3,'significant'))], ...
%     '', ...
%     ['\itWe = ' num2str(round(We(4,1),3,'significant'))], ...
%     '', ...
%     ['\itWe = ' num2str(round(We(5,1),2,'significant'))], ...
%     '', ...
%     ['\itWe = ' num2str(round(We(6,1),2,'significant'))], ...
%     '', ...
%     ['\itWe = ' num2str(round(We(1,2),3,'significant'))], ...
%     '', ...
%     ['\itWe = ' num2str(round(We(2,2),3,'significant'))], ...
%     '', ...
%     ['\itWe = ' num2str(round(We(3,2),3,'significant'))], ...
%     '', ...
%     ['\itWe = ' num2str(round(We(4,2),3,'significant'))], ...
%     '', ...
%     ['\itWe = ' num2str(round(We(5,2),3,'significant'))], ...
%     '', ...
%     ['\itWe = ' num2str(round(We(6,2),2,'significant'))], ...
%     '', ...
%     ['\itWe = ' num2str(round(We(1,3),3,'significant'))], ...
%     '', ...
%     ['\itWe = ' num2str(round(We(2,3),3,'significant'))], ...
%     '', ...
%     ['\itWe = ' num2str(round(We(3,3),3,'significant'))], ...
%     '', ...
%     ['\itWe = ' num2str(round(We(4,3),3,'significant'))], ...
%     '', ...
%     ['\itWe = ' num2str(round(We(5,3),3,'significant'))], ...
%     '', ...
%     ['\itWe = ' num2str(round(We(6,3),2,'significant'))], ...
%     '', ...

elseif model_var=="liquid_content"
legend(['dry; DAOR=' num2str(DAOR(1))], ...
    ['m=0.02; DAOR=' num2str(DAOR(2))], ...
    ['m=0.04; DAOR=' num2str(DAOR(3))], ...
    ['m=0.06; DAOR=' num2str(DAOR(4))], ...
    ['m=0.08; DAOR=' num2str(DAOR(5))], ...
    ['m=0.1; DAOR=' num2str(DAOR(6))], ...
    'Location','south', ...
    'fontsize',12);
end

figure(2)
xlabel('\it{x/d}');
ylabel('Free surface profile');
% xlim([-0.1 0.1]);
% ylim([-0.1 0.06]);
box on;
set(gca,'LineWidth',1);
set(gca,'FontSize',18);
if model_var=="omega"
legend(['\omega=60^{\circ}/s; DAOR=' num2str(DAOR(1))], ...
    ['\omega=90^{\circ}/s; DAOR=' num2str(DAOR(2))], ...
    ['\omega=120^{\circ}/s; DAOR=' num2str(DAOR(3))], ...
    ['\omega=150^{\circ}/s; DAOR=' num2str(DAOR(4))], ...
    'Location','southeast', ...
    'fontsize',12);
elseif model_var=="mu_s"
legend(['\mu_s=0.1; DAOR=' num2str(DAOR(1))], ...
    ['\mu_s=0.5; DAOR=' num2str(DAOR(2))], ...
    ['\mu_s=0.8; DAOR=' num2str(DAOR(3))], ...
    'Location','southeast', ...
    'fontsize',12);
elseif model_var=="mu_r"
legend(['\mu_r=0.001; DAOR=' num2str(DAOR(1))], ...
    ['\mu_r=0.01; DAOR=' num2str(DAOR(2))], ...
    ['\mu_r=0.1; DAOR=' num2str(DAOR(3))], ...
    ['\mu_r=1; DAOR=' num2str(DAOR(4))], ...
    'Location','southeast', ...
    'fontsize',12);
elseif model_var=="gamma"
legend(['dry; DAOR=' num2str(DAOR(1))], ...
    '', ...
    ['\gamma=0.0365N/m; DAOR=' num2str(DAOR(2))], ...
    '', ...
    ['\gamma=0.073N/m; DAOR=' num2str(DAOR(3))], ...
    '', ...
    ['\gamma=0.146N/m; DAOR=' num2str(DAOR(4))], ...
    '', ...
    'Location','southeast', ...
    'fontsize',12);
elseif model_var=="liquid_content"
legend(['dry; DAOR=' num2str(DAOR(1))], ...
    ['m=0.02; DAOR=' num2str(DAOR(2))], ...
    ['m=0.04; DAOR=' num2str(DAOR(3))], ...
    ['m=0.06; DAOR=' num2str(DAOR(4))], ...
    ['m=0.08; DAOR=' num2str(DAOR(5))], ...
    ['m=0.1; DAOR=' num2str(DAOR(6))], ...
    'Location','southeast', ...
    'fontsize',12);
end
set(gca,'FontName','Times New Roman');
%% ===========Exp profile plot==============
% % % extract calibration and Container center info
% maskdir='C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis\pivlab_masks\';
% % if filenum(fi,1)<156
% %     CCR=load([maskdir 'C0149-C0155_CR.mat']);
% %     calibration_coe=load([maskdir 'calibration_coe_c0149-c0155.mat']);
% % elseif filenum(fi,1)>155 && filenum(fi,1)<177
%     CCR=load([maskdir 'C0157-C0172_CR.mat']);
%     calibration_coe=load([maskdir 'calibration_coe_c0157-c0172.mat']);
% % elseif filenum(fi,1)>176
% %     CCR=load([maskdir 'C0177-C0182_CR.mat']);
% %     load([maskdir 'calibration_coe_c0177-c0182.mat']);
% % end
% CCreal=CCR.CC*calibration_coe.pixeltorealdimension;
% CCreal(1,2)=(CCreal(1,2)-1080*calibration_coe.pixeltorealdimension)*-1;
% 
% figure(4)
% exp_flowdepth=load('E:\EXP RoDrtest\Exp\sharpened\C0157\C0157_strr-coe80_strainrate_depth_2.mat');
% % plot(exp.surface_x_real/dp,exp.surface_y_real/dp,'k-');
% exp=load('E:\EXP RoDrtest\Exp\sharpened\C0157\C0157_particleboundary_2.mat');
% accu_up_surf_points_exp=[];
% bx=exp.bx;
% by=exp.by;
% for i=1:size(bx,1)
%     accu_up_surf_points_exp=[accu_up_surf_points_exp;[bx{i,1},by{i,1}]];
% end
% accu_up_surf_points_exp(:,1)=accu_up_surf_points_exp(:,1)-CCreal(1);
% accu_up_surf_points_exp(:,2)=accu_up_surf_points_exp(:,2)-CCreal(2);
% [xfit_range2_exp,inBetween_exp,xfit_range,surf_mid]=errorshade(accu_up_surf_points_exp);
% 
% plot(xfit_range/dp,surf_mid/dp+1, ...
%             'r-', ...
%             'LineWidth',1);
% fill(xfit_range2_exp/dp, inBetween_exp/dp+1, ...
%     'r', ...
%     'FaceAlpha',0.3,'LineStyle','none');
%% ===========================================
% % % extract calibration and Container center info
% maskdir='C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis\pivlab_masks\';
% % if filenum(fi,1)<156
%     CCR=load([maskdir 'C0149-C0155_CR.mat']);
%     calibration_coe=load([maskdir 'calibration_coe_c0149-c0155.mat']);
% % elseif filenum(fi,1)>155 && filenum(fi,1)<177
% %     CCR=load([maskdir 'C0157-C0172_CR.mat']);
% %     calibration_coe=load([maskdir 'calibration_coe_c0157-c0172.mat']);
% % elseif filenum(fi,1)>176
% %     CCR=load([maskdir 'C0177-C0182_CR.mat']);
% %     load([maskdir 'calibration_coe_c0177-c0182.mat']);
% % end
% CCreal=CCR.CC*calibration_coe.pixeltorealdimension;
% CCreal(1,2)=(CCreal(1,2)-1080*calibration_coe.pixeltorealdimension)*-1;
% 
% figure(4)
% exp_flowdepth=load('E:\EXP RoDrtest\Exp\sharpened\C0152\C0152_strr-coe80_strainrate_depth_2.mat');
% % plot(exp.surface_x_real/dp,exp.surface_y_real/dp,'k-');
% exp=load('E:\EXP RoDrtest\Exp\sharpened\C0152\C0152_particleboundary_2.mat');
% accu_up_surf_points_exp=[];
% bx=exp.bx;
% by=exp.by;
% for i=1:size(bx,1)
%     accu_up_surf_points_exp=[accu_up_surf_points_exp;[bx{i,1},by{i,1}]];
% end
% accu_up_surf_points_exp(:,1)=accu_up_surf_points_exp(:,1)-CCreal(1);
% accu_up_surf_points_exp(:,2)=accu_up_surf_points_exp(:,2)-CCreal(2);
% [xfit_range2_exp,inBetween_exp,xfit_range,surf_mid]=errorshade(accu_up_surf_points_exp);
% 
% plot(xfit_range/dp,surf_mid/dp+2, ...
%             'k-', ...
%             'LineWidth',1);
% fill(xfit_range2_exp/dp, inBetween_exp/dp+2, ...
%     'k', ...
%     'FaceAlpha',0.3,'LineStyle','none');



    %% colormap
%     colormap(flipud(gray));
%     whiteTOblue=[flip([0:0.01:1]'),flip([0:0.01:1]'),flip([0.58:0.0042:1]')];
%     newmap=[whiteTOblue;jet(2000)];
%     colormap(newmap);
%     hold on
%     
%     set(gca, 'XScale', 'log');
%     set(gca, 'YScale', 'log');
% 
%     ax=gca;
% %     ax.XAxis.Exponent = 0;
% %     xtickformat('%.4f')
%     xlabel('\Delta\phi')
%     ylabel('T')
%     colorbar
%% tmp legend
% figure(1)
% legend('\omega=15\itrpm', ...
% '', ...
% '', ...
% '', ...
% '', ...
% '', ...
% '', ...
% '', ...
% '\omega=20\itrpm', ...
% '', ...
% '', ...
% '', ...
% '', ...
% '', ...
% '', ...
% '', ...
% 'Location','northwest', ...
% 'fontsize',12);

%% function ------------------- surface std shade ----------------
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