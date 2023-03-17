clc;
clear;
%%
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\dataanalysis\');
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis\function');
[num,txt,raw]=xlsread(['C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\parametric_study.xlsx'],1);
keypara=num(4:end,:);
density_threshold=0.6;
% modelnum=[31,33,35,37,39,41]';
modelnum=[45]';
%%
for modeli=1:size(modelnum)
datafile=(['g:\RoDrtest\RoDr_alphastart0-1\model' num2str(modelnum(modeli,1))]);
load([datafile '\model' num2str(modelnum(modeli,1)) '_posi.mat']);
load([datafile '\model' num2str(modelnum(modeli,1)) '_velocity.mat']);
load([datafile '\model' num2str(modelnum(modeli,1)) '_id.mat']);
[Rtxt,Ctxt]=size(centroids);
DEMtimestep=10^-7;
DEMstepinterval=50000;
DAORt_den=[];
DAORt_ras=[];
goft_den=[];
goft_raw=[];
revolution=[];
DAOR=[];
DAOR_gof=[];
omega=keypara(modelnum(modeli,1),5);
d_p=keypara(modelnum(modeli,1),7);
V_p=4/3*pi*(d_p/2)^3;
Rp=d_p/2;
Rcon=0.1;%container radius
Centercir=[0,0];
container_d=0.03;
phi_sample_meshsize=1*d_p;
sample_range=[50:10:Rtxt-10]';
% strain_rate_threshold=1000;
e_sample_r=10;
accu_steps=10;
pocesstime=[1:accu_steps]'/accu_steps;
strain_rate_cell=cell(accu_steps,1);
%% ======= prepare the mesh region to extract data and do analysis =========
meshbound=[-0.12 0.12];
xgridpf=linspace(meshbound(1),meshbound(2),(meshbound(2)-meshbound(1))/phi_sample_meshsize)';
zgridpf=linspace(meshbound(1),meshbound(2),(meshbound(2)-meshbound(1))/phi_sample_meshsize)';
[Rxgri,Cxgri]=size(xgridpf);
[Rzgri,Czgri]=size(zgridpf);
Xpf=repmat(xgridpf(1:Rxgri-1,1),1,Rxgri-1);
Zpf=repmat(zgridpf(1:Rxgri-1,1)',Rxgri-1,1);


%%
% sample_range=[300]';
    for samplei=size(sample_range,1)-40
%         1:size(sample_range,1)
    up_surf_p_posiaccu=[];
    revolution(samplei,1)=sample_range(samplei)*DEMstepinterval*DEMtimestep*omega/360;% for model 25-30 sample_range-24
    strain_rate_accu=zeros(Rxgri-2,Rxgri-2);
    srms=size(strain_rate_accu,1);%strain rate matrix size
        %% prepare data for timestep accumulation
        delta_vx=zeros(Rxgri-2,Rxgri-2);
        delta_vz=zeros(Rxgri-2,Rxgri-2);
        delta_posix=zeros(Rxgri-2,Rxgri-2);
        delta_posiz=zeros(Rxgri-2,Rxgri-2);
        accu_average_counts=zeros(Rxgri-2,Rxgri-2);%count how many frames the cell has particles
        %%
        for esri=1:accu_steps
            pocesstime(esri,1)
            
                tmpx_step1=[];
                tmpz_step1=[];
                tmpvel_step1=[];
                tmpid_step1=[];
                idpovel_step2=[];% id,position,velocity record of step 2
                tmpx_step2=[];
                tmpz_step2=[];
                tmpvel_step2=[];
                tmpid_step2=[];
%             :e_sample_r %each sampling range

%         atom_posi_eachstep=cell2mat(centroids(sample_range(samplei)+esri-5,1));
%         atom_vel_eachstep=cell2mat(velocity(sample_range(samplei)+esri-5,1));
%         tmpx_step1=[tmpx_step1;atom_posi_eachstep(:,1)];%accumultive x coordinates from several steps
%         tmpz_step1=[tmpz_step1;atom_posi_eachstep(:,3)];
%         atom_vel_eachstep=cell2mat(velocity(sample_range(samplei)+esri-5,1));
%         tmpvel_step1=[tmpvel_step1;atom_vel_eachstep];

        % %---------mark particles timestep 1-----------
        
        atom_posi_step1=cell2mat(centroids(sample_range(samplei)+esri-1,1));
        atom_vel_step1=cell2mat(velocity(sample_range(samplei)+esri-1,1));
        atom_posi_step2=cell2mat(centroids(sample_range(samplei)+esri,1));
        atom_vel_step2=cell2mat(velocity(sample_range(samplei)+esri,1));
        
        tmpx_step1=[tmpx_step1;atom_posi_step1(:,1)];%accumultive x coordinates from several steps
        tmpz_step1=[tmpz_step1;atom_posi_step1(:,3)];
        tmpx_step2=[tmpx_step2;atom_posi_step2(:,1)];
        tmpz_step2=[tmpz_step2;atom_posi_step2(:,3)];
        
        tmpvel_step1=[tmpvel_step1;atom_vel_step1];
        tmpvel_step2=[tmpvel_step2;atom_vel_step2];
        
        atom_id_step1=cell2mat(id(sample_range(samplei)+esri-1,1));
        atom_id_step2=cell2mat(id(sample_range(samplei)+esri,1));
        tmpid_step1=[tmpid_step1;atom_id_step1];
        tmpid_step2=[tmpid_step2;atom_id_step2];

        %% convexhull method
        % [k,av]=convhull(tmpx_step1,tmpz_step1,'simplify',true);
        % plot(tmpx_step1(k),tmpz_step1(k),'r-',tmpx_step1,tmpz_step1,'b*')

        %% boundary method
        % k=boundary(tmpx_step1,tmpz_step1,0.8);
        % plot(tmpx_step1(k),tmpz_step1(k),'r-',tmpx_step1,tmpz_step1,'b*')

        %% top surface method
        % boundL=min(tmpx_step1);
        % boundR=max(tmpx_step1);
        % LRrange=boundR-boundL;
        % divLR=LRrange/10;

%         end
        idpovel_step2=[tmpid_step2,tmpx_step2,tmpz_step2,tmpvel_step2(:,1),tmpvel_step2(:,3)];
        idpovel_step2=sortrows(idpovel_step2,1);
        veltotalscaler=sqrt(sum(tmpvel_step1.^2,2));
        %% ============= plot 2D particles ===============
%         figure(2)
%         particle_plot=scatter(tmpx_step1,tmpz_step1);
        %% =================dynamic Angle of Repose: packing fraction and total velocity======================
%         tic
%         meshbound=[-0.12 0.12];
%         xgridpf=linspace(meshbound(1),meshbound(2),(meshbound(2)-meshbound(1))/phi_sample_meshsize)';
%         zgridpf=linspace(meshbound(1),meshbound(2),(meshbound(2)-meshbound(1))/phi_sample_meshsize)';
%         [Rxgri,Cxgri]=size(xgridpf);
%         [Rzgri,Czgri]=size(zgridpf);
        Nincell=zeros(size(xgridpf,1)-1);
        veltincell=zeros(size(xgridpf,1)-1);%scaler velocity in cell
        velxincell_step1=zeros(size(xgridpf,1)-1);
        velzincell_step1=zeros(size(xgridpf,1)-1);
        posixincell_step1=zeros(size(xgridpf,1)-1);
        posizincell_step1=zeros(size(xgridpf,1)-1);
        idincell_step1=cell(size(xgridpf,1)-1,1);% mark id of particles in cell
        for xgi=1:Rxgri-1
                for zgi=1:Rzgri-1
%                     xrange=[];
%                     zrange=[];
%                     xrange=intersect(find(tmpx_step1<xgridpf(xgi+1,1)),find(tmpx_step1>xgridpf(xgi,1)));
%                     zrange=intersect(find(tmpz_step1<zgridpf(zgi+1,1)),find(tmpz_step1>zgridpf(zgi,1)));
                    xrangea=find(tmpx_step1<xgridpf(xgi+1,1));
                    xrangeb=find(tmpx_step1>xgridpf(xgi,1));
                    zrangea=find(tmpz_step1<zgridpf(zgi+1,1));
                    zrangeb=find(tmpz_step1>zgridpf(zgi,1));
                    xrange=xrangea(ismembc(xrangea,xrangeb));
                    zrange=zrangea(ismembc(zrangea,zrangeb));
%                     inters_id=intersect(xrange,zrange);
                    inters_rownum=xrange(ismembc(xrange,zrange));
                    Nincell(xgi,zgi)=size(inters_rownum,1);
                    if inters_rownum>0
                        accu_average_counts(xgi,zgi)=accu_average_counts(xgi,zgi)+1;
                        veltincell(xgi,zgi)=mean(veltotalscaler(inters_rownum,1));
                        idincell_step1{xgi,zgi}=tmpid_step1(inters_rownum);
                        posixincell_step1(xgi,zgi)=mean(tmpx_step1(inters_rownum));
                        posizincell_step1(xgi,zgi)=mean(tmpz_step1(inters_rownum));
                        velxincell_step1(xgi,zgi)=mean(tmpvel_step1(inters_rownum,1));
                        velzincell_step1(xgi,zgi)=mean(tmpvel_step1(inters_rownum,3));
                    else
                        idincell_step1{xgi,zgi}=[];
                    end
                end
        end
%         Xpf=repmat(xgridpf(1:Rxgri-1,1),1,Rxgri-1);
%         Zpf=repmat(zgridpf(1:Rxgri-1,1)',Rxgri-1,1);
        phi_incell=V_p.*Nincell/(abs(xgridpf(2)-xgridpf(1))*abs(zgridpf(2)-zgridpf(1))*container_d*e_sample_r);
        % the circle to extract upper surface particles
        Rcir_uppersurface=Rcon-2*d_p;
        
%         pos_boundparti=[tmpx_step1(k),tmpz_step1(k)];
%         dist_boundparti_Cencir=sqrt(tmpx_step1(k).^2+tmpz_step1(k).^2);
%         up_surf_p=intersect(k(dist_boundparti_Cencir<Rcir_uppersurface),k(pos_boundparti(:,1)>0.04));
%         up_surf_p_posiaccu=[up_surf_p_posiaccu; tmpx_step1(up_surf_p) tmpz_step1(up_surf_p)];%position accumulation
%         plot(tmpx_step1(up_surf_index(:,targi)),tmpz_step1(up_surf_index(:,targi)),'m-',tmpx_step1,tmpz_step1,'b*');
%         plot(tmpx_step1(up_surf_index),tmpz_step1(up_surf_index),'*');
%         hold on
% toc
        %% ===============obtain velocity and position in step2====================
        velxincell_step2=zeros(size(idincell_step1,1),size(idincell_step1,2));
        velzincell_step2=zeros(size(idincell_step1,1),size(idincell_step1,2));
        posixincell_step2=zeros(size(idincell_step1,1),size(idincell_step1,2));
        posizincell_step2=zeros(size(idincell_step1,1),size(idincell_step1,2));
        for step2i=1:size(idincell_step1,1)
            for step2j=1:size(idincell_step1,2)
                idtmp=cell2mat(idincell_step1(step2i,step2j));
                if ~isempty(idtmp)
                posixincell_step2(step2i,step2j)=mean(idpovel_step2(idtmp,2));
                posizincell_step2(step2i,step2j)=mean(idpovel_step2(idtmp,3));
                velxincell_step2(step2i,step2j)=mean(idpovel_step2(idtmp,4));
                velzincell_step2(step2i,step2j)=mean(idpovel_step2(idtmp,5));
                else
                    continue
                end
            end
        end
        
        %% ==============strain rate for each cell====================
%------------ version 1 -----------------
%         delta_vx=(abs(diff(posixincell_step2,1,2))-abs(diff(posixincell_step1,1,2)))/1/DEMstepinterval/DEMtimestep;
%         delta_vz=(abs(diff(posizincell_step2))-abs(diff(posizincell_step1)))/1/DEMstepinterval/DEMtimestep;
% %         delta_posix=diff(posixincell_step2,1,2)-diff(posixincell_step1,1,2);
% %         delta_posiz=diff(posizincell_step2)-diff(posizincell_step1);
%         delta_posix=abs(diff(posixincell_step1,1,2));
%         delta_posiz=abs(diff(posizincell_step1));
%         strain_rate=0.5*(delta_vx(1:size(delta_vx,2),1:size(delta_vx,2))./delta_posiz(1:size(delta_vx,2),1:size(delta_vx,2))+delta_vz(1:size(delta_vx,2),1:size(delta_vx,2))./delta_posix(1:size(delta_vx,2),1:size(delta_vx,2)));
%         strain_rate(isnan(strain_rate))=0;
% %         strain_rate(strain_rate<-strain_rate_threshold)=-strain_rate_threshold;
% %         strain_rate(strain_rate>strain_rate_threshold)=strain_rate_threshold;
%         strain_rate_accu=strain_rate_accu+strain_rate;
%         strain_rate_cell{esri,1}=strain_rate;
%------------ version 1 end -----------------

%------------ version 2 -----------------
        delta_vx_tmp=(posixincell_step2-posixincell_step1)/1/DEMstepinterval/DEMtimestep;
        delta_vx=delta_vx+delta_vx_tmp(1:Rxgri-2,1:Rxgri-2);
        delta_vz_tmp=(posizincell_step2-posizincell_step1)/1/DEMstepinterval/DEMtimestep;
        delta_vz=delta_vz+delta_vz_tmp(1:Rxgri-2,1:Rxgri-2);
%         delta_posix=diff(posixincell_step2,1,2)-diff(posixincell_step1,1,2);
%         delta_posiz=diff(posizincell_step2)-diff(posizincell_step1);
        delta_posix_tmp=abs(diff(posixincell_step1,1,2));
        delta_posix=delta_posix+delta_posix_tmp(1:Rxgri-2,1:Rxgri-2);
        delta_posiz_tmp=abs(diff(posizincell_step1));
        delta_posiz=delta_posiz+delta_posiz_tmp(1:Rxgri-2,1:Rxgri-2);
%         strain_rate=0.5*(delta_vx(1:Rxgri-2,1:Rxgri-2)./delta_posiz(1:Rxgri-2,1:Rxgri-2)+delta_vz(1:Rxgri-2,1:Rxgri-2)./delta_posix(1:Rxgri-2,1:Rxgri-2));
%         strain_rate(isnan(strain_rate))=0;
%         strain_rate(strain_rate<-strain_rate_threshold)=-strain_rate_threshold;
%         strain_rate(strain_rate>strain_rate_threshold)=strain_rate_threshold;
%         strain_rate_accu=strain_rate_accu+strain_rate;
%         strain_rate_cell{esri,1}=strain_rate;
%------------ version 2 end -----------------

%------------ version 3 -----------------
%         delta_vx_tmp=diff(posixincell_step2)/1/DEMstepinterval/DEMtimestep;
%         delta_vx=delta_vx+delta_vx_tmp(1:Rxgri-2,1:Rxgri-2);
%         delta_vz=delta_vz+diff(posizincell_step2,1,2)/1/DEMstepinterval/DEMtimestep;
% 
%         delta_posix=abs(diff(posixincell_step1,1,2));
%         delta_posiz=abs(diff(posizincell_step1));
% %         strain_rate=0.5*(delta_vx(1:Rxgri-2,1:Rxgri-2)./delta_posiz(1:Rxgri-2,1:Rxgri-2)+delta_vz(1:Rxgri-2,1:Rxgri-2)./delta_posix(1:Rxgri-2,1:Rxgri-2));
% %         strain_rate(isnan(strain_rate))=0;
% %         strain_rate(strain_rate<-strain_rate_threshold)=-strain_rate_threshold;
% %         strain_rate(strain_rate>strain_rate_threshold)=strain_rate_threshold;
% %         strain_rate_accu=strain_rate_accu+strain_rate;
% %         strain_rate_cell{esri,1}=strain_rate;
%------------ version 3 end -----------------
        end
        %% average accumulated data
        delta_vx=delta_vx./accu_average_counts;%delta_vx becomes nan here, it does not affect the results
        delta_vz=delta_vz./accu_average_counts;
        delta_posix=delta_posix./accu_average_counts;
        delta_posiz=delta_posiz./accu_average_counts;
%         strain_rate_accu=strain_rate_accu/accu_steps;
        strain_rate_accu=0.5*(delta_vx(1:Rxgri-3,1:Rxgri-3)./delta_posiz(2:Rxgri-2,2:Rxgri-2)+delta_vz(1:Rxgri-3,1:Rxgri-3)./delta_posix(2:Rxgri-2,2:Rxgri-2));
        %% data to be saved
%         strain_rate_accu
        Xmatrix=Xpf(1:srms,1:srms);
        Zmatrix=Zpf(1:srms,1:srms);
        %% plot
        strain_rate_accu_tmp=strain_rate_accu;
        strain_rate_threshold=100000;
        strain_rate_accu_tmp(strain_rate_accu_tmp<-strain_rate_threshold)=-strain_rate_threshold;
        strain_rate_accu_tmp(strain_rate_accu_tmp>strain_rate_threshold)=strain_rate_threshold;
        figure(2);
        [c,h]=contourf(Xpf(2:srms,2:srms),Zpf(2:srms,2:srms),strain_rate_accu_tmp,200);

%         delta_vx=(diff(posixincell_step2)-diff(posixincell_step1))/2/50000/10^-7;
%         delta_vz=(diff(posizincell_step2)-diff(posizincell_step1))/2/50000/10^-7;
%         delta_posix=diff(posixincell_step2)-diff(posixincell_step1);
%         delta_posiz=diff(posizincell_step2)-diff(posizincell_step1);
%         strain_rate=0.5*(delta_vx+delta_vz);
%         strain_rate(isnan(strain_rate))=0;
%         strain_rate(strain_rate<-0.002)=0;
%         strain_rate(strain_rate>0.002)=0;
%         figure(1);
%         contourf(Xpf(1:end-1,:),Zpf(1:end-1,:),strain_rate);


%         delta_vx=(posixincell_step2-posixincell_step1)/1/50000/10^-7;
%         delta_vz=(posizincell_step2-posizincell_step1)/1/50000/10^-7;
%         delta_posix=(posixincell_step2-posixincell_step1);
%         delta_posiz=(posizincell_step2-posizincell_step1);
%         strain_rate=0.5*(delta_vx./delta_posiz+delta_vz./delta_posix);
%         strain_rate(isnan(strain_rate))=0;
%         strain_rate(strain_rate<-700)=0;
%         strain_rate(strain_rate>700)=0;
%         figure(1);
%         contourf(Xpf(:,:),Zpf(:,:),strain_rate);

%         set(gca,'ColorScale','log');
        set(h, 'edgecolor','none');
        colormap(jet);
        colorbar;
        
        %% ------------extract boundary line position phi=0.225-------------------
%         [boundlinetmp]=contours(Xpf,Zpf,phi_incell,[0.01 0.01]);%use "contours" to not plot the contour
%         boundline=[boundlinetmp(1,2:end);boundlinetmp(2,2:end)]';
%         dist_boundparti_Cencir=sqrt(boundline(:,1).^2+boundline(:,2).^2);
%         up_surf_index=intersect(find(dist_boundparti_Cencir<Rcir_uppersurface),find(boundline(:,1)>-0.1));
%         % % ------------------AOR of differential concecutive points--------
% %         diffangle=atan(diff(boundline(up_surf_index,2))./diff(boundline(up_surf_index,1)));
% %         for i=1:size(diffangle,1)-4
% %         diffangle_ave(i,1)=flip(mean(diffangle(i:i+4,1)));
% %         end
%         
%         %% --------------------plot-------------------
%         figure(2);
%         hold on
%     %     contourf(Xpf,Zpf,phi_incell,50,'LineColor','none');
%         contour(Xpf,Zpf,phi_incell,[0.55 0.55]);
% 
%         contour(Xpf,Zpf,phi_incell,[0.5 0.5]);
%         contour(Xpf,Zpf,phi_incell,[0.4 0.4]);
%         contour(Xpf,Zpf,phi_incell,[0.3 0.3]);
%         contour(Xpf,Zpf,phi_incell,[0.225 0.225]);
%         contour(Xpf,Zpf,phi_incell,[0.1 0.1]);
%         contour(Xpf,Zpf,phi_incell,[0.05 0.05]);
%         contour(Xpf,Zpf,phi_incell,[0.01 0.01]);
%         caxis([0.1 0.7]);
%         legend('\phi=0.55','\phi=0.5','\phi=0.4','\phi=0.3','\phi=0.225','\phi=0.1','\phi=0.05','\phi=0.01');
%         set(gca,'FontSize',15);
%         hold on;
%         plot(boundline(up_surf_index,1),boundline(up_surf_index,2),'ko');
% %         figure(3)
% %         plot(boundline(up_surf_p(1:end-1),1),diffangle)
% %         hold on
% %         plot(boundline(up_surf_p([1:size(diffangle_ave,1)]'+2)),diffangle_ave);
% %         ylabel('radian')
% %         xlabel('distance to left boundary')
% %         xlim([-0.11 0.11]);
%         %% ---------------fitting DAOR-------------------
%         % Set up fittype and options.
%         ft = fittype( 'a*x+b', 'independent', 'x', 'dependent', 'y' );
%         opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%         opts.Display = 'Off';
%         opts.StartPoint = [0.395894484205228 0.936450101590675];
%         fitangle_a=[];
%         fitangle_b=[];
%         fitdatamidx=[];
%         fitdatamidy=[];
%         num_points_to_fit=20;
%         for fiti=1:size(up_surf_index,1)-num_points_to_fit+1
%         [fitresult, gof] = fit(boundline(up_surf_index(fiti:fiti+num_points_to_fit-1),1),boundline(up_surf_index(fiti:fiti+num_points_to_fit-1),2),ft,opts);
%         fitdatamidx(fiti,1)=mean(boundline(up_surf_index(fiti:fiti+num_points_to_fit-1),1));
%         fitdatamidy(fiti,1)=mean(boundline(up_surf_index(fiti:fiti+num_points_to_fit-1),2));
%         fitangle_a(fiti,1)=atan(fitresult.a);%dynamic angle of repose in radius
%         fitangle_b(fiti,1)=fitresult.b;%dynamic angle of repose in degree
%         gof_fitangle(fiti,1)=gof.rsquare;
%         end
%         DAOR(samplei,1)=max(fitangle_a);     
%         DAOR_gof(samplei,1)=gof_fitangle(fitangle_a==max(fitangle_a));
%         figure(2);
% %         hold on;
% %         plot(boundline(up_surf_index(10:end-4),1),boundline(up_surf_index(10:end-4),2),'ko');
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
% %         plot(boundline(up_surf_index([1:end-num_points_to_fit+1]+2),1),fitangle_a);
%         %% --------plot max angle point-------------
%         coord_maxanglex=fitdatamidx(fitangle_a==max(fitangle_a));
%         coord_maxangley=fitdatamidy(fitangle_a==max(fitangle_a));
%         plot(coord_maxanglex,coord_maxangley,'go','MarkerSize',10,'LineWidth',1);
%         
%         %% ----------plot perpendicular line of the fitting line
%         currslope=tan(max(fitangle_a));
%         perpslope=-1/currslope;
%         perp_b=fitdatamidy(fitangle_a==max(fitangle_a))-(perpslope*fitdatamidx(fitangle_a==max(fitangle_a)));
% %         figure(2); hold on;
% %         yft_perp=perpslope*xft+perp_b;
% %         plot(xft,yft_perp,'g-','LineWidth',1);
%         %% ----------plot intersection between perpendicular line and circle
% %         figure(2);hold on;
%         [interx,intery] = linecirc(perpslope,perp_b,Centercir(1),Centercir(2),Rcon);
% %         plot(interx(1),intery(1),'y*','MarkerSize',10,'LineWidth',1);
%         dist_OI=sqrt((coord_maxanglex-interx(1)).^2+(coord_maxangley-intery(1)).^2);
%         %% -----plot rotated data according to AOR point---------------
% %         coordtmp=[reshape(xextacu,[],1) (reshape(yextacu,[],1)-1080)*-1];% (data-1080)*-1 means a reverse in y direction
%         rota_datacoord= rotation([tmpx_step1 tmpz_step1],[coord_maxanglex,coord_maxangley],-DAOR(samplei,1)/pi*180);
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
%         %% ---------plot sampling region---------
%         figure(7);hold on;
%         atomrange=d_p*5;%the assumed pixel for one particle
%         S1(1,1)=coord_maxanglex+atomrange;S1(1,2)=coord_maxangley;
%         S2(1,1)=coord_maxanglex-atomrange;S2(1,2)=coord_maxangley;
%         S3(1,1)=coord_maxanglex+atomrange;S3(1,2)=coord_maxangley-dist_OI;
%         S4(1,1)=coord_maxanglex-atomrange;S4(1,2)=coord_maxangley-dist_OI;
%         S_coord=[S1;S2;S4;S3];
%         plot(S_coord(:,1),S_coord(:,2),'o','MarkerSize',10,'LineWidth',1,'MarkerFaceColor',[0.9290 0.6940 0.1250]);
%         plot(polyshape(S_coord(:,1),S_coord(:,2)));
%         %% -----extract velocity profile--------------
%         
%         vel_dircorr=ones(size(tmpvel_step1,1),1);%velocity direction correction
%         vel_dircorr(tmpvel_step1(:,1)>0)=-1;
%         vel_totacu=sqrt(tmpvel_step1(:,1).^2+tmpvel_step1(:,3).^2)/sqrt(9.81*d_p).*vel_dircorr;
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
%         %% extract velocity field in gravity frame
%         topbound_particles=boundline(up_surf_index,:);       
%         % % ---x and z boundary---------
%         lb=min(topbound_particles(:,1));%left boundary
%         rb=max(topbound_particles(:,1));
% 
%         % % ---- subdivide particles to column domain -------
%         ww=d_p*1;%window width
%         subd_domains=ceil(abs(lb-rb)/ww)';
%         drlbound=linspace(lb,rb,subd_domains)';%domain right and lef boundary arrays
%         drlbound_midpoint=drlbound(1:end-1)+ww/2;
%         sol_liq_interface=zeros(subd_domains,2);
%         cust_grad_color=colormap(copper(size(drlbound,1)-1));
%         % % --- plot velocity countour ----
% %         tmpx_matrix=repmat(tmpx_step1,1,size(tmpz_step1,1));
% %         tmpz_matrix=repmat(tmpz_step1',size(tmpx_step1,1),1);
% %         vel_matrix=repmat(velscaler,1,size(tmpz_step1,1));
% %         contour(tmpx_matrix,tmpz_matrix,vel_matrix);
%         figure(10);
% %         velscaler=abs(tmpx_step1);
%         velscaler=veltotalscaler;
% %         velplot=scatter(tmpx_step1,tmpz_step1,[],velscaler);
%         veltcontour=contourf(Xpf,Zpf,veltincell);
%         for dblri=1:size(drlbound,1)-1
%             subrldx=[];
%             subrldz=[];
%             subrlvel=[];
%             subdx_index=[];
%             subdx_index=intersect(find(tmpx_step1>drlbound(dblri)),find(tmpx_step1<drlbound(dblri+1)));
%             subrldx=tmpx_step1(subdx_index);
%             subrldz=tmpz_step1(subdx_index);
%             subrlvel=velscaler(subdx_index);
%      
%             % % ---- subdivide column domains to layerly domains to obtain total velocity layerly -------
%             bb=min(subrldz);%bottom boundary in each column
%             tb=max(subrldz);  
%             wd=d_p*1;%window depth
%             wnum=ceil(abs(bb-tb)/wd)';%total window number
%             dtbbound=linspace(bb,tb,wnum)';
%             dtbbound_midpoint=dtbbound(1:end-1)+wd/2;
%             subtbvel=zeros(size(dtbbound,1)-1,1);
%              for dbtbi=1:size(dtbbound,1)-1
%                  subdz_index=[];
%                  subdz_index=intersect(find(subrldz>dtbbound(dbtbi)),find(subrldz<dtbbound(dbtbi+1)));
%                  subtbvel(dbtbi,1)=mean(subrlvel(subdz_index));
%              end
%              figure(10);hold on;
%              max_subtbvel=max(subtbvel);
%              min_subtbvel=min(subtbvel);
%              sol_liq_interface(dblri,:)=[max(subrldz) dtbbound_midpoint(subtbvel==min_subtbvel)];
%              plot(repmat(drlbound_midpoint(dblri,1),1,2),sol_liq_interface(dblri,:),'*','MarkerSize',10,'color',cust_grad_color(dblri,:));
%              figure(11)%single column domain velocity profile
%              plot(subtbvel,dtbbound(1:end-1),'color',cust_grad_color(dblri,:));hold on;
%         end
%         
    %% ==============dynamic Angle of Repose: fit raw surface data=================
%     % Set up fittype and options.
%     ft = fittype( 'a*x+b', 'independent', 'x', 'dependent', 'y' );
%     opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%     opts.Display = 'Off';
%     opts.StartPoint = [0.395894484205228 0.936450101590675];
%     % Fit model to data.
%     % [fitresult, gof] = fit( tmpx_step1(up_surf_p),tmpz_step1(up_surf_p), ft, opts );
%     [fitresult, gof] = fit( up_surf_p_posiaccu(:,1),up_surf_p_posiaccu(:,2), ft, opts );
%     DAOR_raw=atan(fitresult.a)/pi*180;%dynamic angle of repose in degree
%     gof_raw=gof.rsquare;
% 
%     % plot(up_surf_p_posiaccu(:,1),up_surf_p_posiaccu(:,2),'ko');
%     % xlim([0 0.12]);
%     % ylim([-0.1 0.05]);
%     % hold on
%     % xft=[-0.1:0.01:0.1];
%     % yft=fitresult.a*xft+fitresult.b;
%     % plot(xft,yft,'r-','LineWidth',2);
%     DAORt_fitraw(targi,1)=DAOR_raw;
%     goft_raw(targi,1)=gof_raw;


    %% =================dynamic Angle of Repose: probability density======================
%     xgrid=linspace(-0.12,0.12,50);
%     ygrid=linspace(-0.12,0.12,50);
%     [x1,y1] = meshgrid(xgrid, ygrid);
%     % Perform kernel density estimate
%     % [x y] is actual data, xi is the desired grid points to evaluate
%     % f is an estimate of the density, ep(:,1) is the X location, ep(:,2) is the y location
%     xi = [x1(:) y1(:)];
%     [f,ep]=ksdensity([tmpx_step1(:,1) tmpz_step1(:,1)],xi,'Bandwidth',0.01); % remove the outputs to see a 3D plot of the distribution
%     f=f/max(f);
%     % format data in matrix for contourf and plot
%     X = reshape(ep(:,1),length(xgrid),length(ygrid));
%     Y = reshape(ep(:,2),length(xgrid),length(ygrid));
%     Z = reshape(f,length(xgrid),length(ygrid));
%     % ----------------plot----------------
%     %     figure(1)
% %     hold on
% %     X(Z<=density_threshold)=nan;
% %     Y(Z<=density_threshold)=nan;
% %     Z(Z<=density_threshold)=nan;
%     figure(1)
%     fig=contourf(X,Y,Z,50,'LineColor','none');
%     caxis([0 1]);
%     xlim([-0.12 0.12]);
%     ylim([-0.12 0.12]);
% 
%     % ------------------for fitting-------------
%     Xtmp=X(Z>density_threshold);
%     Ytmp=Y(Z>density_threshold);
%     Ztmp=Z(Z>density_threshold)*100;
%     Ztmp=round(Ztmp);
%     positmp=[];
%     for Xi=1:size(Xtmp,1)
%         positmp=[positmp;repmat([Xtmp(Xi,1) Ytmp(Xi,1)],Ztmp(Xi,1),1)];
%     end
% 
%     % ------------------Set up fittype and options-------------
%     ft = fittype( 'a*x+b', 'independent', 'x', 'dependent', 'y' );
%     opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%     opts.Display = 'Off';
%     opts.StartPoint = [0.395894484205228 0.936450101590675];
%     % Fit model to data.
%     % [fitresult, gof] = fit( tmpx_step1(up_surf_p),tmpz_step1(up_surf_p), ft, opts );
%     [fitresult, gof] = fit( positmp(:,1),positmp(:,2), ft, opts );
%     DAOR_den=atan(fitresult.a)/pi*180;%dynamic angle of repose in degree
%     gof_den=gof.rsquare;
%     % hold on
%     % xft=[-0.1:0.01:0.1];
%     % yft=fitresult.a*xft+fitresult.b;
%     % plot(xft,yft,'r-','LineWidth',2);
%     DAORt_fitdensity(targi,1)=DAOR_den;
%     goft_den(targi,1)=gof_den;
    
    end
    %% ===========save results to mat file=============
%    cd(datafile);
%    datafilename=['model' num2str(modelnum(modeli,1)) '_strainrate_accu100_gammainf.mat'];
%    save(datafilename,'strain_rate_accu');
%    save(datafilename,'Xmatrix','-append');
%    save(datafilename,'Zmatrix','-append');
% %    clear strain_rate_accu Xmatrix Zmatrix
   
       %% ===========save figure=============
%    cd('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Paper figures\strain_rate_figure');
%    saveas(gca,['strainrate_accu100_gammainf_model' num2str(modelnum(modeli,1)) '.png']);
%    close all;
   
end

    
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