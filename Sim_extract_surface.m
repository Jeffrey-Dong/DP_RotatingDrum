clc;
clear;
%%
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\dataanalysis\');
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\2d_packing\data_analysis_code\functions\');
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis\function');
[num,txt,raw]=xlsread(['C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\parametric_study.xlsx'],1);
keypara=num(4:end,:);
density_threshold=0.6;
% modelnum=[31,33,35,37,39,41]';
% modelnum=[33,53,54,59:69]';
% modelnum=[77:79]';
% modelnum=[87:97,99:102]';
% modelnum=[81:84]';
% modelnum=[70]';
% modelnum=[103]';
% modelnum=[104,105]';
% modelnum=[106:111]';
% modelnum=[112:114]';
% modelnum=[90:97,99:114,116]';
% modelnum=[103,112,119,129:130]';
% modelnum=[122:123,124,126:128]';
% modelnum=[141:150]';
% modelnum=[139]';
modelnum=[140,151,152]';%stop at 140
%%
for modeli=1:size(modelnum,1)
%% prepare data ----method1---------
datafile=(['E:\RoDrtest\RoDr_alphastart0-1\model' num2str(modelnum(modeli,1))]);
% load([datafile '\model' num2str(modelnum(modeli,1)) '_posi.mat']);
% load([datafile '\model' num2str(modelnum(modeli,1)) '_velocity.mat']);
% load([datafile '\model' num2str(modelnum(modeli,1)) '_id.mat']);
% [Rtxt,Ctxt]=size(centroids);
%% prepare data ----method2---------
Foldersdir_pairwise = (['D:\RoDrtest\RoDr_model' num2str(modelnum(modeli,1)) '_pairwise']);
subFolders_pairwise = dir(Foldersdir_pairwise);
subFolders_pairwise=subFolders_pairwise(3:end);
if modelnum(modeli,1)<118
    Foldersdir_atominfo = (['E:\RoDrtest\RoDr_alphastart0-1\model' num2str(modelnum(modeli,1))]);
    txt_read_atominfo=dir([Foldersdir_atominfo '\*.csv']);
    [txt_read_atominfo]=sortnamebysequence(txt_read_atominfo);
    [Rtxt,Ctxt]=size(txt_read_atominfo);
else
    Foldersdir_atominfo = (['D:\RoDrtest\model' num2str(modelnum(modeli,1)) '_atominfo']);
    txt_read_atominfo=dir([Foldersdir_atominfo '\*.txt']);
    [txt_read_atominfo]=sortnamebysequence(txt_read_atominfo);
    [Rtxt,Ctxt]=size(txt_read_atominfo);
end
liquid_content=keypara(modelnum(modeli,1),6);
%%

DEMtimestep=10^-7;
DEMstepinterval=50000;
DAORt_den=[];
DAORt_ras=[];
goft_den=[];
goft_raw=[];
revolution=[];
DAOR=[];
DAOR_gof=[];
up_surf_points={};
omega=keypara(modelnum(modeli,1),5);
d_p=keypara(modelnum(modeli,1),7);
V_p=4/3*pi*(d_p/2)^3;
Rp=d_p/2;
if modelnum(modeli,1)<118
    Rcon=0.1;%container radius
else
    Rcon=0.093;%container radius
end
Centercir=[0,0];
container_d=0.03;
phi_sample_meshsize=1*d_p;
samplenumber=30;

%% vary startstep based on wet and dry condition
if liquid_content>0
    startstep=500;
else
    startstep=300;
end

sampleinterval=round((Rtxt-startstep)/samplenumber);

%% 
sample_range=[startstep:sampleinterval:Rtxt]';
%% ======= prepare the mesh region to extract data and do analysis =========
meshbound=[-0.12 0.12];
xgridpf=linspace(meshbound(1),meshbound(2),(meshbound(2)-meshbound(1))/phi_sample_meshsize)';
zgridpf=linspace(meshbound(1),meshbound(2),(meshbound(2)-meshbound(1))/phi_sample_meshsize)';
[Rxgri,Cxgri]=size(xgridpf);
[Rzgri,Czgri]=size(zgridpf);
Xpf=repmat(xgridpf(1:Rxgri-1,1),1,Rxgri-1);
Zpf=repmat(zgridpf(1:Rxgri-1,1)',Rxgri-1,1);

%%
e_sample_r=size(sample_range,1)-1;
    for samplei=1:e_sample_r
    disp(['model' num2str(modelnum(modeli,1)) '_step' num2str(samplei/e_sample_r)]);
    tmpx=[];
    tmpz=[];
    tmpvel=[];
    
        %% prepare data for timestep accumulation
%         tic
    for esri=1:7 %each sampling range
        revolution(esri,1)=sample_range(samplei)*DEMstepinterval*DEMtimestep*omega/360;% for model 25-30 sample_range-24
        if modelnum(modeli,1)<118
            [atomnum,~,~]=xlsread( ...
            [Foldersdir_atominfo '\' txt_read_atominfo(sample_range(samplei)-4+esri).name],1);
            idest=atomnum(:,1);
            x=atomnum(:,18);
            y=atomnum(:,19);
            z=atomnum(:,20);
            radius_est=atomnum(:,12);
            vx=atomnum(:,6);
            vy=atomnum(:,7);
            vz=atomnum(:,8);
            VliqonVatom=atomnum(:,17);
            vel=sqrt(vx.^2+vy.^2+vz.^2);
        else
            cd([Foldersdir_atominfo]);
            %     :size(samplei,2) %file number
            % [id,type,x,y,z,vx,vy,vz,fx,fy,fz,radius,liquidcontent,fpenx,vorovolume,omegax,omegay,omegaz]=textread(txt_read(i).name,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','headerlines',9);
            [idest,~,x,y,z,~,~,~,vx,vy,vz,~,~,~,~,~,~,radius_est,VliqonVatom]=textread( ...
                txt_read_atominfo(sample_range(samplei)-4+esri).name,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','headerlines',9);
            vel=sqrt(vx.^2+vy.^2+vz.^2);
        end

%         atom_posi_eachstep=cell2mat(centroids(sample_range(samplei+esri-1),1));
%         atom_vel_eachstep=cell2mat(velocity(sample_range(samplei+esri-1),1));

        tmpx=[tmpx;x];%accumultive x coordinates from several steps
        tmpz=[tmpz;z];
%         atom_vel_eachstep=cell2mat(velocity(sample_range(samplei+esri-1),1));
%         tmpvel=[tmpvel;atom_vel_eachstep];

        %% convexhull method
        % [k,av]=convhull(tmpx,tmpz,'simplify',true);
        % plot(tmpx(k),tmpz(k),'r-',tmpx,tmpz,'b*')

        %% boundary method
        % k=boundary(tmpx,tmpz,0.8);
        % plot(tmpx(k),tmpz(k),'r-',tmpx,tmpz,'b*')

        %% top surface method
        % boundL=min(tmpx);
        % boundR=max(tmpx);
        % LRrange=boundR-boundL;
        % divLR=LRrange/10;

        end
%         veltotalscaler=sqrt(sum(tmpvel.^2,2));
%         toc
        %% ============= plot 2D particles ===============
%         figure(2)
%         particle_plot=scatter(tmpx,tmpz);
        
%% =================dynamic Angle of Repose: packing fraction and total velocity======================
% tic

        Nincell=zeros(size(xgridpf,1)-1);
        veltincel=zeros(size(xgridpf,1)-1);%scaler velocity in cell
        for xgi=1:Rxgri-1
            for zgi=1:Rzgri-1
                xrange=[];
                zrange=[];
                xrangea=find(tmpx<xgridpf(xgi+1,1));
                xrangeb=find(tmpx>xgridpf(xgi,1));
                xrange=xrangea(ismember(xrangea,xrangeb));
                zrangea=find(tmpz<zgridpf(zgi+1,1));
                zrangeb=find(tmpz>zgridpf(zgi,1));
                zrange=zrangea(ismember(zrangea,zrangeb));
                Nincell(xgi,zgi)=size(xrange(ismember(xrange,zrange)),1);
%                 veltincel(xgi,zgi)=mean(veltotalscaler(intersect(xrange,zrange),1));
            end
        end

        phi_incell=V_p.*Nincell/(abs(xgridpf(2)-xgridpf(1))*abs(zgridpf(2)-zgridpf(1))*container_d*esri);
        % the circle to extract upper surface particles
        Rcir_uppersurface=Rcon-2*d_p;
        

% toc
        %% ------------extract boundary line position phi=0.225-------------------
        [boundlinetmp]=contours(Xpf,Zpf,phi_incell,[0.05 0.05]);
        boundlinetmp=boundlinetmp';
        circle_index=find(boundlinetmp(:,2)==max(boundlinetmp(:,2)));
        boundlinetmp1=boundlinetmp(circle_index+1:circle_index+max(boundlinetmp(:,2)),:);
%         figure(3);hold on;scatter(boundlinetmp1(:,1),boundlinetmp1(:,2),[],'r.');
        dist_boundparti_Cencir=sqrt(boundlinetmp1(:,1).^2+boundlinetmp1(:,2).^2);
        up_surf_index=find(dist_boundparti_Cencir<Rcir_uppersurface);
        up_surf_points{samplei,1}=boundlinetmp1(up_surf_index,:);
        % % ------------------AOR of differential concecutive points--------
%         diffangle=atan(diff(boundline(up_surf_index,2))./diff(boundline(up_surf_index,1)));
%         for i=1:size(diffangle,1)-4
%         diffangle_ave(i,1)=flip(mean(diffangle(i:i+4,1)));
%         end
        
           
    end
    % % ----------------save up surface points as .mat file------------
cd(datafile);
datafilename=['model' num2str(modelnum(modeli,1)) '_up_surface_points'];
save(datafilename,'up_surf_points');
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