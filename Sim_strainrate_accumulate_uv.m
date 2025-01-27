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
% modelnum=[81:84]';
% modelnum=[49:52]';
% modelnum=[70]';
% modelnum=[103]';
% modelnum=[104,105]';
% modelnum=[106:111]';
% modelnum=[119,121:124,126:131]';
% modelnum=[134]';
% modelnum=[141:150]';
modelnum=[120,125,140,151,152]';
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
omega=keypara(modelnum(modeli,1),5);
d_p=keypara(modelnum(modeli,1),7);
V_p=4/3*pi*(d_p/2)^3;
Rp=d_p/2;
Rcon=0.1;%container radius
Centercir=[0,0];
container_d=0.03;
phi_sample_meshsize=1*d_p;
%% sampling range
e_sample_r=2;
sample_range=[500:e_sample_r:Rtxt-10]';
% strain_rate_threshold=1000;
accu_steps=size(sample_range,1);
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
    for samplei=1
%         :size(sample_range,1) %the start timestep of accumulation
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
        accu_average_counts=zeros(Rxgri-1,Rxgri-1);%count how many frames the cell has particles
        velxincell_step1_total=zeros(size(xgridpf,1)-1);
        velzincell_step1_total=zeros(size(xgridpf,1)-1);
        %%
        tmpx_step1=[];
        tmpz_step1=[];
        tmpvel_step1=[];
%         tmpid_step1=[];
        for esri=1:accu_steps
            disp(['model' num2str(modelnum(modeli,1)) '_' num2str(round(pocesstime(esri,1),4))]);
            
            if modelnum(modeli,1)<118
                revolution(esri,1)=sample_range(samplei+esri-1)*DEMstepinterval*DEMtimestep*omega/360;% for model 25-30 sample_range-24
                [atomnum,~,~]=xlsread( ...
                [Foldersdir_atominfo '\' txt_read_atominfo(sample_range(samplei+esri-1)).name],1);
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
                    txt_read_atominfo(sample_range(samplei+esri-1)).name,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','headerlines',9);
                vel=sqrt(vx.^2+vy.^2+vz.^2);
            end
            tmpx_step1=[tmpx_step1;x];%accumultive x coordinates from several steps
            tmpz_step1=[tmpz_step1;z];
            tmpvel_step1=[tmpvel_step1;[vx,vy,vz]];
%             tmpid_step1=[tmpid_step1;idest];

        end
        %% convexhull method
        % [k,av]=convhull(tmpx_step1,tmpz_step1,'simplify',true);
        % plot(tmpx_step1(k),tmpz_step1(k),'r-',tmpx_step1,tmpz_step1,'b*')

        %% boundary method
        % k=boundary(tmpx_step1,tmpz_step1,0.8);
        % plot(tmpx_step1(k),tmpz_step1(k),'r-',tmpx_step1,tmpz_step1,'b*')

        %% top surface method
%         % boundL=min(tmpx_step1);
%         % boundR=max(tmpx_step1);
%         % LRrange=boundR-boundL;
%         % divLR=LRrange/10;

%         veltotalscaler=sqrt(sum(tmpvel_step1.^2,2));
        %% ============= plot 2D particles ===============
%         figure(2)
%         particle_plot=scatter(tmpx_step1,tmpz_step1);
        %% =================allocate particles to each cell and calculate Packing fraction and mean velocity======================
        tic
        Nincell=zeros(size(xgridpf,1)-1);
        veltincell=zeros(size(xgridpf,1)-1);%scaler velocity in cell
        meanvelxincell=zeros(size(xgridpf,1)-1);
        meanvelzincell=zeros(size(xgridpf,1)-1);
        posixincell_step1=zeros(size(xgridpf,1)-1);
        posizincell_step1=zeros(size(xgridpf,1)-1);
        idincell_step1=cell(size(xgridpf,1)-1,1);% mark id of particles in cell
        for xgi=1:Rxgri-1
                for zgi=1:Rzgri-1
                    xwindowa=find(tmpx_step1<xgridpf(xgi+1,1));
                    xwindowb=find(tmpx_step1>xgridpf(xgi,1));
                    zwindowa=find(tmpz_step1<zgridpf(zgi+1,1));
                    zwindowb=find(tmpz_step1>zgridpf(zgi,1));
                    xwindow=xwindowa(ismembc(xwindowa,xwindowb));
                    zwindow=zwindowa(ismembc(zwindowa,zwindowb));
%                     inters_id=intersect(xrange,zrange);
                    inters_rownum=xwindow(ismembc(xwindow,zwindow));
                    Nincell(xgi,zgi)=size(inters_rownum,1);
                    if Nincell(xgi,zgi)>(accu_steps*0.05)
%                         accu_average_counts(xgi,zgi)=accu_average_counts(xgi,zgi)+1;
%                         veltincell(xgi,zgi)=mean(veltotalscaler(inters_rownum,1));
%                         idincell_step1{xgi,zgi}=tmpid_step1(inters_rownum);
%                         posixincell_step1(xgi,zgi)=mean(tmpx_step1(inters_rownum));
%                         posizincell_step1(xgi,zgi)=mean(tmpz_step1(inters_rownum));
                        meanvelxincell(xgi,zgi)=mean(tmpvel_step1(inters_rownum,1));
                        meanvelzincell(xgi,zgi)=mean(tmpvel_step1(inters_rownum,3));
                        tmpx_step1(inters_rownum,:)=[];
                        tmpz_step1(inters_rownum,:)=[];
                        tmpvel_step1(inters_rownum,:)=[];
                    else
%                         idincell_step1{xgi,zgi}=[];
                    end
                end
        end
        toc
        end
    %% ===========save results=============
       cd(datafile);
%        datafilename=['model' num2str(modelnum(modeli,1)) '_velmatrix' num2str(accu_steps) 'frames.mat'];% add "accu_steps" into names
       datafilename=['model' num2str(modelnum(modeli,1)) '_velmatrix_frames.mat'];
       save(datafilename,'Xpf');
       save(datafilename,'Zpf','-append');
    %    save(datafilename,'DAORt_fitraw','-append');
       save(datafilename,'meanvelxincell','-append');
       save(datafilename,'meanvelzincell','-append');
       clear Xpf Zpf meanvelxincell meanvelzincell
       
    end
    
