clc;
clear;

addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\dataanalysis\');
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\2d_packing\data_analysis_code\functions\');
[num,txt,raw]=xlsread(['C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\parametric_study.xlsx'],1);
keypara=num(4:end,:);
density_threshold=0.6;

pldist_threshold_coe_dry=1.00000001;%this is only for creating pairlist using pair distance method when dry condition is found with 0 liquid content

clustermethod="angle";
% angle_threshold_array=[pi/3,pi/4,pi/6,pi/9]';
angle_threshold_array=[pi/3]';
% clustermethod="dotvel";
dotvel_coe=3;
dp=0.002;
% modelnum=[31,33,35,37,39,41]';
% modelnum=[49:52]';
% modelnum=[99:102]';
% modelnum=[97,99:102]';
% modelnum=[86:94]';
% modelnum=[81:84,77]';
% modelnum=[103:111]';
% modelnum=[113:114]';
% modelnum=[119,121:123,124,126:128]';
% modelnum=[119,121:123]';
% modelnum=[141:150]';
modelnum=[119,120,134,121:123,124,125,151,126:128,139,140,152,141:143]';
% modelnum=[120,134,125,140,151,152]';
% modelnum=[119,124,139]';
%%
for angi=1:size(angle_threshold_array,1)
    angle_threshold=angle_threshold_array(angi,1);
%% model number
for modeli=1:size(modelnum,1)
%% clear to save memory when jump to the next case
clearvars -except dp modeli modelnum num txt raw keypara density angi angle_threshold angle_threshold_array pldist_threshold_coe_dry clustermethod dotvel_coe;

%% prepare data
% tic
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
omega=keypara(modelnum(modeli,1),5);
d_p=keypara(modelnum(modeli,1),7);
liquid_content=keypara(modelnum(modeli,1),6);
% toc
%% load strain rate boundary
boundary_info=load(['E:\RoDrtest\RoDr_alphastart0-1\model' num2str(modelnum(modeli,1)) ...
    '\model' num2str(modelnum(modeli,1)) '_strrbound_strrcoe7_surface.mat']);
full_bound_filter=[boundary_info.strr_botbound; ...
                   boundary_info.strr_botbound(end,1)+5*dp,boundary_info.strr_botbound(end,2); ...
                   boundary_info.strr_botbound(end,1)+5*dp,boundary_info.strr_botbound(end,2)+10*dp; ...
                   boundary_info.strr_botbound(1,1)-5*dp,boundary_info.strr_botbound(end,2)+10*dp; ...
                   boundary_info.strr_botbound(1,1)-5*dp,boundary_info.strr_botbound(1,2); ...
                   boundary_info.strr_botbound(1,1),boundary_info.strr_botbound(1,2)
                   ];
% figure(2)
% plot(boundary_info.strr_botbound(:,1),boundary_info.strr_botbound(:,2),'k-')
% hold on
% plot(boundary_info.strr_topbound(:,1),boundary_info.strr_topbound(:,2),'o-')
% plot(boundary_info.freesurfacepoints(:,1),boundary_info.freesurfacepoints(:,2),'r-')
% plot(full_bound_filter(:,1),full_bound_filter(:,2),'*-');
%% atom infor from mat reading
% load([atominfo_datafile '\model' num2str(modelnum(modeli,1)) '_posi.mat']);
% load([atominfo_datafile '\model' num2str(modelnum(modeli,1)) '_velocity.mat']);
% load([atominfo_datafile '\model' num2str(modelnum(modeli,1)) '_id.mat']);
% load([atominfo_datafile '\model' num2str(modelnum(modeli,1)) '_radius.mat']);

%% prepare folder info for pairinfor reading from txt reading
% txt_read_pairwise=dir([Foldersdir_pairwise '\*.txt']);
% [txt_read_pairwise]=sortnamebysequence(txt_read_pairwise);
% [Rtxt_pair,Ctxt_pair]=size(txt_read_pairwise);

%% prepair stepnumber for reading pair.txt data
% tic
for sfipair=1:size(subFolders_pairwise,1)
    subFolders_names{sfipair,1}=string(subFolders_pairwise(sfipair).name);
     namesplit=split(subFolders_names{sfipair,1},'.');
     if namesplit(2,1)=="vtk"
         break
     else
         modelNo=str2num(namesplit(1,1));
         stepnumber(sfipair,1)=modelNo;
     end
end
     stepnumber=sort(stepnumber);

% toc
%% Prepare container for atom and contact pair info
% sample_step=[300:20:Rtxt_atom]';
% if Rtxt>=650
%     sample_step=[400:10:650]';
% else
    sample_step=[500:10:Rtxt]';
% end
atomid_clusterid=[];
cluster_id_size=[];
atomid_clustersize=[];
cluster_x=[];
max_cluster_y=[];
cluster_z=[];
cluster_count=0;
cluster_binrecord=[];
cluster_sizerecord=[];
pair_vel_flowzone_pairmark_record=[];
id_centroids_radius_vxyz_record=[];
dbs_cluster_size_record=[];
dbs_cluster_id_record=[];

%% clear pair list 

%%
for i=1:size(sample_step,1)
    disp(['model' num2str(modelnum(modeli,1)) '_step' num2str(sample_step(i,1))])
    %% atom information from csv reading
%     tic
    if modelnum(modeli,1)<118
        revolution(esri,1)=sample_range(samplei+esri-1)*DEMstepinterval*DEMtimestep*omega/360;% for model 25-30 sample_range-24
        [atomnum,~,~]=xlsread( ...
        [Foldersdir_atominfo '\' txt_read_atominfo(sample_step(i)).name],1);
        id=atomnum(:,1);
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
        % id type x y z ix iy iz vx vy vz fx fy fz omegax omegay omegaz radius f_surfaceLiquidContent[0]
        [id,~,x,y,z,~,~,~,vx,vy,vz,fx,fy,fz,~,~,~,radius_est,VliqonVatom]=textread( ...
            txt_read_atominfo(sample_step(i)).name,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','headerlines',9);
        vel=sqrt(vx.^2+vy.^2+vz.^2);
    end
    id_centroids_radius_vxyz=[id,x,y,z,radius_est,vx,vy,vz,fx,fy,fz];
    id_centroids_radius_vxyz=sortrows(id_centroids_radius_vxyz,1);
    centroids_est=id_centroids_radius_vxyz(:,2:4);%each step
    id_est=id_centroids_radius_vxyz(:,1);%each step
    vxyz_est=id_centroids_radius_vxyz(:,6:8);
    fxyz_est=id_centroids_radius_vxyz(:,9:11);
    ftotal=sqrt(fxyz_est(:,1).^2+fxyz_est(:,2).^2+fxyz_est(:,3).^2);
    id_centroids_radius_vxyz_record{i,1}=id_centroids_radius_vxyz;
%     toc
    %% pair information: create pairlist based on pair distance
    tic
    index_abovestrrbound= inpolygon(centroids_est(:,1),centroids_est(:,3),full_bound_filter(:,1),full_bound_filter(:,2));
    centroids_est_strrfilted=centroids_est(index_abovestrrbound,:);
    id_est_strrfilted=id_est(index_abovestrrbound,:);
    centroids_est_filtered=centroids_est_strrfilted;
    id_est_filtered=id_est_strrfilted;
%     centroids_est_filterted=centroids_est(ftotal>=mean(ftotal));
%     id_est_filtered=id_est(ftotal>=mean(ftotal));
%     centroids_est_filterted=centroids_est;
%     id_est_filtered=id_est;
    if liquid_content<0
        pairlistmethod='pldist';
        disp('prepare pairlist based on parir distance')
        %% ---------------method 1 using getpairlist function-------not efficient------
%         path(path,'C:\Users\jeffr\Dropbox (Sydney Uni)\projects\2d_packing\data_analysis_code');
%             if i==1 %this if satement avoid repeating creating pairlist since total pairlist reamain constant among each case
%                 pl=getpairlist(size(idest,1));
%             else
%                 %do nothing
%             end
%     
%         % % build the distance list
%     %     for pli=1:length(pl)
%     %         pl(pli,3)=norm([centroids(pl(pli,1),1),centroids(pl(pli,1),2)]-[centroids(pl(pli,2),1),centroids(pl(pli,2),2)]);
%     %     end
%         pldist=pdist(centroids_est)';
%         pl_tmp1=pl;
%         pl_tmp1(:,3)=pldist;% pair distance

        %% ------------------method 2 filtering neighbour first then assign pair distance-----------
%         tic
        pl=[];
%         centroids_est_vxzless0=centroids_est(vxyz_est(:,1)<0 & vxyz_est(:,3)<0,:);
%         id_est_vxzless0=id_est(vxyz_est(:,1)<0 & vxyz_est(:,3)<0,:);
%         centroids_est_vxzless0=centroids_est;
%         id_est_vxzless0=id_est;
        for pli=1:size(centroids_est_filtered,1)
            xbound=[centroids_est_filtered(pli,1)-1.05*d_p centroids_est_filtered(pli,1)+1.05*d_p];
            ybound=[centroids_est_filtered(pli,2)-1.05*d_p centroids_est_filtered(pli,2)+1.05*d_p];
            zbound=[centroids_est_filtered(pli,3)-1.05*d_p centroids_est_filtered(pli,3)+1.05*d_p];
            id_inthisloop=id_est_filtered(pli+1:end,1);
            index_x1=id_inthisloop(centroids_est_filtered(pli+1:end,1)>xbound(1));
            index_x2=id_inthisloop(centroids_est_filtered(pli+1:end,1)<xbound(2));
            index_y1=id_inthisloop(centroids_est_filtered(pli+1:end,2)>ybound(1));
            index_y2=id_inthisloop(centroids_est_filtered(pli+1:end,2)<ybound(2));
            index_z1=id_inthisloop(centroids_est_filtered(pli+1:end,3)>zbound(1));
            index_z2=id_inthisloop(centroids_est_filtered(pli+1:end,3)<zbound(2));
            index_x=index_x1(ismember(index_x1,index_x2));
            index_y=index_y1(ismember(index_y1,index_y2));
            index_z=index_z1(ismember(index_z1,index_z2));
            index_xy=index_x(ismember(index_x,index_y));
            index_xyz=index_xy(ismember(index_xy,index_z));
            pl=[pl;[ones(size(index_xyz,1),1)*id_est_filtered(pli,1), index_xyz]];
        end
        pl_tmp1=pl;
        pl_tmp1(:,3)=sqrt(sum((centroids_est(pl(:,1),:)-centroids_est(pl(:,2),:)).^2,2));% pair distance
%         toc
        pl_tmp1(:,4)=id_centroids_radius_vxyz(pl_tmp1(:,1),5)+id_centroids_radius_vxyz(pl_tmp1(:,2),5);% r1+r2
        %%
        pl_tmp1(pl_tmp1(:,3)>(pl_tmp1(:,4)*pldist_threshold_coe_dry),:)=[];%remove large gap distance particles based on coe*(r1+r2) wet condition
    %     pl_tmp1(pl_tmp1(:,3)>(pl_tmp1(:,4)*1.02),:)=[];%remove large gap distance particles based on coe*(r1+r2) dry condition
        id1=pl_tmp1(:,1);
        id2=pl_tmp1(:,2);
     
%% pair information: create pairlist based on pair force
    else
        pairlistmethod='plforce';
        disp('prepare pairlist based on parir force')
    %     cd([Foldersdir_pairwise]);
        [~,x1,y1,z1,x2,y2,z2,id1,id2,~,fntx,fnty,fntz,fnx,fny,fnz,ftx,fty,ftz]=textread( ...
            [Foldersdir_pairwise '\' num2str(stepnumber(sample_step(i,1),1)) '.txt'],'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','headerlines',9);
        % % % [contacts,xi,yi,~,xj,yj,~,id1,id2,id3,fntx,fnty,fntz,fnx,fny,fnz,ftx,fty,ftz]
        pairfn=sqrt(fnx.^2+fny.^2+fnz.^2);
        pairdist=sqrt((x1-x2).^2+(y1-y2).^2+(z1-z2).^2);
        pl_tmp1_prepare=[id1,id2];
        pl_tmp1_prepare(:,3)=pairdist;% pair distance
        pl_tmp1_prepare(:,4)=id_centroids_radius_vxyz(pl_tmp1_prepare(:,1),5)+id_centroids_radius_vxyz(pl_tmp1_prepare(:,2),5);% r1+r2
        pl_tmp1_prepare(:,5)=pairfn;
%         pl_tmp1_prepare_vxindex=id_centroids_radius_vxyz(pl_tmp1_prepare(:,1),6);
%         pl_tmp1_prepare_vzindex=id_centroids_radius_vxyz(pl_tmp1_prepare(:,1),8);
        pl_tmp1_prepare=pl_tmp1_prepare(ismember(pl_tmp1_prepare(:,1),id_est_filtered)|ismember(pl_tmp1_prepare(:,2),id_est_filtered),:);
%         pl_tmp1=pl_tmp1_prepare;
        pl_tmp1_abovefnmean=pl_tmp1_prepare(pl_tmp1_prepare(:,5)>=mean(pl_tmp1_prepare(:,5)),:);%f_pair>f_mean
        pl_tmp1=pl_tmp1_abovefnmean;
        pl_tmp1_record{i,1}=pl_tmp1;
        fmean{i,1}=mean(pl_tmp1_prepare(:,5));
%         figure(2)
%         unique_id_pl_tmp1=unique([pl_tmp1_abovefnmean(:,1);pl_tmp1_abovefnmean(:,2)]);
%         centroids_est_tmp=centroids_est(unique_id_pl_tmp1,:);
%         scatter3(centroids_est_tmp(:,1),centroids_est_tmp(:,2),centroids_est_tmp(:,3))
%         axis equal
    end
    toc

end

%% save to mat file

   atominfo_datafilename=['E:\RoDrtest\RoDr_alphastart0-1\model' num2str(modelnum(modeli,1)) '\model' num2str(modelnum(modeli,1)) '_fmean.mat'];

   save(atominfo_datafilename,'pl_tmp1_record','-v7.3');
       clear pl_tmp1_record
   save(atominfo_datafilename,'fmean','-append');
       clear fmean

end
end
% figure(10)
% cluster_id_size(cluster_id_size(:,2)>50,:)=[];
% histogram(cluster_id_size(:,2),"Normalization","Probability");
% set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')

%% tmp check
% % your data
% a = [1,2,3,4,5,6]';
% b = [10,11,12,13,14,15]';
% Interval = [a b];
% % number to check
% c = 8;
% % this finds the index of he rows(2) that have x in between 
% idx = find(c > Interval(:,1) & c < Interval(:,2));
% % number of intervals with positive check
% numIdx = sum(c > Interval(:,1) & c < Interval(:,2));
%% test graph error
% figure(100)
% s = [1 1 1 1 1 2 2 7 7 9 3 3 1 4 4 5 8];
% t = [2 3 4 5 7 6 7 5 9 6 6 10 10 10 8 8 9];
% % weights = [1 1 1 1 3 3 2 4 1 6 2 8 8 9 3 2 10 12 15 16];
% G = graph(s,t,10);
% x = [0 0.5 -0.5 -0.5 0.5 0 1.5 0 2 -1.5];
% y = [0 0.5 0.5 -0.5 -0.5 2 0 -2 0 0];
% z = [5 3 3 3 3 0 1 0 0 1];
% plot(G,'XData',x,'YData',y,'ZData',z);

%% test atan2
% a=[1 0 0];
% b=[-1,1,0];
% c=atan2(norm(cross(a,b)),dot(a,b));

%% check biconnect cut
% % s = [1 1 2 2 3 4 4 5 6 6 7 7 8];
% % t = [2 3 3 4 4 5 7 6 7 10 8 9 9];
% s = [1 1 2 2 3 4 5 6 6 7 7 8 11 11];
% t = [2 3 3 4 4 11 6 7 10 8 9 9 5 7];
% G = graph(s,t);
% p = plot(G,'LineWidth',2);
% p.EdgeCData =  biconncomp(G);
% 
% bincell = biconncomp(G, 'OutputForm', 'cell');
% n = length(bincell);
% 
% figure(2)
% for ii = 1:n
%     subplot(3,2,ii)
%     plot(subgraph(G, bincell{ii}), 'NodeLabel', bincell{ii});
% end
%% testing parallel and perpendicular vectors
function [n1parak,n2parak,n1perpk,n2perpk]=n12_para_perpen_k(n1,n2,k)
%% uncomment this part for testing
% B = [1 1 0]; % you make this happen. B should be [1 x 3] for the following to work
% v = [-1 -1 -1]; % you make this happen. Same for v
% v_par_norm = dot(v,B)/norm(B);
% v_par = v_par_norm*B/norm(B);
% v_perp = v - v_par_norm*B/norm(B);
%% this part is for function
%----------n1 and k-----------
n1_projecton_k = dot(n1,k)/norm(k); %scaler
n1parak = n1_projecton_k*k/norm(k); 
n1perpk = n1 - n1parak;
%----------n2 and k-----------
n2_projecton_k = dot(n2,k)/norm(k);
n2parak = n2_projecton_k*k/norm(k);
n2perpk = n2 - n2parak;
end

function createfigure(X1, Y1, Z1, Size1, Color1, X2, Y2, Z2, Size2, Color2)
%CREATEFIGURE(X1, Y1, Z1, Size1, Color1, X2, Y2, Z2, Size2, Color2)
%  X1:  vector of scatter3 x data
%  Y1:  vector of scatter3 y data
%  Z1:  vector of scatter3 z data
%  SIZE1:  vector of scatter3 size data
%  COLOR1:  vector of scatter3 color data
%  X2:  vector of scatter3 x data
%  Y2:  vector of scatter3 y data
%  Z2:  vector of scatter3 z data
%  SIZE2:  vector of scatter3 size data
%  COLOR2:  vector of scatter3 color data

%  Auto-generated by MATLAB on 23-Jan-2023 22:14:58

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
axis off
hold(axes1,'on');

% Create scatter3
scatter3(X1,Y1,Z1,Size1,Color1,'MarkerEdgeColor','none',...
    'MarkerFaceColor','flat');

% Create scatter3
scatter3(X2,Y2,Z2,Size2,Color2,'MarkerEdgeColor','none',...
    'MarkerFaceColor','flat');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[-0.08 0.1]);
view(axes1,[0.754294656644959 -1.4026523406981]);
box(axes1,'on');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'BoxStyle','full','CLim',[1 10000],'ClippingStyle','rectangle',...
    'ColorScale','log','Colormap',...
    [0.18995 0.07176 0.23217;0.19483 0.08339 0.26149;0.19956 0.09498 0.29024;0.20415 0.10652 0.31844;0.2086 0.11802 0.34607;0.21291 0.12947 0.37314;0.21708 0.14087 0.39964;0.22111 0.15223 0.42558;0.225 0.16354 0.45096;0.22875 0.17481 0.47578;0.23236 0.18603 0.50004;0.23582 0.1972 0.52373;0.23915 0.20833 0.54686;0.24234 0.21941 0.56942;0.24539 0.23044 0.59142;0.2483 0.24143 0.61286;0.25107 0.25237 0.63374;0.25369 0.26327 0.65406;0.25618 0.27412 0.67381;0.25853 0.28492 0.693;0.26074 0.29568 0.71162;0.2628 0.30639 0.72968;0.26473 0.31706 0.74718;0.26652 0.32768 0.76412;0.26816 0.33825 0.7805;0.26967 0.34878 0.79631;0.27103 0.35926 0.81156;0.27226 0.3697 0.82624;0.27334 0.38008 0.84037;0.27429 0.39043 0.85393;0.27509 0.40072 0.86692;0.27576 0.41097 0.87936;0.27628 0.42118 0.89123;0.27667 0.43134 0.90254;0.27691 0.44145 0.91328;0.27701 0.45152 0.92347;0.27698 0.46153 0.93309;0.2768 0.47151 0.94214;0.27648 0.48144 0.95064;0.27603 0.49132 0.95857;0.27543 0.50115 0.96594;0.27469 0.51094 0.97275;0.27381 0.52069 0.97899;0.27273 0.5304 0.98461;0.27106 0.54015 0.9893;0.26878 0.54995 0.99303;0.26592 0.55979 0.99583;0.26252 0.56967 0.99773;0.25862 0.57958 0.99876;0.25425 0.5895 0.99896;0.24946 0.59943 0.99835;0.24427 0.60937 0.99697;0.23874 0.61931 0.99485;0.23288 0.62923 0.99202;0.22676 0.63913 0.98851;0.22039 0.64901 0.98436;0.21382 0.65886 0.97959;0.20708 0.66866 0.97423;0.20021 0.67842 0.96833;0.19326 0.68812 0.9619;0.18625 0.69775 0.95498;0.17923 0.70732 0.94761;0.17223 0.7168 0.93981;0.16529 0.7262 0.93161;0.15844 0.73551 0.92305;0.15173 0.74472 0.91416;0.14519 0.75381 0.90496;0.13886 0.76279 0.8955;0.13278 0.77165 0.8858;0.12698 0.78037 0.8759;0.12151 0.78896 0.86581;0.11639 0.7974 0.85559;0.11167 0.80569 0.84525;0.10738 0.81381 0.83484;0.10357 0.82177 0.82437;0.10026 0.82955 0.81389;0.0975 0.83714 0.80342;0.09532 0.84455 0.79299;0.09377 0.85175 0.78264;0.09287 0.85875 0.7724;0.09267 0.86554 0.7623;0.0932 0.87211 0.75237;0.09451 0.87844 0.74265;0.09662 0.88454 0.73316;0.09958 0.8904 0.72393;0.10342 0.896 0.715;0.10815 0.90142 0.70599;0.11374 0.90673 0.69651;0.12014 0.91193 0.6866;0.12733 0.91701 0.67627;0.13526 0.92197 0.66556;0.14391 0.9268 0.65448;0.15323 0.93151 0.64308;0.16319 0.93609 0.63137;0.17377 0.94053 0.61938;0.18491 0.94484 0.60713;0.19659 0.94901 0.59466;0.20877 0.95304 0.58199;0.22142 0.95692 0.56914;0.23449 0.96065 0.55614;0.24797 0.96423 0.54303;0.2618 0.96765 0.52981;0.27597 0.97092 0.51653;0.29042 0.97403 0.50321;0.30513 0.97697 0.48987;0.32006 0.97974 0.47654;0.33517 0.98234 0.46325;0.35043 0.98477 0.45002;0.36581 0.98702 0.43688;0.38127 0.98909 0.42386;0.39678 0.99098 0.41098;0.41229 0.99268 0.39826;0.42778 0.99419 0.38575;0.44321 0.99551 0.37345;0.45854 0.99663 0.3614;0.47375 0.99755 0.34963;0.48879 0.99828 0.33816;0.50362 0.99879 0.32701;0.51822 0.9991 0.31622;0.53255 0.99919 0.30581;0.54658 0.99907 0.29581;0.56026 0.99873 0.28623;0.57357 0.99817 0.27712;0.58646 0.99739 0.26849;0.59891 0.99638 0.26038;0.61088 0.99514 0.2528;0.62233 0.99366 0.24579;0.63323 0.99195 0.23937;0.64362 0.98999 0.23356;0.65394 0.98775 0.22835;0.66428 0.98524 0.2237;0.67462 0.98246 0.2196;0.68494 0.97941 0.21602;0.69525 0.9761 0.21294;0.70553 0.97255 0.21032;0.71577 0.96875 0.20815;0.72596 0.9647 0.2064;0.7361 0.96043 0.20504;0.74617 0.95593 0.20406;0.75617 0.95121 0.20343;0.76608 0.94627 0.20311;0.77591 0.94113 0.2031;0.78563 0.93579 0.20336;0.79524 0.93025 0.20386;0.80473 0.92452 0.20459;0.8141 0.91861 0.20552;0.82333 0.91253 0.20663;0.83241 0.90627 0.20788;0.84133 0.89986 0.20926;0.8501 0.89328 0.21074;0.85868 0.88655 0.2123;0.86709 0.87968 0.21391;0.8753 0.87267 0.21555;0.88331 0.86553 0.21719;0.89112 0.85826 0.2188;0.8987 0.85087 0.22038;0.90605 0.84337 0.22188;0.91317 0.83576 0.22328;0.92004 0.82806 0.22456;0.92666 0.82025 0.2257;0.93301 0.81236 0.22667;0.93909 0.80439 0.22744;0.94489 0.79634 0.228;0.95039 0.78823 0.22831;0.9556 0.78005 0.22836;0.96049 0.77181 0.22811;0.96507 0.76352 0.22754;0.96931 0.75519 0.22663;0.97323 0.74682 0.22536;0.97679 0.73842 0.22369;0.98 0.73 0.22161;0.98289 0.7214 0.21918;0.98549 0.7125 0.2165;0.98781 0.7033 0.21358;0.98986 0.69382 0.21043;0.99163 0.68408 0.20706;0.99314 0.67408 0.20348;0.99438 0.66386 0.19971;0.99535 0.65341 0.19577;0.99607 0.64277 0.19165;0.99654 0.63193 0.18738;0.99675 0.62093 0.18297;0.99672 0.60977 0.17842;0.99644 0.59846 0.17376;0.99593 0.58703 0.16899;0.99517 0.57549 0.16412;0.99419 0.56386 0.15918;0.99297 0.55214 0.15417;0.99153 0.54036 0.1491;0.98987 0.52854 0.14398;0.98799 0.51667 0.13883;0.9859 0.50479 0.13367;0.9836 0.49291 0.12849;0.98108 0.48104 0.12332;0.97837 0.4692 0.11817;0.97545 0.4574 0.11305;0.97234 0.44565 0.10797;0.96904 0.43399 0.10294;0.96555 0.42241 0.09798;0.96187 0.41093 0.0931;0.95801 0.39958 0.08831;0.95398 0.38836 0.08362;0.94977 0.37729 0.07905;0.94538 0.36638 0.07461;0.94084 0.35566 0.07031;0.93612 0.34513 0.06616;0.93125 0.33482 0.06218;0.92623 0.32473 0.05837;0.92105 0.31489 0.05475;0.91572 0.3053 0.05134;0.91024 0.29599 0.04814;0.90463 0.28696 0.04516;0.89888 0.27824 0.04243;0.89298 0.26981 0.03993;0.88691 0.26152 0.03753;0.88066 0.25334 0.03521;0.87422 0.24526 0.03297;0.8676 0.2373 0.03082;0.86079 0.22945 0.02875;0.8538 0.2217 0.02677;0.84662 0.21407 0.02487;0.83926 0.20654 0.02305;0.83172 0.19912 0.02131;0.82399 0.19182 0.01966;0.81608 0.18462 0.01809;0.80799 0.17753 0.0166;0.79971 0.17055 0.0152;0.79125 0.16368 0.01387;0.7826 0.15693 0.01264;0.77377 0.15028 0.01148;0.76476 0.14374 0.01041;0.75556 0.13731 0.00942;0.74617 0.13098 0.00851;0.73661 0.12477 0.00769;0.72686 0.11867 0.00695;0.71692 0.11268 0.00629;0.7068 0.1068 0.00571;0.6965 0.10102 0.00522;0.68602 0.09536 0.00481;0.67535 0.0898 0.00449;0.66449 0.08436 0.00424;0.65345 0.07902 0.00408;0.64223 0.0738 0.00401;0.63082 0.06868 0.00401;0.61923 0.06367 0.0041;0.60746 0.05878 0.00427;0.5955 0.05399 0.00453;0.58336 0.04931 0.00486;0.57103 0.04474 0.00529;0.55852 0.04028 0.00579;0.54583 0.03593 0.00638;0.53295 0.03169 0.00705;0.51989 0.02756 0.0078;0.50664 0.02354 0.00863;0.49321 0.01963 0.00955;0.4796 0.01583 0.01055],...
    'DataAspectRatio',[1 1 1],'XLimitMethod','tight','YLimitMethod','tight',...
    'ZLimitMethod','tight');
end
