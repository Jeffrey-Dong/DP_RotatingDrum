clc;
clear;

addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\dataanalysis\');
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\2d_packing\data_analysis_code\functions\');
[num,txt,raw]=xlsread(['C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\parametric_study.xlsx'],1);
keypara=num(4:end,:);
density_threshold=0.6;

pldist_threshold_coe_dry=1.05;%this is only for creating pairlist using pair distance method when dry condition is found with 0 liquid content

clustermethod="angle";
angle_threshold=pi/36;
% clustermethod="dotvel";
dotvel_coe=3;

% modelnum=[31,33,35,37,39,41]';
% modelnum=[49:52]';
% modelnum=[99:102]';
% modelnum=[97,99:102]';
% modelnum=[86:94]';
% modelnum=[81:84,77]';
% modelnum=[103:111]';
% modelnum=[113:114]';
% modelnum=[119,121:123,124,126:128,134]';
modelnum=[119,134,121:123]';
% modelnum=[134]';
%%
for modeli=1:size(modelnum,1)
%% clear to save memory when jump to the next case
clearvars -except modeli modelnum num txt raw keypara density angle_threshold pldist_threshold_coe_dry clustermethod dotvel_coe;

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
if liquid_content==0
    %do nothing
else
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
end
% toc
%% Prepare container for atom and contact pair info
% sample_step=[300:20:Rtxt_atom]';
% if Rtxt>=650
%     sample_step=[400:10:650]';
% else
    sample_step=[450:10:Rtxt]';
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
        [id,~,x,y,z,~,~,~,vx,vy,vz,~,~,~,~,~,~,radius_est,VliqonVatom]=textread( ...
            txt_read_atominfo(sample_step(i)).name,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','headerlines',9);
        vel=sqrt(vx.^2+vy.^2+vz.^2);
    end
    id_centroids_radius_vxyz=[id,x,y,z,radius_est,vx,vy,vz];
    id_centroids_radius_vxyz=sortrows(id_centroids_radius_vxyz,1);
    centroids_est=id_centroids_radius_vxyz(:,2:4);%each step
    id_est=id_centroids_radius_vxyz(:,1);%each step
    vxyz_est=id_centroids_radius_vxyz(:,6:8);
    id_centroids_radius_vxyz_record{i,1}=id_centroids_radius_vxyz;
%     toc
    %% pair information: create pairlist based on pair distance
    tic
    if liquid_content==0
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
        centroids_est_vxzless0=centroids_est(vxyz_est(:,1)<0 & vxyz_est(:,3)<0,:);
        id_est_vxzless0=id_est(vxyz_est(:,1)<0 & vxyz_est(:,3)<0,:);
        for pli=1:size(centroids_est_vxzless0,1)
            xbound=[centroids_est_vxzless0(pli,1)-1.5*d_p centroids_est_vxzless0(pli,1)+1.5*d_p];
            ybound=[centroids_est_vxzless0(pli,2)-1.5*d_p centroids_est_vxzless0(pli,2)+1.5*d_p];
            zbound=[centroids_est_vxzless0(pli,3)-1.5*d_p centroids_est_vxzless0(pli,3)+1.5*d_p];
            id_inthisloop=id_est_vxzless0(pli+1:end,1);
            index_x1=id_inthisloop(centroids_est_vxzless0(pli+1:end,1)>xbound(1));
            index_x2=id_inthisloop(centroids_est_vxzless0(pli+1:end,1)<xbound(2));
            index_y1=id_inthisloop(centroids_est_vxzless0(pli+1:end,2)>ybound(1));
            index_y2=id_inthisloop(centroids_est_vxzless0(pli+1:end,2)<ybound(2));
            index_z1=id_inthisloop(centroids_est_vxzless0(pli+1:end,3)>zbound(1));
            index_z2=id_inthisloop(centroids_est_vxzless0(pli+1:end,3)<zbound(2));
            index_x=index_x1(ismember(index_x1,index_x2));
            index_y=index_y1(ismember(index_y1,index_y2));
            index_z=index_z1(ismember(index_z1,index_z2));
            index_xy=index_x(ismember(index_x,index_y));
            index_xyz=index_xy(ismember(index_xy,index_z));
            pl=[pl;[ones(size(index_xyz,1),1)*id_est_vxzless0(pli,1), index_xyz]];
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
        pl_tmp1=[id1,id2];
        pl_tmp1(:,3)=pairdist;% pair distance
        pl_tmp1(:,4)=id_centroids_radius_vxyz(pl_tmp1(:,1),5)+id_centroids_radius_vxyz(pl_tmp1(:,2),5);% r1+r2
%         vx=id_centroids_radius_vxyz(pl_tmp1(:,1),6);
%         vz=id_centroids_radius_vxyz(pl_tmp1(:,1),8);
%         pl_tmp1(vx>0|vz>0,:)=[];
        % pairfnt=sqrt(fntx.^2+fnty.^2+fntz.^2);
    % pairft=sqrt(ftx.^2+fty.^2+ftz.^2);
    end
    toc
%% graph all pair network based on contact pair force
% id1_2_x1_y1_z1=[id1,id2,x1,y1,z1];
% id_centre=[id,x,y,z];
% id_centre=sortrows(id_centre,1);
% id_centre(find(id_centre(:,1)>max([id1;id2]),1):end,:)=[];
% 
% id1_2_tmp=[id1,id2];
% switch_index=(id1_2_tmp(:,1)>id1_2_tmp(:,2));
% id1_2=id1_2_tmp;
% id1_2(switch_index,1)=id1_2_tmp(switch_index,2);
% id1_2(switch_index,2)=id1_2_tmp(switch_index,1);
% id1_2=sortrows(id1_2,1);
% 
% Gallpair=graph(id1_2(:,1),id1_2(:,2));
% tmp_x1=id_centre(id1_2(:,1),2);
% tmp_y1=id_centre(id1_2(:,1),3);
% tmp_z1=id_centre(id1_2(:,1),4);
% tmp_x2=id_centre(id1_2(:,2),2);
% tmp_y2=id_centre(id1_2(:,2),3);
% tmp_z2=id_centre(id1_2(:,2),4);
% tmpdist=sqrt((tmp_x1-tmp_x2).^2+(tmp_y1-tmp_y2).^2+(tmp_z1-tmp_z2).^2);
% plot(Gallpair,'XData',id_centre(:,2),'YData',id_centre(:,3),'ZData',id_centre(:,4));

%% graph all pair network based on distance based neighbour
% figure(1)
% Gallpair=graph(id1,id2);
% plot(Gallpair,'XData',centroids_est(:,1),'YData',centroids_est(:,2),'ZData',centroids_est(:,3));
% axis equal;
% ylim([-0.1 0.1])
    %% correct id etc.
%     [id,centroids,radius_corrected,pairid_force_dis_sumradij,Fij,Wij,Bij,ave_force] = ...
%         correct_id_and_obtain_matrix_3D_RotationDrum(id,x,y,z,radius,id1,id2,pairfn,VliqonVatom,vx,vy,vz,pairdist_directcal);
    % % centr a99*99
    % 966666666666oids=zeros(size(x,1),3);
    % centroids=[];
    % centroids(:,1)=x;
    % centroids(:,2)=y;
    % centroids(:,3)=z;
%     for j=1:size(x1,1)
%     normalise_pairfn=pairfn;
%      a=[x1(j,1); x2(j,1)];
%      b=[y1(j,1); y2(j,1)];
%      c=[z1(j,1); z2(j,1)];
%      figure(1);hold on;
%     plot3(a,b,c);
%     end
%% Coordination number
% crit_l=1.0001;
%     [Zc,Cn_ep,Cn_ep_neibid]=Z_positionbase_rotationdrum(radius_est,[x,y,z],crit_l,idest);
% % figure(1)
% % scatter3(x,y,z,[],Cn_ep,"filled");
% % title('Coordination number')
% % figure(2)
% % scatter3(x,y,z,[],vel,"filled");
% % title('velocity')
% % figure(3)
% % plot(Cn_ep,vel,'k.');
% % xlabel('Coordination number');
% % ylabel('Particle velocity');
% % xlim([-1 9]);
% % figure(4)
% % histogram(Cn_ep,'Normalization','probability');
% % xlabel('Coordination number');
% % ylabel('Probability');
% % xlim([-1 9]);

%% import cell format velocity

%     load([Foldersdir_atominfo '\model58_velmatrix1529frames']);
% load([Foldersdir_atominfo '\model46_velmatrix400frames']);
load(['E:\RoDrtest\RoDr_alphastart0-1\model' num2str(modelnum(modeli,1)) ...
    '\model' num2str(modelnum(modeli,1)) '_velmatrix_frames']);
%     figure(1)
%     contourf(Xpf,Zpf,velxincell_step1_total_ave)
for vmi=1:size(id_est,1)
    xtmp=centroids_est(id_est==vmi,1);
    ztmp=centroids_est(id_est==vmi,3);
    vm_x_index=find(Xpf(:,1)>xtmp,1)-1;%locate the index of bigger one for the interval [smaller,bigger]
    vm_z_index=find(Zpf(1,:)>ztmp,1)-1;
    vm_ep(vmi,1:2)=[meanvelxincell(vm_x_index,vm_z_index) meanvelzincell(vm_x_index,vm_z_index)];%mean velocity for each particle
end

%% import DAOR and obtain n_flow
% load([Foldersdir_atominfo '\model55_DAOR_gof']);
% DAOR_steady=mean(DAOR(revolution(:,1)>100,1));
% v_flow_angle=DAOR_steady-pi;
% n_flow=[cos(v_flow_angle) 0 sin(v_flow_angle)]; %the unit vector of average velocity in the flow region, which corresponds to to DAOR

%% Obtain velocity fluctuation for each particle
%--------------veloctiy fluctuation and its unit vector---------------
n_ep=zeros(size(id_est,1),1);
for idi=1:size(id_est,1)
    v_fluctuation=[vxyz_est(idi,1)-vm_ep(idi,1) vy(idi) vxyz_est(idi,3)-vm_ep(idi,2)];
    n_ep(idi,1:3)=v_fluctuation;
%     /norm(v_fluctuation);
end
%% --------------screen all pairs for clustering using angles between each pair of velocity fluctuation unit vectors---------
if clustermethod=="angle"
    %%
    tic
    pair_vel_angle_pairmark=zeros(size(pl_tmp1,1),1);
    pair_vel_flowzone_pairmark=zeros(size(pl_tmp1,1),1);
    for pairi=1:size(pl_tmp1,1)
        index_id1=find(id_est==id1(pairi,1));
        index_id2=find(id_est==id2(pairi,1));
        n1=[n_ep(index_id1,1) n_ep(index_id1,2) n_ep(index_id1,3)];
        n2=[n_ep(index_id2,1) n_ep(index_id2,2) n_ep(index_id2,3)];
        pair_vel_angle_pairmark(pairi,1) = atan2(norm(cross(n1,n2)),dot(n1,n2));
        %% extract flow zone based on velocity
%         if vx(index_id1)<0 && vz(index_id1)<0 && vx(index_id2)<0 && vz(index_id2)<0
%             pair_vel_flowzone_pairmark(pairi,1)=1;
%         end
        %% extract pair correlation based on velocity angle
        if pair_vel_angle_pairmark(pairi,1)<=angle_threshold
            pair_vel_angle_pairmark(pairi,2)=1;
        else
            pair_vel_angle_pairmark(pairi,2)=0;
        end
    end
    pl_tmp2=pl_tmp1;
    pl_tmp2(pair_vel_angle_pairmark(:,2)==0,:)=[];
    % pl_tmp2(pl_tmp2(:,3)>pl_tmp2(:,4)*1.032,:)=[];
    tmpidxyz=id_centroids_radius_vxyz;
    % tmpidxyz(non_corre_id,:)=[];
    toc
elseif clustermethod=="dotvel"
%% --------screen all pairs for clustering using velocity dot product------------
    % correlation each pair
    pair_vel_corre=zeros(size(pl_tmp1,1),2);
    for pairi=1:size(pl_tmp1,1)
        index_id1=find(id==id1(pairi,1));
        index_id2=find(id==id2(pairi,1));
        if vx(index_id1)<0 && vz(index_id1)<0 && vx(index_id2)<0 && vz(index_id2)<0
        n1=[n_ep(index_id1,1) n_ep(index_id1,2) n_ep(index_id1,3)];
        n2=[n_ep(index_id2,1) n_ep(index_id2,2) n_ep(index_id2,3)];
        posi1=[x(index_id1,1) y(index_id1,1) z(index_id1,1)];
        posi2=[x(index_id2,1) y(index_id2,1) z(index_id2,1)];
        k=(posi1-posi2)/norm(posi1-posi2);
    %     An1k=atan2(norm(cross(n1,k)),dot(n1,k));%angle between n1 and k
    %     An2k=atan2(norm(cross(n2,k)),dot(n1,k));%angle between n2 and k
        [n1parak,n2parak,n1perpk,n2perpk]=n12_para_perpen_k(n1,n2,k);
        pair_vel_corre(pairi,1) = dot(n1parak,n2parak);
    %     .*dot(n1perpk,n2perpk);
        pair_vel_corre(pairi,2) = dot(n1perpk,n2perpk);
    %     if pair_vel_angle_pairmark(pairi,1)<=pi/9 && vx(id==id1(pairi,1))<0 && vz(id==id1(pairi,1))<0
    %         pair_vel_angle_pairmark(pairi,2)=1;
    %     else
    %         continue
    %     end
        else
            pair_vel_corre(pairi,1)=0;
            pair_vel_corre(pairi,2)=0;
        end
    end
    %%
    velo_fluc_coore_ep=zeros(size(id,1),2);
    for idi=1:size(id,1)
        indextmp=find(pl_tmp1(:,1)==idi | pl_tmp1(:,2)==idi);
        N_this_cluster=size(unique([pl_tmp1(indextmp,1),pl_tmp1(indextmp,2)]),1);
        velo_fluc_coore_ep(idi,1)=sum(pair_vel_corre(indextmp,1))/N_this_cluster;%parallel
        velo_fluc_coore_ep(idi,2)=sum(pair_vel_corre(indextmp,2))/N_this_cluster;%perpendicular
    end
    % % ====================extract high correlated particles from pl===========
    % the id of particles that have correlation less than "correlation_threshold" and "nan"
    correlation_threshold=(omega/180*pi*0.1)^2*dotvel_coe;
    non_corre_id=[id_centroids_radius_vxyz(isnan(velo_fluc_coore_ep(:,1)),1);
        id_centroids_radius_vxyz(velo_fluc_coore_ep(:,1)<correlation_threshold,1);
        id_centroids_radius_vxyz(isnan(velo_fluc_coore_ep(:,2)),1);
        id_centroids_radius_vxyz(velo_fluc_coore_ep(:,2)<correlation_threshold,1)];
    non_corre_id=unique(non_corre_id);
    % the index of particles that correlation is smaller than "correlation_threshold" in the pair list
    index_corre_pair=[];
    for nclui=1:size(non_corre_id,1)
        index_corre_pair=[index_corre_pair;find(pl_tmp1(:,1)==non_corre_id(nclui,1));find(pl_tmp1(:,2)==non_corre_id(nclui,1))];
    end
    index_corre_pair=unique(index_corre_pair);
    % --------------------remove the particles that correlation is smaller than "correlation_threshold" in the pair list---------
    pl_tmp2=pl_tmp1;
    pl_tmp2(index_corre_pair,:)=[];
    pl_tmp2(pl_tmp2(:,3)>pl_tmp2(:,4)*1.032,:)=[];
    tmpidxyz=id_centroids_radius_vxyz;
    tmpidxyz(non_corre_id,:)=[];
end
%% create graph using residuel pairlist
Gprecluster=graph(pl_tmp2(:,1),pl_tmp2(:,2));
% % max id number in residule connected pairs
% figure(4);
% idmax_residue=max([pl_tmp2(:,1);pl_tmp2(:,2)]);
% plot(Gprecluster,'XData',id_centroids_radius_vxyz(1:idmax_residue,2), ...
%                  'YData',id_centroids_radius_vxyz(1:idmax_residue,3), ...
%                  'ZData',id_centroids_radius_vxyz(1:idmax_residue,4));
% axis equal;
% ylim([-0.1 0.1]);

%% extract high correlated particles from centroids and id
%    centroids_tmp=centroids;
%    centroids_tmp(isnan(velo_fluc_coore_ep(:,1)) | velo_fluc_coore_ep(:,1)<0.1,:)=[];
%    velo_fluc_coore_ep(isnan(velo_fluc_coore_ep(:,1)) | velo_fluc_coore_ep(:,1)<0.1,:)=[];
% %    velo_fluc_coore_ep(isnan(velo_fluc_coore_ep(:,2)),:)=[];
%     
% %     figure(3);scatter3(centroids_tmp(:,1),centroids_tmp(:,2),centroids_tmp(:,3),[],velo_fluc_coore_ep(:,1));axis equal;
% %     figure(4);scatter3(centroids_tmp(:,1),centroids_tmp(:,2),centroids_tmp(:,3),[],velo_fluc_coore_ep(:,2));axis equal;
%     figure(5);scatter3(centroids_tmp(:,1),centroids_tmp(:,2),centroids_tmp(:,3),[],velo_fluc_coore_ep(:,1));axis equal;
% %     figure(6);scatter3(centroids_tmp(:,1),centroids_tmp(:,2),centroids_tmp(:,3),[],velo_fluc_coore_ep(:,1)-velo_fluc_coore_ep(:,2));axis equal;

    %% velocity correlation method: correlation between all pair contact pairs
%     pair_vel_angle_clustermark=zeros(size(fnx,1),1);
%     for pairi=1:size(fnx,1)
%         index_id1=find(id==id1(pairi));
%         index_id2=find(id==id2(pairi));
%         v1=[vx(index_id1) vy(index_id1) vz(index_id1)];
%         v2=[vx(index_id2) vy(index_id2) vz(index_id2)];
%         pair_vel_angle_clustermark(pairi,1) = atan2(norm(cross(v1,v2)),dot(v1,v2));
%         if pair_vel_angle_clustermark(pairi,1)<=pi/18
%             pair_vel_angle_clustermark(pairi,2)=1;
%         else
%             continue
%         end
%     end


 %% cluster
% c=clusterdata([x,y,z],20);
% figure(5)
% scatter3(x,y,z,[],c,"filled");
% figure(5)
% scatter3(x1,y1,z1,[],pair_vel_angle_clustermark(:,2),"filled");
%% function

%% for temperary test on angle between two vectors
% P1=[1,1,1];P2=[1,1,1];
% a = atan2(norm(cross(P1,P2)),dot(P1,P2));
% a/pi*180

%% cluster identification based on graph theory and the angle between pair velocity vector
% pairid=[id1,id2];
% cluster_pairid=pairid(pair_vel_angle_pairmark(:,2)==1,:);
% G=graph(cluster_pairid(:,1),cluster_pairid(:,2));
% figure(6)
% plot(G)

%% Apply conncomp function
[bin,binsize] = conncomp(Gprecluster,'Type','weak');
bincell=conncomp(Gprecluster,'OutputForm','cell');
n=length(bincell');
% % idx = binsize(bin) == max(binsize);%obtain max clusters
% idx = binsize(bin) > 0;%obtain max clusters
% SG = subgraph(Gprecluster, idx);
% figure(7)
% plot(SG);

%% Apply biconncomp function
% bincell = biconncomp(Gprecluster, 'OutputForm', 'cell');
% n = length(bincell);
% cluster_binrecord{i,1}=bincell;
% 
% % % plot each cluster
% % figure(7)
% % for ii = 1:n
% %     subplot(20,20,ii)
% %     plot(subgraph(SG, bincell{ii}), 'NodeLabel', bincell{ii});
% % end

%% release cluster id from bincell for ploting
cluster_bin={};
count=0;
for clusi=1:n
    if size(bincell{1,clusi}',1)>1 || size(bincell{1,clusi}',1)==1 && id_centroids_radius_vxyz(bincell{1,clusi}',6)<0 && id_centroids_radius_vxyz(bincell{1,clusi}',8)<0
        count=count+1;
        cluster_id_size(count,1:2)=[clusi size(bincell{1,clusi}',1)];
        cluster_bin{count,1}=bincell{1,clusi};
        
%     elseif 
%         count=count+1;
%         cluster_id_size(count,1:2)=[clusi size(bincell{1,clusi}',1)];
%-------------for plot---------------
    cluster_id=ones(size(bincell{1,clusi}',1),1)*clusi;
    atomid_clusterid=[atomid_clusterid;bincell{1,clusi}',cluster_id];
    atomid_clustersize=[atomid_clustersize;bincell{1,clusi}',ones(size(bincell{1,clusi}',1),1)*size(bincell{1,clusi}',1)];
    else
        %continue
    end
    
end

cluster_binrecord{i,1}=cluster_bin;
cluster_sizerecord=[cluster_sizerecord;cluster_id_size];
% %------------for plot-------------
for clusj=1:size(atomid_clusterid,1)
cluster_x(clusj,1)=centroids_est(id_est==atomid_clusterid(clusj,1),1);
cluster_y(clusj,1)=centroids_est(id_est==atomid_clusterid(clusj,1),2);
cluster_z(clusj,1)=centroids_est(id_est==atomid_clusterid(clusj,1),3);
end

% figure(1)
% histogram(cluster_sizerecord(:,2),"Normalization","Probability")
% set(gca, 'XScale', 'log')
% set(gca, 'YScale', 'log')
% xlabel('cluster size')
% ylabel('probability')
%% 
figure(8);
fig8=scatter3(x(:,1),y(:,1),z(:,1),0.01,'b.');
hold on;
fig8=scatter3(cluster_x,cluster_y,cluster_z,[],atomid_clusterid(:,2),"filled");
axis equal
title('cluster id')
colorbar;
% % saveas(fig8, 'cluster_id_example.png'); % Save as PNG to avoid jpeg artifacts.
% % close figure 8;
%% 
figure(9)
fig9=scatter3(x(:,1),y(:,1),z(:,1),0.1,'b.');
hold on;
fig9=scatter3(cluster_x,cluster_y,cluster_z,[],atomid_clustersize(:,2),"filled");
axis equal
title('cluster size')
colorbar; 
% saveas(fig9, 'cluster_size_example.png'); % Save as PNG to avoid jpeg artifacts.
% close figure 9;
%% plot single clusters
% figure(10)
% max_cluster_x=[];
% max_cluster_y=[];
% max_cluster_z=[];
% tmp_atomid_clustersize=atomid_clustersize;
% tmp_atomid_clustersize(tmp_atomid_clustersize(:,2)>300,:)=[];
% % the id of particles in the largest cluster
% [~,id_maxcluster]=ismember(tmp_atomid_clustersize(tmp_atomid_clustersize(:,2)==max(tmp_atomid_clustersize(:,2)),1),id);
% max_cluster_x(:,1)=x(id_maxcluster);
% max_cluster_y(:,1)=y(id_maxcluster);
% max_cluster_z(:,1)=z(id_maxcluster);
% scatter3(max_cluster_x,max_cluster_y,max_cluster_z,"MarkerFaceColor",[0 1 0]);

end

%% save to mat file
   if clustermethod=="dotvel"
   atominfo_datafilename=['E:\RoDrtest\RoDr_alphastart0-1\model' num2str(modelnum(modeli,1)) '\model' num2str(modelnum(modeli,1)) '_clustermethod_' char(clustermethod) '_coe' num2str(dotvel_coe) '.mat'];
   elseif clustermethod=="angle"
   atominfo_datafilename=['E:\RoDrtest\RoDr_alphastart0-1\model' num2str(modelnum(modeli,1)) '\model' num2str(modelnum(modeli,1)) '_clustermethod_' char(clustermethod) '_Angle' num2str(angle_threshold/pi*180) '.mat'];
   end
   save(atominfo_datafilename,'cluster_binrecord','-v7.3');
   clear cluster_binrecord
   save(atominfo_datafilename,'sample_step','-append');
   clear sample_step
   save(atominfo_datafilename,'cluster_sizerecord','-append');
   clear cluster_sizerecord
   save(atominfo_datafilename,'id_centroids_radius_vxyz_record','-append');
   clear id_centroids_radius_vxyz_record


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