clc;
clear;
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\dataanalysis\');
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\2d_packing\data_analysis_code\functions\');
Foldersdir_pairwise = ('E:\RoDrtest\RoDr_model58_pairwise');
Foldersdir_atominfo = ('E:\RoDrtest\RoDr_alphastart0-1\model58');
%% pairwise
txt_read_pairwise=dir([Foldersdir_pairwise '\*.txt']);
[txt_read_pairwise]=sortnamebysequence(txt_read_pairwise);
[Rtxt_pair,Ctxt_pair]=size(txt_read_pairwise);
%% atom
txt_read_atominfo=dir([Foldersdir_atominfo '\*.csv']);
[txt_read_atominfo]=sortnamebysequence(txt_read_atominfo);
[Rtxt_atom,Ctxt_atom]=size(txt_read_atominfo);
% radius=0.002;
%% Prepare atom and contact pair info
% sample_step=[300:20:Rtxt_atom]';
sample_step=[400]';
atomid_clusterid=[];
cluster_id_size=[];
atomid_clustersize=[];
cluster_x=[];
max_cluster_y=[];
cluster_z=[];
cluster_count=0;
%%
for i=1:size(sample_step,1)
    %% atom and pair information
    % cd([Foldersdir_pairwise '\' subFolders_pairwise(sfipair).name]);
    [~,x1,y1,z1,x2,y2,z2,id1,id2,~,fntx,fnty,fntz,fnx,fny,fnz,ftx,fty,ftz]=textread( ...
        [Foldersdir_pairwise '\' txt_read_pairwise(sample_step(i)).name],'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','headerlines',9);
    % % % [contacts,xi,yi,~,xj,yj,~,id1,id2,id3,fntx,fnty,fntz,fnx,fny,fnz,ftx,fty,ftz]
    pairfn=sqrt(fnx.^2+fny.^2+fnz.^2);
    pairdist_directcal=sqrt((x1-x2).^2+(y1-y2).^2+(z1-z2).^2);
    % pairfnt=sqrt(fntx.^2+fnty.^2+fntz.^2);
    % pairft=sqrt(ftx.^2+fty.^2+ftz.^2);

    [num,txt,raw]=xlsread( ...
        [Foldersdir_atominfo '\' txt_read_atominfo(sample_step(i)).name],1);
    id=num(:,1);
    x=num(:,18);
    y=num(:,19);
    z=num(:,20);
    radius=num(:,12);
    vx=num(:,6);
    vy=num(:,7);
    vz=num(:,8);
    VliqonVatom=num(:,17);
    vel=sqrt(vx.^2+vy.^2+vz.^2);
%     [id,type,x,y,z,ix,iy,iz,vx,vy,vz,fx,fy,fz,omegax,omegay,omegaz,radius,VliqonVatom]=textread( ...
%         txt_read_atominfo(samplei(1,i)).name,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','headerlines',9);
    %%

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
    %% Network analysis gap factor step1
%     tic
%     louvain_repeats=20;
%     resolution=0.5:0.4:2;
% %      :0.1:1;
%      rho=1;
%      [g,com,Q,Qmax,rc,centroids1,Stotal,Stotal1,n]=Network_analysis( ...
%          centroids,pairid_force_dis_sumradij,louvain_repeats,resolution,rho,Bij,Wij,Fij);
% % % === Function: Network_analysis, input explanation ===== Network_analysis(centroids,pairid_force_dis,louvain_repeats,resolution,rho)
% % %     g = gap_factor
% % %     com = particle id in each community
% % %     Q = weithting function
% % %     rc = size(quantity of particles) of each community
% % %     Stotal1 = max Q community allocation
% % %     n = total community numbers
%         gap_factor(i,1:length(resolution))=g;
% %      plot(resolution,gap_factor(i,1:length(resolution)),'Color',color{1,file(1,i)});
% %      hold on;
%      toc
     %% Network analysis gap factor step2
%      % % find max gap factor and further obtain information; 
%      % % this can further be transfered to tool box function
%      optimal_resolution(i,1)=resolution(1,gap_factor==max(gap_factor));
%      optimal_resolution(i,2)=find(gap_factor==max(gap_factor));%corresponding column
%      color_rc=zeros(length(centroids),1);
%      color_fave=zeros(length(centroids),1);
%      for comi=1:n
%          member=cell2mat(com(comi,optimal_resolution(i,2)));
%          color_rc(member(:,1),1)=rc(comi,optimal_resolution(i,2));
%          color_fave(member(:,1),1)=ones(length(member),1)*aveforce_in_community(member,Fij);
%          color_size(member(:,1),1)=ones(length(member),1)*length(member);
%      end
%      size=25;
% %      scatter(centroids1(:,1),centroids(:,2),size,color_fave,'filled');%scatter fave in each c
%      scatter(centroids1(:,1),centroids(:,2),size,color_size,'filled');%scatter size in each c
%      scatter(centroids1(:,1),centroids(:,2),size,color_rc,'filled');%scatter rc in each c
%      colorbar;
%      colormap jet;
%      axis equal;
%      xlim([0 0.1]);

  %% Coordination number
crit_l=1.0001;
    [Zc,Cn_ep,Cn_ep_neibid]=Z_positionbase_rotationdrum(radius,[x,y,z],crit_l,id);
% figure(1)
% scatter3(x,y,z,[],Cn_ep,"filled");
% title('Coordination number')
% figure(2)
% scatter3(x,y,z,[],vel,"filled");
% title('velocity')
% figure(3)
% plot(Cn_ep,vel,'k.');
% xlabel('Coordination number');
% ylabel('Particle velocity');
% xlim([-1 9]);
% figure(4)
% histogram(Cn_ep,'Normalization','probability');
% xlabel('Coordination number');
% ylabel('Probability');
% xlim([-1 9]);

%% import cell format velocity
if i==1
    load([Foldersdir_atominfo '\model58_velmatrix1529frames']);
    figure(1)
    contourf(Xpf,Zpf,velxincell_step1_total_ave)
    for vmi=1:size(id,1)
    xtmp=x(id==vmi);
    ztmp=z(id==vmi);
    vm_x_index=find(Xpf(:,1)>xtmp,1)-1;%locate the index of bigger one for the interval [smaller,bigger]
    vm_z_index=find(Zpf(1,:)>ztmp,1)-1;
    vm_ep(vmi,1:2)=[velxincell_step1_total_ave(vm_x_index,vm_z_index) velzincell_step1_total_ave(vm_x_index,vm_z_index)];%mean velocity for each particle
    end
else
    continue
end
%% import DAOR and obtain n_flow
% load([Foldersdir_atominfo '\model55_DAOR_gof']);
% DAOR_steady=mean(DAOR(revolution(:,1)>100,1));
% v_flow_angle=DAOR_steady-pi;
% n_flow=[cos(v_flow_angle) 0 sin(v_flow_angle)]; %the unit vector of average velocity in the flow region, which corresponds to to DAOR

%% velocity correlation method: correlation between each particle and the unit vector of DAOR
    n_ep=zeros(size(id,1),1);
    for idi=1:size(id,1)
        v_fluctuation=[vx(idi)-vm_ep(idi,1) vy(idi) vz(idi)-vm_ep(idi,2)];
        n_ep(idi,1:3)=v_fluctuation;
%         /norm(v_fluctuation);
    end

    pair_vel_angle_clustermark=zeros(size(fnx,1),1);
    for pairi=1:size(fnx,1)
        index_id1=find(id==id1(pairi));
        index_id2=find(id==id2(pairi));
        n1=[n_ep(index_id1,1) n_ep(index_id1,2) n_ep(index_id1,3)];
        n2=[n_ep(index_id2,1) n_ep(index_id2,2) n_ep(index_id2,3)];
        pair_vel_angle_clustermark(pairi,1) = atan2(norm(cross(n1,n2)),dot(n1,n2));
        if pair_vel_angle_clustermark(pairi,1)<=pi/6 && vx(id==id1(pairi))<0 && vz(id==id1(pairi))<0
            pair_vel_angle_clustermark(pairi,2)=1;
        else
            continue
        end
    end


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

%% cluster identification based on graph theory
pairid=[id1,id2];
cluster_pairid=pairid(pair_vel_angle_clustermark(:,2)==1,:);
G=graph(cluster_pairid(:,1),cluster_pairid(:,2));
figure(6)
plot(G)
%% remove single & double particle clusters
[bin,binsize] = conncomp(G,'Type','weak');
idx = binsize(bin) == max(binsize);
SG = subgraph(G, idx);
figure(7)
plot(SG)
%% Apply biconncomp function
bincell = biconncomp(SG, 'OutputForm', 'cell');
n = length(bincell);

% plot each cluster
figure(8)
for ii = 1:n
    subplot(20,20,ii)
    plot(subgraph(G, bincell{ii}), 'NodeLabel', bincell{ii});
end

%% release cluster id from bincell
for clusi=1:n
    if size(bincell{1,clusi}',1)>2 && size(bincell{1,clusi}',1)<5000
    cluster_count=cluster_count+1;
    cluster_id_size(cluster_count,1:2)=[clusi size(bincell{1,clusi}',1)];
    cluster_id=ones(size(bincell{1,clusi}',1),1)*clusi;
    atomid_clusterid=[atomid_clusterid;bincell{1,clusi}',cluster_id];
    atomid_clustersize=[atomid_clustersize;bincell{1,clusi}',ones(size(bincell{1,clusi}',1),1)*size(bincell{1,clusi}',1)];
    else
        continue
    end
end
% clusterid=unique(clusterid(:,1),"rows");
for clusj=1:size(atomid_clusterid,1)
cluster_x(clusj,1)=x(id==atomid_clusterid(clusj,1));
cluster_y(clusj,1)=y(id==atomid_clusterid(clusj,1));
cluster_z(clusj,1)=z(id==atomid_clusterid(clusj,1));
end

figure(8)
scatter3(cluster_x,cluster_y,cluster_z,[],atomid_clusterid(:,2),"filled");
title('cluster id')
figure(9)
scatter3(cluster_x,cluster_y,cluster_z,[],atomid_clustersize(:,2),"filled");
title('cluster size')

%% plot single clusters
figure(10)
max_cluster_x=[];
max_cluster_y=[];
max_cluster_z=[];
tmp_atomid_clustersize=atomid_clustersize;
tmp_atomid_clustersize(tmp_atomid_clustersize(:,2)>300,:)=[];
% the id of particles in the largest cluster
[~,id_maxcluster]=ismember(tmp_atomid_clustersize(tmp_atomid_clustersize(:,2)==max(tmp_atomid_clustersize(:,2)),1),id);
max_cluster_x(:,1)=x(id_maxcluster);
max_cluster_y(:,1)=y(id_maxcluster);
max_cluster_z(:,1)=z(id_maxcluster);
scatter3(max_cluster_x,max_cluster_y,max_cluster_z,"MarkerFaceColor",[0 1 0]);

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

