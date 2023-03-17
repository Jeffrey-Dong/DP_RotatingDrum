clc;
clear;
%%
pivdata=load('D:\EXP RoDrtest\Exp\for test\PIVlab_doublemask_sharpened_pass6432.mat');
cd('D:\EXP RoDrtest\Exp\for test')
xextacu=zeros(size(cell2mat(pivdata.x(1,1)),1),size(cell2mat(pivdata.x(1,1)),2));%acu=accumulate
yextacu=zeros(size(cell2mat(pivdata.x(1,1)),1),size(cell2mat(pivdata.x(1,1)),2));
voriginal_extacu=zeros(size(cell2mat(pivdata.x(1,1)),1),size(cell2mat(pivdata.x(1,1)),2));
uoriginal_extacu=zeros(size(cell2mat(pivdata.x(1,1)),1),size(cell2mat(pivdata.x(1,1)),2));

for i=1:size(pivdata.x,1)
    xext=cell2mat(pivdata.x(i,1));
    yext=cell2mat(pivdata.y(i,1));
    voriginal_ext=cell2mat(pivdata.v_original(i,1));
    uoriginal_ext=cell2mat(pivdata.u_original(i,1));
    vel=sqrt(voriginal_ext.^2+voriginal_ext.^2);
%     contour(xext,yext,vel);
%     figure(i+2)
%     quiver(xext,yext,uoriginal_ext,voriginal_ext,3);
%     set(gca, 'ydir', 'reverse');
%     axis off
   xextacu=xextacu+xext;
   yextacu=yextacu+yext;
   voriginal_extacu=voriginal_extacu+voriginal_ext;
   uoriginal_extacu=uoriginal_extacu+uoriginal_ext;
end
   %%
   xextacu=xextacu/size(pivdata.x,1);
   yextacu=yextacu/size(pivdata.x,1);
   voriginal_extacu=voriginal_extacu/size(pivdata.x,1);
   uoriginal_extacu=uoriginal_extacu/size(pivdata.x,1);
    figure(4);
    orifigure=imread('D:\EXP RoDrtest\Exp\for test\original_photos\C0055.MP4_20220225_094202.186.jpg');
    imshow(orifigure);
    hold on
    quiver(xextacu,yextacu,uoriginal_extacu,voriginal_extacu,3);
    set(gca, 'ydir', 'reverse');
    axis off
    
    %% -----plot rotated data according to AOR point---------------
    coordtmp=[reshape(xextacu,[],1) (reshape(yextacu,[],1)-1080)*-1];% (data-1080)*-1 means a reverse in y direction
    rota_datacoord= rotation(coordtmp,[coord_maxanglex,coord_maxangley],rotation_angle);
    after_rotx=reshape(rota_datacoord(:,1),size(xextacu,1),size(xextacu,2));
    after_roty=(reshape(rota_datacoord(:,2),size(xextacu,1),size(xextacu,2)));
    figure(7);quiver(after_rotx,after_roty,uoriginal_extacu,voriginal_extacu,3);
    hold on; plot(coord_maxanglex,coord_maxangley,'go','MarkerSize',10,'LineWidth',1);
%     set(gca, 'ydir', 'reverse');
    
    %% -----extract velocity profile--------------
    d_p=5.4;
    vel_matrix2col_v=reshape(voriginal_extacu,[],1);
    vel_dircorr=ones(size(vel_matrix2col_v,1),1);%velocity direction correction
    vel_dircorr(vel_matrix2col_v(:,1)<0)=-1;
    vel_totacu=reshape(sqrt(voriginal_extacu.^2+uoriginal_extacu.^2),[],1).*vel_dircorr;
    vel_isrx_index=intersect(find(rota_datacoord(:,1)<S1(1)),find(rota_datacoord(:,1)>S2(1)));%extract data in sampling region
    vel_isry_index=intersect(find(rota_datacoord(:,2)<S1(2)),find(rota_datacoord(:,2)>S3(2)));%extract data in sampling region
    vel_isr_index=intersect(vel_isrx_index,vel_isry_index);
    vel_isr_coord=rota_datacoord(vel_isr_index,:);
    vel_isr=vel_totacu(vel_isr_index,:);
%     figure(8);quiver(after_rotx,after_roty,uoriginal_extacu,voriginal_extacu,3);
    figure(9)
    subplot(2,1,1);
    scatter(vel_isr_coord(:,2),vel_isr_coord(:,1),[],vel_isr);
    xlim([50 550]);
    axis equal;
    wsv=6;%window size vertical; unit = pixels
    wsh=S1(1)-S2(1);
    movingsteps=ceil((S1(2)-S3(2))/wsv);
    vel_el=zeros(movingsteps,1);%velocity each layer
    middleline=zeros(movingsteps,1);
    for velisri=1:movingsteps
        upperbound=S3(2)+velisri*wsv;
        bottombound=S3(2)+(velisri-1)*wsv;
        middleline(velisri,1)=(upperbound+bottombound)*0.5;
        veltmp_extract=[];
        veltmp_extract=vel_isr(intersect(find(vel_isr_coord(:,2)>bottombound),find(vel_isr_coord(:,2)<upperbound)));
        veltmp_extract(isnan(veltmp_extract))=[];
        vel_el(velisri,1)=mean(veltmp_extract);
          
    end
    subplot(2,1,2);
    plot(middleline,vel_el);
    
    