% analysis electro_analysis for multicases
% clc;
% clear;
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\2d_packing\data_analysis_code\functions\');
%% Read excel file
% d_file=(['C:\Users\jeffr\Dropbox (Sydney Uni)\projects\2d_packing\electro conductivity\results\6Jul2020\all pair contact gra01\']);
% % d_file=(['C:\Users\jeffr\Dropbox (Sydney Uni)\projects\2d_packing\electro conductivity\results\6Jul2020\model46_gra01\PCP calculate using pair contact\']);
% % d_file=(['D:\electro_conduct\tmp\model46_gra01\model46 all pair contact\']);
% % d_file=(['C:\Users\jeffr\Dropbox (Sydney Uni)\projects\2d_packing\electro conductivity\results\22Jun2020\hydrostatic pressure 395Pa\all pcp by contact pair\']);
% d=dir([d_file '*.xlsx']);
% % % sort name by step sequence
% [d]=sortnamebysequence(d);
% file=1:size(d,1);
%% load color and marker
% tmpcolor1=load('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\2d_packing\data_analysis_code\customise_colorbar\color');
% tmpcolor2=struct2cell(tmpcolor1);
% custo_colormap=tmpcolor2{1,1};
custo_colormap={'k','r','b','g','m','c'}';
matlab_marker_list=load('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\2d_packing\data_analysis_code\functions\matlab_marker_list');
matlab_marker_list=struct2cell(matlab_marker_list);
matlab_marker_list=matlab_marker_list{1,1};

%%
gra_dir_cell={'','_gra22-5deg','_gra45deg'};
Mstress=[4,5,6,7]';
Sratio={'1','1-01','1-03','1-06','1-1','1-2','1-3'};
Ns={'1','3','5','6','7','8','9'};
Ns_data=[0.1 0.3 0.5 0.6 0.7 0.8 0.9]';
GR={'2e-6','2e-7','2e-8','2e-9','2e-10'};
GR_data=[2E-6;2E-7;2E-8;2E-9;2e-10];
fri={'fri0','fri0-01','fri0-05','fri0-1'};

%%
for viii=1%gra_dir; gravity direction
for vii=2%Mstress
for vi=6%Size ratio
for v=3%Ns
for iv=4%GR
for iii=1%fri
    caseNo=(1:60)';
for ii=70
%% Read txt file
% _set-sw_groeve1
% _settlethenswell
% _iniD0-0025
% _iniD_0_0025
% _CoeRes0-6
% _shxA005d
% _shsameD4-005
% _Sratio1-06
%  '_relax'
% '_gra45deg'
% '_gra22-5deg'
%-----N5000 binary
% read_filename1=(['D:\D_d4\D_d1-6_' num2str(caseNo(ii,1)) '_Ns0' char(Ns(1,iiiii)) '_CoeRes0-6_GR' char(GR(1,iiii)) '_' char(fri(1,iii)) '_set-sw_groeve1_shsameD4-005_Sratio' char(Sratio(1,iiiiii))]);
%-----N5000 binary with different degree of Mean stress
% read_filename1=(['D:\D_d4\D_d1-6_' num2str(caseNo(ii,1)) '_Ns0' char(Ns(1,v)) '_CoeRes0-6_GR' char(GR(1,iv)) '_' char(fri(1,iii)) '_set-sw_groeve1_sh0_Sratio' char(Sratio(1,vi)) '_Mstre-' num2str(Mstress(vii,1)) char(gra_dir_cell(1,viii))]);
% sameD4-005
%-----N4978 fullycrystallized different degree of Mean stress
% read_filename1=(['D:\D_d4\D_d1-6_1_fullycrys_void0_fri02_Mstre-5']);
% _fri02
%-----N5000 mono
% read_filename1=(['D:\D_d4\D_d1-6_' num2str(caseNo(ii,1)) '_CoeRes0-6_GR' char(GR(1,iiii)) '_' char(fri(1,iii)) '_set-sw_groeve1_shsameD4-005']);
%-----N700 mono
% read_filename1=(['D:\D_d4\D_d4_' num2str(caseNo(ii,1)) '_CoeRes0-6_GR' char(GR(1,iiii)) '_' char(fri(1,iii)) '_set-sw_groeve1_shxA005d']);
%-----N1109 binary
% read_filename1=(['D:\D_d4\D_d3_3-2_' num2str(caseNo(ii,1)) '_Ns03_CoeRes0-6_GR2e-9_fri0_set-sw_groeve1_shsameD4-005_Sratio1-06']);
%-----N2065 mono
% read_filename1=(['D:\D_d4\D_d2-3_' num2str(caseNo(ii,1)) '_CoeRes0-6_GR' char(GR(1,iiii)) '_' char(fri(1,iii)) '_set-sw_groeve1_shxA005d']);
%-----for test
% read_filename1=(['D:\tmptest\D_d4_' num2str(caseNo(ii,1)) '_CoeRes0-6_GR' char(GR(1,iiii)) '_' char(fri(1,iii)) '_set-sw_groeve1_shxA02d']);
% read_filename1=(['C:\Users\jeffr\Dropbox (Sydney Uni)\projects\2d_packing\electro conductivity\results\25Mar2021\N50']);

% datominfo=dir([read_filename1 '_atominfo\*.txt']);
% dpairwise=dir([read_filename1 '_pairwise\*.txt']);
% dwallwise=dir([read_filename1 '_wallwise\*.txt']);
% [datominfo]=sortnamebysequence(datominfo);
% [dpairwise]=sortnamebysequence(dpairwise);
% [dwallwise]=sortnamebysequence(dwallwise);
d=dir([read_filename1 '_atominfo\*.txt']);
[d]=sortnamebysequence(d);
 

%% reference filenumber
% ref1=load('D:\D_d4\HydroP10000_fileNo_d1-6mono');
% ref2=struct2cell(ref1);
% ref3=ref2{1,1};
% ref_filenum=cell2mat(ref3(iiiiii,iiiii));

%%
% for i=ref_filenum(ii,1)
for i=size(d,1)
% [num,txt,raw]=xlsread(['C:\Users\jeffr\Dropbox (Sydney Uni)\projects\2d_packing\Xu\19Jun220\WonH' num2str(WonH(1,i)) 'AveCN' num2str(CN(1,i)) '.xlsx'],2);
%% Read excel file
% [num1,num2,num3,id,x,y,radius,id1,id2,force,PCP,CN,step_curr]=read_excel_data(d_file,d);

%% Read 
cd([read_filename1 '_atominfo']);
[id,~,x,y,~,~,~,~,~,~,~,~,~,~,~,~,~,radius]=textread(d(file(1,i)).name,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','headerlines',9);
cd([read_filename1 '_pairwise']);
[~,~,~,~,~,~,~,id1,id2,~,~,~,~,fnx,fny,~,~,~,~]=textread(d(file(1,i)).name,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','headerlines',9);
% [contacts,xi,yi,~,xj,yj,~,id1,id2,~,~,~,~,fnx,fny,~,~,~,~]
cd([read_filename1 '_wallwise']);
[~,~,~,id_wallcontact,~,~,~,wfnx,wfny,~]=textread(d(file(1,i)).name,'%f%f%f%f%f%f%f%f%f%f','headerlines',9);

force=sqrt(fnx.^2+fny.^2);

%% Prepare info
[id_corre,centroids,radius,pairid_force_dis,Fij,Wij,Bij,ave_force] = correct_id_and_obtain_matrix(id,x,y,radius,id1,id2,force);

%% contact area based on hertzian contact (1 to export)
% Ei=1.1e11;%Young's modulus
% Ej=Ei;
% nui=0.2;%Poisson's Ratio
% nuj=nui;
% Estar=(Ei*Ej)/(Ei*(1-nuj^2)+Ej*(1*nui^2));
% for cta=1:size(pairid_force_dis,1)
%     radi=radius(pairid_force_dis(cta,1));
%     radj=radius(pairid_force_dis(cta,2));
% Radiistar=radi*radj/(radi+radj);
% end
% C_a=pi*(3.*pairid_force_dis(:,3).*Radiistar/(4*Estar)).^(2/3);%contact area
% % pair_resistance=(1./(0.1*(pairid_force_dis(:,3)./1681).^(-2/3)));
% mean_C_a(ii)=mean(C_a);

%% Coordination Number (2 to export)
% % Average mechanical coordination number
% Z_m(ii,1)=size(pairid_force_dis,1)*2/size(centroids,1);
% CN_me=[];
% for CNi=1:size(centroids,1)
%     CNa=pairid_force_dis(pairid_force_dis(:,1)==CNi,:);
%     CNb=pairid_force_dis(pairid_force_dis(:,2)==CNi,:);
%     CNcombine=[CNa;CNb];
%     CN_me(CNi,1)=size(CNcombine,1);%mechanical coodination number of each particle
%     
% end
%     cell_CN_me{ii,1}=CN_me;

% % plot particles
% scatter(centroids(:,1),centroids(:,2),10,CN_me,'filled');%scatter fave in each c
% 
%     scatter_plot=gca;
%     axis equal
%     xlim([0 0.1]);
%     ylim([0 0.1]);
%     axis off;
%     set(gcf,'color',[1 1 1]);
% %     set(gcf, 'InvertHardcopy', 'off')
%     colormap parula(7);
%     Color_Bar=colorbar;
%     Color_Bar.Label.String='Coordination Number';
%     Color_Bar.FontSize=12;
%     Color_Bar.Label.FontSize=18;
% 
%     hold on
    
%% bar graph of contact quantity and force
% PDF_binwidth=0.25;
% figure(1)
% f_histo(vi,ii)=histogram(pairid_force_dis(:,3)./mean(pairid_force_dis(:,3)),'normalization','probability','BinWidth',PDF_binwidth,'FaceColor',custo_colormap{vi});
% Values=f_histo(vi,ii).Values;
% BinNumber=f_histo(vi,ii).BinWidth:f_histo(vi,ii).BinWidth:f_histo(vi,ii).BinLimits(1,2);
% close figure 1
% figure(2)
% 
% f_pdf(vi,ii)=plot(BinNumber,Values,[custo_colormap{vi} ]);
% % matlab_marker_list{vi}
% set(gca, 'YScale', 'log');
% % set(gca, 'XScale', 'log')
% ylim([1e-5 2]);
% xlim([0 20]);
% xlabel('F_n/<F_n>');
% ylabel('Probability');
% % legend([f_pdf(1,1),f_pdf(1,2),f_pdf(5,1),f_pdf(6,1)],'Bina; \eta=1.1; C_s=0.5; Fri 0','Bina; \eta=1.2; C_s=0.5 Fri 0' ,'mono; Fully crystal; Fri 0', 'mono; Fully crystal; Fri 0.2','fontsize',12);

% b=a.BinEdges;
% b=b(2:length(b));
% CN3d=zeros(1,length(a.Values));
% CN3d(1,:)=PCP(i,1);
% figure(2)
% plot3(b,CN3d,a.Values);
% hold on;


% figure(2);
% a.Normalization='probability';
% hold on;

%% cdf
% figure(2);
% b=cdfplot(num2(:,3));
% hold on;

%% pdf
% [f,xi,bw]=ksdensity(pairid_force_dis(:,3)./mean(pairid_force_dis(:,3)));
% f_distribution(viii,ii)=plot(xi,f,[custo_colormap{viii}]);
% hold on;
% set(gca, 'YScale', 'log')
% ylim([1e-6 2])
% xlim([0 20]);
% figure(3);
%     path(path,'C:\Users\jeffr\Dropbox (Sydney Uni)\projects\2d_packing\data_analysis_code\functions');
%     [pcp,vno,V,avno,oripar,orivec1]=crystallization(centroids,radius,0.9,0.001);

% CN3d=zeros(1,length(f));
% CN3d(1,:)=pcp;
% plot3(xi,CN3d,f);
% hold on;

%% packing fraction in total (8 to export)
% % addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\2d_packing\data_analysis_code\functions\');
% % [total_Packing_fra]=packing_fraction_intotal(centroids,radius,0.01);
% % fcover=sum(abs(wfny(wfny<0)));
% % if total_Packing_fra>0.8
% % if fcover>500
% % 
%     [FaT,StrT,FaT_angF1x,FaT_angF1y,StrT_angF1x,StrT_angF1y,FaT_Prin,StrT_Prin]=fabric_tensor(pairid_force_dis,centroids,radius);
% % 
%       Fa1_Fa2_A(ii,:)=FaT_Prin;
%       Str1_Str2_A(ii,:)=StrT_Prin;
% %     FaTcell{i,ii}=FaT;
% %     StrTcell{i,ii}=StrT;
% %     StrT_Prin_cell{i,ii}=StrT_Prin;
% % 
%       hydro_P(ii,1)=(StrT(1,1)+StrT(2,2))/2;
% % 
% %     if hydro_P>50000000
% %     filenum(ii,1)=i;
% %     mean_stress(ii,1)=hydro_P(ii,1);
% % 
%       [pcp(ii,1),vno,V,avno,oripar,orivec1,id_cryornot,pidfdis,id_cry_phi6l]=crystallization(centroids,radius,0.75,min(radius)*0.4,pairid_force_dis);
%       [pcp_m(ii,1),vno,V,avno,oripar,orivec1,id_cryornot,pidfdis,id_cry_phi6l]=crystallization(centroids,radius,0.75,0,pairid_force_dis);
% %     
% %     radius_small=radius(radius==min(radius));
% %     volfra_s(ii,1)=sum(pi*radius_small.^2)/sum(pi*radius.^2);
% %     packingfra(ii,1)=total_Packing_fra;
% % 
%       mean_d(ii,1)=mean(radius)*2;
%       FaT_angleF1x(ii,1)=FaT_angF1x;
%       StrT_angleF1x(ii,1)=StrT_angF1x;
% %     
% %     f_cover(ii,1)=fcover;
% %     f_cover_id=id_wallcontact(wfny<0)-1;
% %     break
% %     else
% %         continue
% %     end
% % else
% %     continue
% % end

%% packing fraction along x and y
%     centroids(:,1)=num1(:,2);
%     centroids(:,2)=num1(:,3);
%     radius=num1(:,4);
%     path(path,'C:\Users\jeffr\Dropbox (Sydney Uni)\projects\2d_packing\data_analysis_code\functions');
%     jeff_yellow=[0.9290, 0.6940, 0.1250];
%     jeff_orange=[0.8500, 0.3250, 0.0980];
%     color={jeff_orange 'blue' 'green' 'yellow' 'red' 'black'};
%     [Px,xonX,Py,yonY]=packingfraction_layerly(centroids,radius);
%     figure(1);
%     plot(xonX,Px,'Color',color{1,file(1,i)});
%     hold on;
%     figure(2);
%     plot(Py,yonY,'Color',color{1,file(1,i)});
%     hold on;
   
%% pair distribution function
%     [gr,mel]=pair_distribution(centroids,radius);
%     plot(mel,gr);
%     grtmp(ii,:)=gr;
%     meltmp(ii,:)=mel;
    % % plot 3d 
%     CN3d=zeros(1,length(mel));
%     CN3d(1,:)=pcp;
%     figure(2)
%     plot3(mel,CN3d,gr,'Color',color{1,file(1,i)});
%     hold on;

    
%% crystallization
% %     addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\2d_packing\data_analysis_code\functions');
% if i>=60
%     [pcp(i,1),vno,V,avno,oripar,orivec1,id_cryornot,pidfdis,id_cry_phi6l,phi6l,phil]=crystallization(centroids,radius,0.75,0,pairid_force_dis);
% else
%     pcp(i,ii)=0;
% end

% % %---------convert radian to degree----------
%     avnodeg=rad2deg(avno);
%     avnodeg(avnodeg<0,1)=avnodeg(avnodeg<0,1)+60;
%     avnodeg=avnodeg-avnodeg;

% % % ------------plot crystallization through voronoi------------
% figure(1)
%     crys_plot=voronoi(centroids(:,1),centroids(:,2));
%     hold on;
%     % plot crystallized voronoi cells
% 
%     patch('Faces',vno,'Vertices',V,'FaceVertexCData',avnodeg,'FaceColor','flat');
%     hold on;
%     %  plot arrow based on angle
% %     quiver(oripar(:,1),oripar(:,2),orivec1(:,1),orivec1(:,2),'Color','r');
%     
%     xlim([min(centroids(:,1))-max(radius) max(centroids(:,1))+max(radius)]);
%     ylim([min(centroids(:,2))-max(radius) max(centroids(:,2))+max(radius)]);
% %     pbaspect([1 1 1]);
%     axis equal
%     axis off
%     xlim([0 0.1]);
%     ylim([0 0.1]);
% %     hold on
%     colorbar;
%     colormap gray;
%     figureall=gca;
%     Color_Bar=colorbar;
%     Color_Bar.Limits=[0,60];
%     Color_Bar.TickLabels={'0\circ','10\circ','20\circ','30\circ','40\circ','50\circ','60\circ'};
%     Color_Bar.Label.String='Grain Angle';
%     Color_Bar.Label.Position=[4.081249952316284,29.182908583928835,0];
%     Color_Bar.Label.Rotation=-90;
%     figureall.FontSize=18;
%     xlim([0 0.1]);
%     ylim([0 0.1]);

% % % ------------plot crystallization through scatter------------
% figure(1)
% crange=(0:size(custo_colormap,1))';
% cry_color_cus=[];
% for ci=1:size(avnodeg,1)
%     for cj=1:size(custo_colormap,1)
%         if and(avnodeg(ci,1)>=crange(cj,1),avnodeg(ci,1)<crange(cj+1,1))
%         cry_color_cus(ci,:)=custo_colormap(cj,:);
%         end
%     end
% end
% scattersize=10;
% cry_color_cus=cry_color_cus./255;
% 
% ctmp=linspace(0,60,61);
% cry_color_hsv=linspace(0,60,120);
% 
% scatter(centroids(id_cry_phi6l,1),centroids(id_cry_phi6l,2),scattersize,avnodeg,'filled');%scatter fave in each c
%     scatter_plot=gca;
%     axis equal
%     xlim([0 0.1]);
%     ylim([0 0.1]);
%     axis off;
%     set(gcf,'color',[1 1 1]);
%     set(gcf, 'InvertHardcopy', 'off')
% 
% %---below only for hsv
%     colormap hsv(12);
% %     figureall=gca;
%     Color_Bar=colorbar;
%     Color_Bar.Limits=[0,60];
%     Color_Bar.TickDirection='out';
% %     Color_Bar.Location='Southoutside';
% %     Color_Bar.Position=[0.036363636363636,0.393650793650794,0.911255411255411,0.052380952380956];
%     Color_Bar.TickLabels={'0\circ','10\circ','20\circ','30\circ','40\circ','50\circ','60\circ'};
%     Color_Bar.Label.String='Grain Angle (\alpha)';
%     Color_Bar.Label.FontSize=18;
% %     Color_Bar.Label.Position=[30.40038689934942,-1.56,0];
% %     Color_Bar.Label.Rotation=0;
% %     figureall.FontSize=22;
% %     figureall.FontWeight='bold';
  

% % % % save figure
% crystallization plot
%     figurename=['D_d1-6_' num2str(caseNo(ii,1)) '_Ns0' char(Ns(1,iiiii)) '_CoeRes0-6_GR' char(GR(1,iiii)) '_' char(fri(1,iii)) '_set-sw_groeve1_shsameD4-005_Sratio' char(Sratio(1,iiiiii))];
%     cd('D:\D_d4\figures_HydroP10000_fileNo_d1-6bina\original hsv');
%     saveas(gcf,[figurename '_cry' '.png']);
%     close figure 1

% % angle distribution plot
%     figure(3)
%     polarhistogram(deg2rad(avnodeg),60,'EdgeColor','b','Normalization','probability');
%     set(gca,'FontSize',18);
%     thetalim([0 60]);
%     rlim([0 1]);
%     polar=gca;
%     polar.RAxis.Label.String = 'Probability Density';
%     polar.ThetaAxis.Label.String = 'Grain Angle';
%     polar.ThetaAxis.Label.Position=[45.837322456043125,1.220819789876846,0];
%     polar.ThetaTickLabel={'0\circ','15\circ','30\circ','45\circ','60\circ'};
%     saveas(gcf,[figurename '_avno' '.png']);
%     close figure 3

% % --------------grain size and No analysis
% [grainNo(ii,1),mean_grainsize(ii,1)]=cryedgrain_analysis(pidfdis,centroids);

%% plot circle
%     figure(2)
    [plot_circle]=plot_circles(centroids,radius);
%     axis equal;
%     axis off;
%     xlim([0 0.1]);
%     ylim([0 0.1]);
%     
%     % hold on;
% 
%     saveas(gcf,[figurename '_circle' '.png']);
%     close figure 2


%% force network
normalise_para=mean(pairid_force_dis(:,3));
% %   plot overall force chain
% %   pairid_force_dis=num2(:,3);
%     cd(['C:\Users\jeffr\Dropbox (Sydney Uni)\projects\2d_packing\electro conductivity\results\6Jul2020\all pair contact gra01\']);

%     figure(5)
    force_chain_network=plot_forcechain(pairid_force_dis,normalise_para,centroids,'k');
%     axis equal;
%     axis off;
%     xlim([0 0.1]);
%     ylim([0 0.1]);

%     title([name_split{1,1} '; ' name_split{2,1} '; ' name_split{3,1}]);
%     saveas(force_chain_network,[name_split{1,1} '_' name_split{2,1} '_' name_split{3,1} '.png']);
%     close figure 1
%     hold on;

% %  find main force chain
%     [mfc_pair,fc_pair,mfc_path,fc_path,mfc_aft,fc_aft,mfc_avef,fc_avef,dist]=find_main_forcechain(centroids,pairid_force_dis,posi_ele,nega_ele,Fij,normalise_para);
 
% % plot the single main force chain(largest force chain);can also used to
% % find percolation path
% % mfc_path=mfc_path;
%     for mfci=1:(length(mfc_path)-1)
%         a=[centroids(mfc_path(mfci,1),1),centroids(mfc_path(mfci+1,1),1)];
%         b=[centroids(mfc_path(mfci,1),2),centroids(mfc_path(mfci+1,1),2)];
%         plot(a,b,'r','LineWidth',Fij(mfc_path(mfci,1),(mfc_path(mfci+1,1)))/normalise_para);
%         hold on;
%     end

% % partition edge by angle
% force_thre=9;
% angle_thre=pi/8;
% angle_mainf=pairid_force_dis(pairid_force_dis(:,3)>=force_thre,:);
% % first partition, put all edge with similar angles into one cluster
% angle_mainf=sortrows(angle_mainf,3,'descend');

% % make sure id(first column)<id(second column)
% for amfi=1:size(angle_mainf,1)
%     if angle_mainf(amfi,1)>angle_mainf(amfi,2)
%         angle_mainf(amfi,[1 2])=angle_mainf(amfi,[2 1]);
%     else
%         continue
%     end
%         
% end

%% force chain extraction
% amftmp=angle_mainf;%angle mainf tmp
% amftmp1=angle_mainf;
% count_chain=0;
% for angi=1:300
% if ~isempty(amftmp)
% % % partition through (1)angle,(2)connectivity
%     ang1=amftmp(1,5);
%     if (pi-angle_thre)>ang1>(angle_thre)
%         angran=[ang1-angle_thre ang1+angle_thre];%angle range of the target force pair +- 45degree(pi/4)
%         cluster1=amftmp(and(amftmp(:,5)>angran(1,1),amftmp(:,5)<angran(1,2)),:);%(1)angle
%         
%         Gtmp=graph(cluster1(:,1),cluster1(:,2));%(2)connectivity
%         bins = conncomp(Gtmp)';%find connected nodes in Gtmp
%         idoc=find(bins(:,1)==bins(amftmp(1,1),1));%id all particles in one chain
%         tmp1=nchoosek(idoc,2);
%         tmp2=amftmp(ismember(amftmp(:,1:2),tmp1,'rows'),:);
%         pairoc=tmp2(and(tmp2(:,5)>angran(1,1),tmp2(:,5)<angran(1,2)),:);% pair list in one chain
%     elseif ang1>=(pi-angle_thre)
%         angran=[ang1-angle_thre,pi;0,angle_thre-(pi-ang1)];
%         cluster1=amftmp(or(and(amftmp(:,5)>angran(1,1),amftmp(:,5)<angran(1,2)),and(amftmp(:,5)>angran(2,1),amftmp(:,5)<angran(2,2))),:);%(1)angle
%         
%         Gtmp=graph(cluster1(:,1),cluster1(:,2));%(2)connectivity
%         bins = conncomp(Gtmp)';
%         idoc=find(bins(:,1)==bins(amftmp(1,1),1));
%         tmp1=nchoosek(idoc,2);
%         tmp2=amftmp(ismember(amftmp(:,1:2),tmp1,'rows'),:);
%         pairoc=tmp2(or(and(tmp2(:,5)>angran(1,1),tmp2(:,5)<angran(1,2)),and(tmp2(:,5)>angran(2,1),tmp2(:,5)<angran(2,2))),:);
%     elseif ang1<=angle_thre
%         angran=[0,ang1+angle_thre;pi-(angle_thre-ang1),pi];
%         cluster1=amftmp(or(and(amftmp(:,5)>angran(1,1),amftmp(:,5)<angran(1,2)),and(amftmp(:,5)>angran(2,1),amftmp(:,5)<angran(2,2))),:);%(1)angle
%         
%         Gtmp=graph(cluster1(:,1),cluster1(:,2));%(2)connectivity
%         bins = conncomp(Gtmp)';
%         idoc=find(bins(:,1)==bins(amftmp(1,1),1));
%         tmp1=nchoosek(idoc,2);
%         tmp2=amftmp(ismember(amftmp(:,1:2),tmp1,'rows'),:);
%         pairoc=tmp2(or(and(tmp2(:,5)>angran(1,1),tmp2(:,5)<angran(1,2)),and(tmp2(:,5)>angran(2,1),tmp2(:,5)<angran(2,2))),:);
%     end
%     
%     %delete the pair that already constitude a chain
%     amftmp(ismember(amftmp(:,1:2),pairoc(:,1:2),'rows'),:)=[];  
%     
%     
% % % (1)obtain chain angle and length(number of nodes); (2)plot chain; (3)recorde chain pair and id
%     if size(pairoc,1)>1
%         count_chain=count_chain+1;
%         % % obtain chain angle
%         xidoc=centroids(idoc(:,1),1);
%         yidoc=centroids(idoc(:,1),2);
%         p=polyfit(xidoc,yidoc,1);
%         if p(1,1)>=0
%             changle(count_chain,1)=atan(p(1,1));%chain angle
%         else
%             changle(count_chain,1)=atan(p(1,1))+pi;
%         end
%         % obtain chain length(number of nodes)
%         chalen(count_chain,1)=size(idoc,1)-1;%chain length
%                 
%         % % plot the chain
% %         figure(angi+1)
% %         plot_forcechain(pairoc,normalise_para,centroids);
% %         xlim([-0.025824946429269,0.125840814229269]);
% %         ylim([-0.000931894074707,0.119068105925293]);
% %         axis equal;
% %         hold on;
%         
%     else
%         continue
%     end
% 
%       
% else
%     break
% end
% end
% chalen_ave(i,1)=mean(chalen);%average chain length
% chaquant(i,1)=count_chain;%total quantity of chain

%% partition edge base on the angle between the direction of potential drop
% % e_angle belong to [0 pi] between edge and x axis
% [edgeangle_xaxis,edgeangle_yaxis,e_angle]=edgeangle_xoryaxis(pairid_force_dis,centroids);

     %% legend
%     mylgd{1,i}=char([name_split(5,1)]);
%       mylgd{1,i}=([char(name_split(3,1))]);
%       ' ;CN=' num2str(CN(i,1))
%       le=legend(mylgd{1,i},'fontsize',15);
     
     %% Network analysis gap factor
%      louvain_repeats=30;
%      resolution=0.1:0.1:3;
%      rho=1;
%      [g,com,Q,Qmax,rc,pairid_force_dis,centroids1,Stotal1,n,Fij]=Network_analysis(num1,num2,louvain_repeats,resolution,rho);
%      gap_factor(i,1:length(resolution))=g;
% %      plot(resolution,gap_factor(i,1:length(resolution)),'Color',color{1,file(1,i)});
% %      hold on;
%      
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
     
    %% find cycles
%     figure(1)
%     G=graph(pairid_force_dis(:,1),pairid_force_dis(:,2));
%     hg=plot(G,'XData',centroids(:,1),'YData',centroids(:,2));
% 
%     obj_cycle=spatialgraph2D(G,hg.XData,hg.YData);
%     [~,cyclecluster]=obj_cycle.polyshape;%cyclecluster contains the list of nodes id in the cycle
%     cyclecluster=cyclecluster';
%     size_each_cyc=[];
%     parfor cyci=1:length(cyclecluster)
%         size_each_cyc(cyci,1)=length(cell2mat(cyclecluster(cyci,1)));%size of each cycle
% %         if size_each_cyc(cyci,1)>3
% %             plot
%     end
%     mean_cycnodes(i,1)=mean(size_each_cyc(:,1));
%     cycles_total(i,1)=length(size_each_cyc);
%     cycles_largerthan3(i,1)=length(find(size_each_cyc>3));
%     cycles_largerthan4(i,1)=length(find(size_each_cyc>4));
%     cycles_largerthan5(i,1)=length(find(size_each_cyc>5));
%     cycles_largerthan6(i,1)=length(find(size_each_cyc>6));
    
     %% electro conductivity (2 to export)
% %--------------------- data from excel -------------    
    % % terminal and ground
%     id_wall=id_wallcontact-1;
%     posi_elev=id_wall(wfny<0);
%     nega_elev=id_wall(wfny>0);
    % % posi left; nega right data from excel
%     posi_eleh=id_wall(wfnx>0);
%     nega_eleh=id_wall(wfnx<0);

% %--------------------- data from txt -------------
%      posi_elev=id_wallcontact(wfny>0,1)-1;
%      nega_elev=id_wallcontact(wfny<0,1)-1;
%      posi_eleh=id_wallcontact(wfnx<0,1)-1;
%      nega_eleh=id_wallcontact(wfnx>0,1)-1;
     
    % % kirchhof current law test
%      [Resistancev(ii,1),Networkv,Icv,Uvv,resis_Networkv,pairid_Rij_currv] = TDNW_function_directcurrent_try(pairid_force_dis,posi_elev,nega_elev,centroids,2);
%      [Resistanceh(ii,1),Networkh,Ich,Uvh,resis_Networkh,pairid_Rij_currh] = TDNW_function_directcurrent_try(pairid_force_dis,posi_eleh,nega_eleh,centroids,2);
     % % make the all force to 1
%      unweighted_forcechain=pairid_force_dis;
%      unweighted_forcechain(:,3)=ones(length(pairid_force_dis),1);
%      [Resistance_unwv(ii,1),Networkvunw,Icvunw,Uvvunw,resis_Networkvunw,pairid_Rij_currvunw] = TDNW_function_directcurrent_try(unweighted_forcechain,posi_elev,nega_elev,centroids,2);
%      [Resistance_unwh(ii,1),Networkhunw,Ichunw,Uvhunw,resis_Networkhunw,pairid_Rij_currhunw] = TDNW_function_directcurrent_try(unweighted_forcechain,posi_eleh,nega_eleh,centroids,2);

     % % test U-I curve
%      input_I=1:100:100000; %A
%      for Ii=1:size(input_I',1)
%          [ResistancevUI(Ii,1),Networkv,Icv,Uvv,resis_Networkv,pairid_Rij_currv] = TDNW_function_directcurrent_try(pairid_force_dis,posi_elev,nega_elev,centroids,input_I(1,Ii));
%          U_terminal(Ii,1)=ResistancevUI(Ii,1)*input_I(1,Ii);
%      end
     
     % % gradually remove force contact below certain threshole
%      max_force=max(pairid_force_dis(:,3));
%      min_force=min(pairid_force_dis(:,3));
%      force_inter=(max_force-min_force)/100;
%      force_thre=min_force:force_inter:max_force;
%     for f_threi=1
%       pairid_force_dis_adhoc=pairid_force_dis;  
%         pairid_force_dis_adhoc(pairid_force_dis_adhoc(:,3)<force_thre(1,f_threi),:)=[];
%      [Resistance(f_threi,1),Network,Ic,Uv,resis_Network] = TDNW_function(pairid_force_dis_adhoc,posi_ele,nega_ele,centroids);
%      figure(f_threi);
%      force_chain_network=plot_forcechain(pairid_force_dis_adhoc,normalise_para,centroids);
%      hold on;
%      [mfc_pair,fc_pair,mfc_path,fc_path,mfc_aft,fc_aft,mfc_avef,fc_avef,dist]=find_main_forcechain(centroids,pairid_force_dis_adhoc,posi_ele,nega_ele,Fij,normalise_para);
%      hold on;
%      plot(centroids(posi_ele(:,1),1),centroids(posi_ele(:,1),2),'bo');
%      hold on;
%      plot(centroids(nega_ele(:,1),1),centroids(nega_ele(:,1),2),'b*');
%      
%      title(['delete fc<' num2str(force_thre(1,f_threi)) '; Resistance=' num2str(Resistance(f_threi,1))]);
%      name=['figure_' num2str(f_threi) '_threshold_' num2str(force_thre(1,f_threi)) '_Resistance_' num2str(Resistance(f_threi,1)) '.png'];
%      saveas(figure(f_threi),name);
% %      close figure
%     end

   %% plot current network
% % pairid_Rij_currv=cell2mat(pairid_Rij_currv_cell(8,1));
% figure(3)
% cus_colormapv=customise_colormap(pairid_Rij_currv(:,4),'jet');
% [current_networkv,maxIv,minIv]=plot_currentnetwork(pairid_Rij_currv,centroids,cus_colormapv);
% set(gca, 'clim', [minIv maxIv]);
% colormap jet
% axis equal;
% axis off;
% title('Current Network (vertical)');  
%     Color_Bar=colorbar;
%     Color_Bar.Label.String='Current(Amp)';
%     
% figure(4)
% cus_colormaph=customise_colormap(pairid_Rij_currh(:,4),'jet');
% [current_networkh,maxIh,minIh]=plot_currentnetwork(pairid_Rij_currh,centroids,cus_colormaph);
% set(gca, 'clim', [minIh maxIh]);
% colormap jet
% axis equal;
% axis off;
% title('Current Network (horizontal)');  
%     Color_Bar=colorbar;
%     Color_Bar.Label.String='Current(Amp)';

%% coarse-graining method
% 1. gaussian weighting

% overall corse graining------------------------------------------
cg_scale=1;
[cgcurrh,cgphi6l,cgphil] = coarse_graining(centroids,pairid_force_dis,radiucs,cg_scale,phi6l,phil,pairid_Rij_currh);

% coarse graining for sample points-----------------------------------
rsp=[4903;4335;3977;583;2020;266;480;4588;2204];%random_sample_points
% cg_s=(0:0.002:0.04)'./(mean(radius));% coarse graining scale
cg_s=(10*mean(radius))';
cgphil_rsp=[];
cgcurrh_rsp=[];
cgforcexx_rsp=[];
cgforceyy_rsp=[];
% for cg_si=1:size(cg_s,1)
for cg_si=1
    [cgcurrh,cgphi6l,cgphil,cgforcexx,cgforceyy] = coarse_graining(centroids,pairid_force_dis,radius,cg_s(cg_si,1),phi6l,phil,pairid_Rij_currh);
    cgphil_rsp(cg_si,:)=cgphil(rsp,1)';
    cgcurrh_rsp(cg_si,:)=cgcurrh(rsp,1)';
    cgforcexx_rsp(cg_si,:)=cgforcexx(rsp,1)';
    cgforceyy_rsp(cg_si,:)=cgforceyy(rsp,1)';
end
% clear cgphil_rsp cgcurrh_rsp cgforcexx_rsp cgforceyy_rsp;

% -----------------------determine coarse-graining scale-------------
figure(1)
matlab_marker_list=load('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\2d_packing\data_analysis_code\functions\matlab_marker_list');
matlab_marker_list=struct2cell(matlab_marker_list);
matlab_marker_list=matlab_marker_list{1,1};
cus_colormap_rsp=rand(size(rsp,1),3);
cg_s_matrix=repmat(cg_s,1,size(rsp,1));
for ploti=1:size(rsp,1)
plot(cg_s_matrix(:,ploti),cgphil_rsp(:,ploti),matlab_marker_list{ploti,1},'Color',cus_colormap_rsp(ploti,:));
hold on
end
leg_name_rsp=num2str(rsp);
legend(leg_name_rsp,'FontSize',12);
xlabel('Corse-graining scale(r)','FontSize',18);
ylabel('Corse-grained crystallization index','FontSize',18);

figure(2)
for ploti=1:size(rsp,1)
plot(cg_s_matrix(:,ploti),cgcurrh_rsp(:,ploti)./2,matlab_marker_list{ploti,1},'Color',cus_colormap_rsp(ploti,:));
hold on
end
legend(leg_name_rsp,'FontSize',12);
xlabel('Corse-graining scale(r)','FontSize',18);
ylabel('Corse-grained and normalised current','FontSize',18);

figure(3)
for ploti=1:size(rsp,1)
plot(cg_s_matrix(:,ploti),cgforcexx_rsp(:,ploti)./2,matlab_marker_list{ploti,1},'Color',cus_colormap_rsp(ploti,:));
hold on
end
legend(leg_name_rsp,'FontSize',12);
xlabel('Corse-graining scale(r)','FontSize',18);
ylabel('Corse-grained force F_x_x','FontSize',18);

figure(4)
for ploti=1:size(rsp,1)
plot(cg_s_matrix(:,ploti),cgforceyy_rsp(:,ploti)./2,matlab_marker_list{ploti,1},'Color',cus_colormap_rsp(ploti,:));
hold on
end
legend(leg_name_rsp,'FontSize',12);
xlabel('Corse-graining scale(r)','FontSize',18);
ylabel('Corse-grained force F_y_y','FontSize',18);

% ---------------------Plot coarse-grained contour--------------------------
figure(5)
scatter(centroids(:,1),centroids(:,2),[],cgcurrh./2,'filled');
axis equal;
axis off;
title({'r=1*radius;horizontal';['model' num2str(ii)]});
    Color_Bar=colorbar;
    Color_Bar.Label.String='Coarse-grained and Normalised current (I_c/I_t)';
    Color_Bar.Label.FontSize=15;
    colormap parula(10)

figure(6)
scatter(cgphil,cgcurrh./2);
xlabel('Coarse-grained crystallization index');
ylabel('Coarse-grained and Normalised current (I_c/I_t)');
title(['model' num2str(ii)]);
    Color_Bar2=colorbar;
    Color_Bar2.Label.String='Coarse-grained and Normalised current (I_c/I_t)';
    Color_Bar2.Label.FontSize=15;
    colormap parula(10)

figure(7)
scatter(centroids(:,1),centroids(:,2),[],cgphil,'filled');
axis equal;
axis off;
title({'Coarse-grained Crystallization index;r=10radius';['model' num2str(ii)]});
    Color_Bar3=colorbar;
    Color_Bar3.Label.String='Coarse-grained Crystallization index';
    Color_Bar3.Label.FontSize=15;
    colormap gray(10)

    figure(8)
scatter(centroids(:,1),centroids(:,2),[],cgforcexx,'filled');
axis equal;
axis off;
title({'Coarse-grained F_x_x;r=10radius';['model' num2str(ii)]});
    Color_Bar3=colorbar;
    Color_Bar3.Label.String='Coarse-grained F_x_x';
    Color_Bar3.Label.FontSize=15;
    colormap gray(10)
    
        figure(9)
scatter(centroids(:,1),centroids(:,2),[],cgforceyy,'filled');
axis equal;
axis off;
title({'Coarse-grained F_y_y;r=10radius';['model' num2str(ii)]});
    Color_Bar3=colorbar;
    Color_Bar3.Label.String='Coarse-grained F_y_y';
    Color_Bar3.Label.FontSize=15;
    colormap gray(10)

    figure(10)
scatter(cgforcexx,cgcurrh./2);
xlabel('Coarse-grained F_x_x');
ylabel('Coarse-grained and Normalised current (I_c/I_t)');
title(['model' num2str(ii)]);
    Color_Bar2=colorbar;
    Color_Bar2.Label.String='Coarse-grained and Normalised current (I_c/I_t)';
    Color_Bar2.Label.FontSize=15;
    colormap parula(10)
    
        figure(11)
scatter(cgforceyy,cgcurrh./2);
xlabel('Coarse-grained F_y_y');
ylabel('Coarse-grained and Normalised current (I_c/I_t)');
title(['model' num2str(ii)]);
    Color_Bar2=colorbar;
    Color_Bar2.Label.String='Coarse-grained and Normalised current (I_c/I_t)';
    Color_Bar2.Label.FontSize=15;
    colormap parula(10)

% ------------------- scatter all points
scatter(centroids(:,1),centroids(:,2));
axis equal;
axis off;

%----------------- mark all points
for mari=1:length(centroids)
text(centroids(mari,1),centroids(mari,2),num2str(mari));
end

%-------------------- mark sample points
rsp=[4903;4335;3977;583;2020;266;480;4588;2204];%random_sample_points
for rspi=1:length(rsp)
scatter(centroids(rsp(rspi,1),1),centroids(rsp(rspi,1),2),'MarkerFaceColor',[1 0 0])
text(centroids(rsp(rspi,1),1),centroids(rsp(rspi,1),2),num2str(rsp(rspi,1)));
hold on
end
    %% Potential contour
%     figure(4)
%     [c, h]=TriScatteredContour(centroids(:,1),centroids(:,2),Uvv(1:size(centroids,1),1),100);%use the last parameter to determine the mesh size
    
    %% save figure
%     cd('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\2d_packing\electro conductivity\results\6Jul2020\all cutoff distance');
% %     title([char(name_split(4,1)) '; PCP' num2str(pcp)]);
%     title([d(file(1,i)).name]);
%     % you created a figure and it is "current". % the following, you could have guessed
% %     set(gcf,'color','none');
% %     set(gca,'color','none'); 
% 
%     saveas(gcf,[d(file(1,i)).name '.png']);
%     close all;
    
end%i
end%ii,mulitcase for error test
end%iii,fri
end%iv,GR

%%
% % GR_fri_meanpcp(iiiiiii,1)=mean(pcp);%iii:GR;iiii:fri
% % GR_fri_stdpcp(iiiiiii,1)=std(pcp);
%     cell_pcp{iiiiiii,1}=pcp;
%     cell_pcp_m{iiiiiii,1}=pcp_m;
% % filenumcell{iiiiii,iiiii}=filenum;
% % volfra_small(iiiiii,iiiii)=mean(volfra_s);
% % GR_fri_mean_Packingfra(iiiiii,iiiii)=mean(packingfra);
% % GR_fri_stdPackingfra(iiiiii,iiiii)=std(packingfra);
% % 
% % GR_fri_mean_hydro_P(iiiiiii,1)=mean(hydro_P);
% % GR_fri_std_hydro_P(iiiiiii,1)=std(hydro_P);
%     cell_hydro_P{iiiiiii,1}=hydro_P;
%     cell_Fa1_Fa2_A{iiiiiii,1}=Fa1_Fa2_A;
%     cell_Str1_Str2_A{iiiiiii,1}=Str1_Str2_A;
% 
% % 
% % GR_fri_mean_f_cover(iiiiii,iiiii)=mean(f_cover);
% % GR_fri_std_f_cover(iiiiii,iiiii)=std(f_cover);
% % f_cover_cell{iiiiii,iiiii}=f_cover;
% %
% % GR_fri_mean_d(iiiiii,iiiii)=mean(mean_d);
% % GR_fri_std_d(iiiiii,iiiii)=std(mean_d);
% % mean_Rcell{iiiiii,iiiii}=mean_R;
% %
% % grcell_mono5000=grtmp;
% % melcell_mono5000=meltmp;
% % 
% % GR_fri_Resistancev(iiiiiii,1)=mean(Resistancev);
% % GR_fri_std_Resistancev(iiiiiii,1)=std(Resistancev);
% % GR_fri_Resistanceh(iiiiiii,1)=mean(Resistanceh);
% % GR_fri_std_Resistanceh(iiiiiii,1)=std(Resistanceh);
%     cell_Resistancev{iiiiiii,1}=Resistancev;
%     cell_Resistanceh{iiiiiii,1}=Resistanceh;
%     cell_Resistance_unwv{iiiiiii,1}=Resistance_unwv;
%     cell_Resistance_unwh{iiiiiii,1}=Resistance_unwh;
    

% % GR_fri_grainNo(iiiiii,iiiii)=mean(grainNo);
% % GR_fri_std_grainNo(iiiiii,iiiii)=std(grainNo);
% % GR_fri_mean_grainsize(iiiiii,iiiii)=mean(mean_grainsize);
% % GR_fri_std_mean_grainsize(iiiiii,iiiii)=std(mean_grainsize);
% % grainNo_cell{iiiiii,iiiii}=grainNo;
% % grainsize_cell{iiiiii,iiiii}=mean_grainsize;
% 
%     cell_mean_diameter{iiiiiii,1}=mean_d;
%     cell_mean_contactarea{iiiiiii,1}=mean_C_a;
%     cell_mean_FaT_angleF1x{iiiiiii,1}=FaT_angleF1x;
%     cell_mean_StrT_angleF1x{iiiiiii,1}=StrT_angleF1x;
%     cell_Z_m{iiiiiii,1}=Z_m;

end%v,Ns
end%vi,Sratio
end%vii,Mstress
end%viii,gravity direction

% le=legend(mylgd,'fontsize',18);
% xlabel('effective CN(CN<2 Removed)','fontsize',18);
% ylabel('Resistance','fontsize',18);
% zlabel('Probability Density','fontsize',18);
% title('horizontal test','fontsize',18);
% xlim([0 50]);

%% save filenumber to mat data file
% cd('D:\D_d4');
% datafilename=['Sratio1_Ns01_N5000_pcp_Resis_HydroPe-5_gra45deg.mat'];
% %_gra22-5deg
% % % stdRadius{1,1}=GR_fri_std_R;
% % save(datafilename,'GR_fri_mean_hydro_P');
% % save(datafilename,'GR_fri_std_hydro_P','-append');
% save(datafilename,'cell_hydro_P');
% 
% % save(datafilename,'GR_fri_meanpcp');
% % save(datafilename,'GR_fri_stdpcp','-append');
% save(datafilename,'cell_pcp','-append');
% save(datafilename,'cell_pcp_m','-append');
% 
% % 
% % save(datafilename,'grcell_mono5000','-append');
% % save(datafilename,'melcell_mono5000','-append');
% 
% % save(datafilename,'GR_fri_Resistanceh','-append');
% % save(datafilename,'GR_fri_std_Resistanceh','-append');
% save(datafilename,'cell_Resistanceh','-append');
% save(datafilename,'cell_Resistance_unwh','-append');
% 
% % save(datafilename,'GR_fri_Resistancev','-append');
% % save(datafilename,'GR_fri_std_Resistancev','-append');
% save(datafilename,'cell_Resistancev','-append');
% save(datafilename,'cell_Resistance_unwv','-append');
% 
% save(datafilename,'cell_mean_diameter','-append');
%  
% save(datafilename,'cell_mean_contactarea','-append');
% 
% save(datafilename,'cell_mean_FaT_angleF1x','-append');
% 
% save(datafilename,'cell_mean_StrT_angleF1x','-append');
% 
% save(datafilename,'cell_Z_m','-append');
% 
% save(datafilename,'cell_Fa1_Fa2_A','-append');
% 
% save(datafilename,'cell_Str1_Str2_A','-append');
% 
% save(datafilename,'cell_CN_me','-append');


% % ,'-append'