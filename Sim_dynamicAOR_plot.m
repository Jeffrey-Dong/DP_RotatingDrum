clc;
clear;
%%
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\dataanalysis\');
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\2d_packing\data_analysis_code\functions\');
[num,txt,raw]=xlsread(['C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\parametric_study.xlsx'],1);
keypara=num(4:end,:);
density_threshold=0.6;
Webernumber=num(4:end,18);
cust_grad_color=colormap(jet(6));
% sample=[49,46,50]';
% sample=[51,46,52]';
% sample=[33,59:63]';
sample=[4]';
for modeli=1:size(sample,1)
    if sample(modeli)<59
        datafile=(['E:\RoDrtest\RoDr_alphastart0-1\model' num2str(sample(modeli,1))]);
    else
        datafile=(['D:\RoDrtest\RoDr_alphastart0-1\model' num2str(sample(modeli,1))]);
    end

load([datafile '\model' num2str(sample(modeli,1)) '_DAOR_gof.mat'])


DEMtimestep=10^-7;
DEMstepinterval=50000;
omega=keypara(modeli,5);

%% plot DAOR with Webernumber
% figure(13)
% % plot(revolution(50)/360,DAOR(50),'*','color',cust_grad_color(modeli,:));
% plot(1/Webernumber(sample(modeli,1)),DAOR(50),'*','color',cust_grad_color(modeli,:));
% legrecord{1,modeli}=['\gamma=' num2str(keypara(sample(modeli,1),11))];
% hold on;

%% plot DAOR with revolution
figure(13)
plot(revolution(1:size(DAOR,1)),DAOR,'*','color',cust_grad_color(modeli,:));
% plot(1/Webernumber(sample(modeli,1)),DAOR(50),'*','color',cust_grad_color(modeli,:));
% legrecord{1,modeli}=['\gamma=' num2str(keypara(sample(modeli,1),11))];
hold on;
   
end
xlabel('1/We');
ylabel('AOR(rad)');
% legend('\omega=60^o/s,Fr=0.053','\omega=90^o/s,Fr=0.079','\omega=120^o/s,Fr=0.106','\omega=150^o/s,Fr=0.132')%model45-48
% legend('\mu_s=0.1','\mu_s=0.5','\mu_s=0.8')%model49,46,50
% legend('\mu_r=0.001','\mu_r=0.01','\mu_r=0.1')%model51,46,52
% legend('\omega=90^o/s,Fr=0.079,1/We=0.826','\omega=120^o/s,Fr=0.106,1/We=0.465','\omega=150^o/s,Fr=0.132,1/We=0.299')%model33,53,54
% legend('\mu_s=0.1','\mu_s=0.5','\mu_s=0.8')%model55,33,56
% legend('\mu_r=0.001','\mu_r=0.01','\mu_r=0.1')%model55,33,56
xlim([0 0.1]);
grid on
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