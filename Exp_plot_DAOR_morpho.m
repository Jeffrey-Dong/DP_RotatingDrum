clc;
clear;
%%
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis')
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis\function')
filedir=('i:\EXP RoDrtest\Exp\sharpened\');
dp=0.002;%particle diameter
filenum=[157:162,177:182]';
[exprec_num,exprec_txt,exprec_raw]=xlsread('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\parametric_study.xlsx',9);
omega=exprec_num(2:7,1);
watercontent=exprec_num(1,1);
%%
DAORmatrix=zeros(size(omega,1),size(watercontent,2));
color=flip(colormap(gray(16)));
for fi=2:size(exprec_raw,1)
    for fj=2:size(exprec_raw,2)
        load([filedir exprec_raw{fi,fj} '\' exprec_raw{fi,fj} '_DAOR_2.mat']);
        DAORmatrix(fi-1,fj-1)=DAOR/pi*180;
        
        %% DAOR vs rotation speed
        figure(1);hold on;
        scatter(omega(fi-1,1)/60,DAORmatrix(fi-1,fj-1),'o','filled','MarkerFaceColor',color(fj*3,:));
        %% DAOR vs We
        figure(2);hold on;
        We=2500*0.5*dp*(0.1*omega(fi-1,1)/60)^2/(4*0.072*0.5);
        scatter(We,DAORmatrix(fi-1,fj-1),'o','filled','MarkerFaceColor',color(fj*3,:));
        %% DAOR vs l*
        figure(3);hold on;
        lstar=2500*4/3*pi*(0.5*dp)^3*(omega(fi-1,1)/180*pi)^2/(2*pi*0.072*1);
        scatter(lstar,DAORmatrix(fi-1,fj-1),'o','filled','MarkerFaceColor',color(fj*3,:));
    end
end
figure(1)
xlabel('Rotation speed (rpm)');
ylabel('Dynamic angle of repose (^\circ)');
set(gca,'FontSize',18);
set(gca,'FontName','Times New Roman');
legend('sphere', 'morphology',Location='northwest');
set(gca,'LineWidth',1);
box on
% xlim([0 160]);

figure(2)
xlabel('We');
ylabel('Dynamic angle of repose (^\circ)');
set(gca,'FontSize',18);
set(gca,'FontName','Times New Roman');
legend('sphere', 'morphology');
set(gca,'LineWidth',1);
box on
% xlim([0 160]);

figure(3)
xlabel('{\itl}^*');
ylabel('Dynamic angle of repose (^\circ)');
set(gca,'FontSize',18);
set(gca,'FontName','Times New Roman');
legend('sphere', 'morphology');
set(gca,'LineWidth',1);
box on
% xlim([0 160]);
