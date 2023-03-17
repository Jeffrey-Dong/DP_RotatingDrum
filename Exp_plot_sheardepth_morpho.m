clc;
clear;
%%
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis')
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis\function')
filedir=('e:\EXP RoDrtest\Exp\sharpened\');
dp=0.002;%particle diameter
filenum=[157:162,177:182]';
[exprec_num,exprec_txt,exprec_raw]=xlsread('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\parametric_study.xlsx',9);
omega=exprec_num(2:7,1);
watercontent=exprec_num(1,1);
%%
flowdepthmatrix=zeros(size(omega,1),size(watercontent,2));
% color=colormap(turbo(16));
color=flip(colormap(gray(16)));
for fi=2:size(exprec_raw,1)
    for fj=2:size(exprec_raw,2)
        if fj<3
            load([filedir exprec_raw{fi,fj} '\' exprec_raw{fi,fj} '_strr-coe80_strainrate_depth_2.mat']);
        else
            load([filedir exprec_raw{fi,fj} '\' exprec_raw{fi,fj} '_strr-coe20_strainrate_depth_2.mat']);
        end
        flowdepthmatrix(fi-1,fj-1)=round(maxdepth_norm,1);
        figure(1);hold on;
%         scatter(omega(fi-1,1)/60,flowdepthmatrix(fi-1,fj-1),'o','filled','MarkerFaceColor',color(fj*2,:));
        plot(omega(fi-1,1)/60,flowdepthmatrix(fi-1,fj-1),'o', ...
            'MarkerFaceColor',color(fj*3,:), ...
            MarkerEdgeColor='none');
    end
end

xlabel('Rotation speed (rpm)');
ylabel('Normalized flow region depth (D/d)');
set(gca,'FontSize',18);
set(gca,'FontName','Times New Roman');
legend('sphere', 'morphology','Location','northwest');
set(gca,'LineWidth',1);
box on
xlim([0 2.6]);