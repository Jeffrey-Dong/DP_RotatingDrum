clc;
clear;
%%
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis')
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis\function')
filedir=('e:\EXP RoDrtest\Exp\sharpened\');
dp=0.002;%particle diameter
filenum=[157,152,161,153]';
[exprec_num,exprec_txt,exprec_raw]=xlsread( ...
    'C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\parametric_study.xlsx',8);
omega=exprec_num(2:5,1);
watercontent=exprec_num(1,2:3);
%%
flowdepthmatrix=zeros(size(omega,1),size(watercontent,2));
% color=colormap(turbo(16));
color=flip(colormap(hot(16)));
for fi=2:5
    for fj=2:4
        load([filedir exprec_raw{fi,fj} '\' exprec_raw{fi,fj} '_strr-coe10_strainrate_depth_2.mat']);
        flowdepthmatrix(fi-1,fj-1)=round(maxdepth_norm,1);
        figure(1);hold on;
%         scatter(omega(fi-1,1)/60,flowdepthmatrix(fi-1,fj-1),'o','filled','MarkerFaceColor',color(fj*2,:));
        if fj==2
            marker="o";
        elseif fj==3
            marker="s";
        elseif fj==4
            marker="^";
        end
        plot(omega(fi-1,1)/6,flowdepthmatrix(fi-1,fj-1),marker, ...
            'MarkerFaceColor','none', ...
            MarkerEdgeColor=color(fi*3,:));
    end
end

xlabel('Rotation speed (rpm)');
ylabel('Normalized flow region depth (D/d)');
set(gca,'FontSize',18);
set(gca,'FontName','Times New Roman');
legend('dry', ...
    '0.02', ...
    '0.05', ...
    'Location','northwest');
set(gca,'LineWidth',1);
box on
xlim([8 27]);
ylim([6 14]);