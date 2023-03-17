clc;
clear;
% close all;
%%
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis')
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis\function')
filedir=('e:\EXP RoDrtest\Exp\sharpened\');
dp=0.002;%particle diameter
Foldersdir_atominfo = (['E:\RoDrtest\RoDr_alphastart0-1\model']);
[num,txt,raw]=xlsread(['C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\parametric_study.xlsx'],1);
keypara=num(4:end,:);
omega=keypara(:,5);
d_p=keypara(:,7);
liquid_content=keypara(:,6);
%% plot cluster size vs variables
% model_var="omega";
% model_var="mu_s";
% model_var="mu_r";
model_var="gamma";
% model_var="liquid_content";
% model_var="dry_mu_r";
% model_var="all";
dbscan="on";
if model_var=="omega"
modelnum=[86:89]';% vary rotation speed
elseif model_var == "mu_s"
modelnum=[90,87,91]';% mu_s
elseif model_var == "mu_r"
modelnum=[92,87,93,94]';% mu_r
elseif model_var == "gamma"
% modelnum=[77,95,96,87,97]';% vary surface tension gamma
% modelnum=[139,124,119]';% omega90
modelnum=[119,120,134,121:123]';% omega90
% modelnum=[124,125,151,126:128,]';% omega120
% modelnum=[139,140,152,141:143]';% omega150
elseif model_var == "liquid_content"
modelnum=[77,87,99:102]';% vary liquid content
elseif model_var == "dry_mu_r"
modelnum=[77,81,82]';
elseif model_var == "all"
% modelnum=[121,126]';
end
%% marker
if modelnum(1,1)==119
    marker="s";
elseif modelnum(1,1)==124
    marker="o";
elseif modelnum(1,1)==139
    marker="^";
else
    marker=".";
end
%% custom colormap
color=colormap(turbo(9));
%%
for modeli=1:size(modelnum,1)
    %% model parameters
    omega=num(modelnum(modeli,1)+3,5);%degree/seconds
    dp=num(modelnum(modeli,1)+3,7);% meter
    mur(modeli,1)=keypara(modelnum(modeli,1),10);
    mus=keypara(modelnum(modeli,1),9);
    surften(modeli,1)=keypara(modelnum(modeli,1),11);
    liquidcontent=keypara(modelnum(modeli,1),6);
    rho=2460;%kg/m3
    WE(modeli,1)=2460*dp/2*(omega/180*pi*0.1)^2/surften(modeli,1);
    %% extract MQ infor

        matdata=load([Foldersdir_atominfo num2str(modelnum(modeli,1)) '\model' num2str(modelnum(modeli,1)) '_MQ.mat']);

        manuel_edges=[1:1e4];
        log_edge=log10(manuel_edges);
    


end

