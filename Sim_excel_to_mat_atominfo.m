clc;
clear;
%%
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\dataanalysis\');
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\2d_packing\data_analysis_code\functions\');
for modeli=103
datafile=(['e:\RoDrtest\RoDr_alphastart0-1\model' num2str(modeli)]);
cd(datafile);
atominfo=dir([datafile '\*.csv']);
[atominfo]=sortnamebysequence(atominfo);
[Rtxt,Ctxt]=size(atominfo);
DAORt=[];

sample_range=[1:Rtxt]';
sample_percentage=[1:round(Rtxt/20):Rtxt]'/Rtxt;
% sample_range=[300]';
%%
for targi=1:Rtxt
    if ismember(targi/Rtxt,sample_percentage)
        disp(['model' num2str(modeli) '_' num2str(targi/Rtxt)]);
    end
% %
% cd('E:\RoDrtest\largerparticlestest');
[num,~,~]=xlsread([datafile '\' atominfo(sample_range(targi)).name],1);%[num,txt,raw]
% scatter3(num(:,18),num(:,19),num(:,20));
% scatter(num(:,18),num(:,20));
% tmpx=num(:,18);
% tmpy=num(:,19);
% tmpz=num(:,20);
centroids{targi,1}=[num(:,18) num(:,19) num(:,20)];
id{targi,1}=num(:,1);
velocity{targi,1}=[num(:,6) num(:,7) num(:,8)];
radius{targi,1}=[num(:,12)];
end
%%
   datafilename=['model' num2str(modeli) '_posi.mat'];
   save(datafilename,'centroids','-v7.3');
   clear centroids
   datafilename=['model' num2str(modeli) '_id.mat'];
   save(datafilename,'id','-v7.3');
   clear id
   datafilename=['model' num2str(modeli) '_velocity.mat'];
   save(datafilename,'velocity','-v7.3');
   clear velocity 
   datafilename=['model' num2str(modeli) '_radius.mat'];
   save(datafilename,'radius','-v7.3');
   clear radius
end
%    clear phi_value_total_model T_value_total_model;