clc;
clear;
% % check boundary extraction quality
filenum=flip([177:182]');
for fi=1:size(filenum,1)
filedir=(['e:\EXP RoDrtest\Exp\sharpened\C0' num2str(filenum(fi)) '\']);
fig1=imread([filedir 'uppersurfaceline_expfig_2.jpg']);
fig2=imread([filedir 'uppersurface_cloudpoints_2.jpg']);
%%
figure(1);hold on;
subplot(3,2,fi);
imshow(fig1);
title(['C0' num2str(filenum(fi,1))]);

%%
figure(2);hold on;
subplot(3,2,fi);
imshow(fig2);
title(['C0' num2str(filenum(fi,1))]);

end