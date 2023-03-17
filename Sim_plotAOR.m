clc;
clear;

[num,txt,raw]=xlsread('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\parametric_study.xlsx',3);
%%
% DAOR=[]
for coli=1:size(num,2)
    DAORtmp=num(:,coli);
    DAOR=DAORtmp(~isnan(DAORtmp));
    plot(DAOR);
    hold on;
    clear DAOR DAORtmp
end