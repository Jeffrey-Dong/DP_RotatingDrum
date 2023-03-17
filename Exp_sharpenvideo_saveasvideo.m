clc;
clear;
%% issues: below works, but only with sharpen method 2, method 1 can create negative value of the clour, however method 2 is too slow with 1.2s/image
% load trees
% b= im2frame(X,map)
% imagesc (b.cdata)
%% read file dir
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\dataanalysis\');
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\2d_packing\data_analysis_code\functions\');
filedir=(['E:\EXP RoDrtest\Exp']);
videoinfo=dir([filedir '\*.mp4']);
[videoinfo]=sortnamebysequence(videoinfo);
[Rtxt,Ctxt]=size(videoinfo);
%%
cd(filedir);
for fi=1
%     :Rtxt
    videoname=videoinfo(fi).name;
    vr = VideoReader(['E:\EXP RoDrtest\Exp\' videoname]);
    videoNo=split(videoname,'.');
    % vr.CurrentTime=0;
    numFrames = ceil(vr.FrameRate*vr.Duration);
    vw = VideoWriter(['E:\EXP RoDrtest\Exp\' videoNo{1}],'Motion JPEG AVI');
    vw.Quality = 100;
    open(vw);
    % vidObj2= VideoWriter('xyz.avi');
    % vidObj2.FrameRate=500;
    % frame20 = read(vidObj,20);
    % % vidObj.CurrentTime;
    for fri=1
%         :10
%         :numFrames
        frame = read(vr,fri);
        frame_gray=rgb2gray(frame);
        %% sharpen image method 1

% %         % image = gpuArray(imread('E:\EXP RoDrtest\Exp\for test\C0055.MP4_20220225_094202.186.jpg'));
%         dimage = im2double(frame_gray);
%         gradient = convn(dimage,ones(3)./9,'same') - convn(dimage,ones(5)./25,'same');
%         amount = 5;
%         sharpened = dimage + amount.*gradient;

% %     %     imshow(imresize([dimage, sharpened],0.7)); %compare original and sharpened image
        %% sharpen image method 2
        sharpened = imsharpen(frame,'Radius',2,'Amount',5);
%         imshow(imresize([frame, sharpened],0.7)); %compare original and sharpened image
    %% write new video
    %     v = VideoWriter('E:\EXP RoDrtest\Exp\vidObj2');
    %     open(v);
    %     while(hasFrame(vidObj))
    %         frame = readFrame(vidObj);

        imshow(sharpened);
%         %%
%         % load mandrill
%         % figure
%         % image(X)
%         map=colormap;
%         % axis off
%         F = struct('cdata',[],'colormap',[]);
%         % for j = 1:8
%         q = 2^(9);
%         [Y,newmap] = imapprox(fig,map,q,'nodither');
%         F = im2frame(Y,newmap);
%         % end
        %%
%         frame_write=im2frame(sharpened);
%         frame_write = getframe(gca);
%         writeVideo(vw,frame_write);
        saveas(gcf,['test' '.png']);
        %     pause(1/vidObj.FrameRate/20);
    %     end
    % close all
    end
    % open(vidObj2);
    % while(hasFrame(vidObj))
    %     frame = readFrame(vidObj);
    % %     imshow(frame);
    %     vidObj2.writeVideo(frame);
    % %     pause(1/vidObj.FrameRate/20);
    % end
    close(vw);
    close all;
end

% load mandrill
% figure
% image(X)
% colormap(map)
% axis off
% F(8) = struct('cdata',[],'colormap',[]);
% for j = 1:8
% q = 2^(9-j);
% [Y,newmap] = imapprox(X,map,q,'nodither');
% F(j) = im2frame(Y,newmap);
% end