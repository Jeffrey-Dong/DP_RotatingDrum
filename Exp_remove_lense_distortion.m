clc;
clear;
% images = imread('D:\EXP RoDrtest\Exp\sharpened\C0156\C0156_00001.jpg');
images = imageDatastore('D:\EXP RoDrtest\Exp\sharpened\C0156');
%%
[imagePoints,boardSize] = detectCheckerboardPoints(images.Files);
%%
squareSize = 2;
worldPoints = generateCheckerboardPoints(boardSize,squareSize);
%%
I = readimage(images,1);  
imageSize = [size(I,1),size(I,2)];
cameraParams = estimateCameraParameters(imagePoints,worldPoints, ...
                                  'ImageSize',imageSize);

%%
J1 = undistortImage(I,cameraParams);

%%
figure; imshowpair(I,J1,'montage');
title('Original Image (left) vs. Corrected Image (right)');