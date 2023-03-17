clc;
clear;

addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis\github_repo');
%%
imagefilepath=('e:\EXP RoDrtest\Exp\sharpened\C0177\');
fig=imread([imagefilepath 'C0177_00001.jpg']);
figure(2)
imshow(fig);
%%
load('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis\pivlab_masks\PIVlab_mask1_c0177-c0182.mat');
x=cell2mat(maskiererx(1,1));
y=cell2mat(maskierery(1,1));

% drum_cen=[mean(x) mean(y)];

p1=[x(1) y(1)];
p2=[x(2) y(2)];
p3=[x(3) y(3)];


%% plot circle and export calibration figure for piv analysis
hold on
[CC,Radius]=CircleThru3Dots(p1,p2,p3);
xlim([0 1920]);
ylim([0 1080]);
% export_fig calibration_c0157-c0162.jpg -native
%% adjust calicoe to obtain the mask to remove the area out of the tube circle
pixeltorealdimension=0.2191/1630;
% calicoe=1.26; %calibration coefficient
% figure(1)
theta=linspace(0,2*pi,1001);
figure(2);hold on;

xcorr2=CC(1)+0.1/pixeltorealdimension*cos(theta);
ycorr2=CC(2)+0.1/pixeltorealdimension*sin(theta);
plot(xcorr2,ycorr2,'r--');

%%
% % ----------make the cell of mask for PIVlab readable format-----------
markouter=[0,0;0,1080;1920,1080;1920,0;0,0];
maskiererx{1,1}=[markouter(:,1);xcorr2'];
maskiererx{1,2}=[markouter(:,1);xcorr2'];
maskierery{1,1}=[markouter(:,2);ycorr2'];
maskierery{1,2}=[markouter(:,2);ycorr2'];
calied_pixelR=0.1/pixeltorealdimension;
% % ----------save mask and circle data
cd('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis\pivlab_masks');
maskfilename=['PIVlab_mask2_c0177-c0182.mat'];
save(maskfilename,'maskiererx');
save(maskfilename,'maskierery','-append');

centrefilename=['C0177-C0182_CR.mat'];
save(centrefilename,'CC');
save(centrefilename,'calied_pixelR','-append');

centrefilename=['calibration_coe_c0177-c0182.mat'];
save(centrefilename,'pixeltorealdimension');

%% extract particle boundary from the original picture or grayed picture and plot with shear rate
% 1. run Exp_strainrate_pivlab
% 2. load the image manuely e.g. C0083_00100
[BW,maskedImage] = Exp_segmentImage_fun_grey(C0083_00100);
xtmp=[1:1920];
ytmp=flip([1:1080]');
xtmprepmat=repmat(xtmp,1080,1);
ytmprepmat=repmat(ytmp,1,1920);
double(BW);
ztmpreshape=reshape(double(BW),[],1);
figure(4)
contourf(xtmprepmat,ytmprepmat,double(BW));
xtmpreshape=reshape(xtmprepmat,[],1);
ytmpreshape=reshape(ytmprepmat,[],1);
CCreal=CC*pixeltorealdimension;
disttoCCreal=sqrt((xtmprepmat*pixeltorealdimension-CCreal(1)).^2+(ytmprepmat*pixeltorealdimension-CCreal(2)).^2);
disttoCCrealreshape=reshape(disttoCCreal,[],1);
index_incircle=find(disttoCCrealreshape<=0.1);
index_particle=find(ztmpreshape==1);
index_particle_incircle=index_incircle(ismember(index_incircle,index_particle));
x_particle_incircle=xtmpreshape(index_particle_incircle)*pixeltorealdimension;
y_particle_incircle=ytmpreshape(index_particle_incircle)*pixeltorealdimension;
figure(5)
plot(x_particle_incircle,y_particle_incircle,'ko');
hold on;
k=boundary(x_particle_incircle,y_particle_incircle,1);
plot(x_particle_incircle(k),y_particle_incircle(k),'r-');
figure(2)
hold on
plot(x_particle_incircle(k),y_particle_incircle(k),'r-');

%%
function [CC,Radius]=CircleThru3Dots(A,B,C)
Ah=A*A';
Bh=B*B';
Ch=C*C';
CC=zeros(size(A));
G=(C(2)-B(2))*A(1)+(A(2)-C(2))*B(1)+(B(2)-A(2))*C(1);
CC(1)=((Bh-Ch)*A(2)+(Ch-Ah)*B(2)+(Ah-Bh)*C(2))/(2*G);
CC(2)=-((Bh-Ch)*A(1)+(Ch-Ah)*B(1)+(Ah-Bh)*C(1))/(2*G);
Radius=sqrt((A-CC)*(A-CC)');
theta=linspace(0,2*pi,101);
x=CC(1)+Radius*cos(theta);
y=CC(2)+Radius*sin(theta);
plot(x,y,'r-')
ABC=[A;B;C];
hold on
plot(ABC(:,1),ABC(:,2),'b.','markersize',20)
plot(CC(1),CC(2),'r.','markersize',20)
grid on
box off
axis equal
end

