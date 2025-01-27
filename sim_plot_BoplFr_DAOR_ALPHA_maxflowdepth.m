clc;
clear;
color=colormap(turbo(8));
%% calculate dimenionless number
rho=2460;%kg/m3
omega=[15;20;25].*6./180*pi;%rpm
D=0.2;%m rotation drum diameter
d=0.002;
gamma=[0;0.0073;0.0183;0.0365;0.073;0.146];%N/m surface tension
g=9.81;%gravity; m/s2
for omei=1:3
    % BOplFR(:,omei)=(gamma+0.1*rho.*omega(omei,1).^2.*D.*d^2)./(rho*g*d^2);
    % BOplFR_exp(:,omei)=(gamma*0.5+0.1*rho.*omega(omei,1).^2.*D.*d^2)./(rho*g*d^2);
    BO(:,omei)=gamma./(rho*9.81*(d)^2);
    BO_exp(:,omei)=BO(:,omei).*1;
    FR(:,omei)=sqrt(omega(omei,1)^2*D/(2*g));
    beta=8;
    BOplFR(:,omei)=sqrt(BO(:,omei))+beta*FR(:,omei);
    BOplFR_exp(:,omei)=sqrt(BO_exp(:,omei))+beta*FR(:,omei);
end
%% figure 1 plot WE alpha
[num_DAOR,~,~]=xlsread(['C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\parametric_study.xlsx'],13);
[num_alpha,~,~]=xlsread(['C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\parametric_study.xlsx'],10);
figure(1);hold on;
plot(nan,nan,'ks','MarkerSize',8,LineWidth=1.5);
plot(nan,nan,'ko','MarkerSize',8,LineWidth=1.5);
plot(nan,nan,'k^','MarkerSize',8,LineWidth=1.5);
for i=1:3
    if i==1
        marker="s";
    elseif i==2
        marker="o";
    elseif i==3
        marker="^";
    end
    
    for j=1:6
    errb=errorbar(1./BOplFR(j,i),-num_alpha(j,i*3-2),num_alpha(j,i*3-1),marker);
%     errb.Marker = marker;
    errb.MarkerSize = 8;
    errb.Color = color(j+1,:);
    errb.CapSize = 8;
    errb.LineWidth =1.5;
    end
end
%% figure 1 settings
xlim([0.1 1.1]);
ylim([2 7]);
set(gca, 'XScale', 'linear');
% set(gca,'XTick',[0.1,1,10,100]);
legend(['\omega = 15 rpm'], ...
    ['\omega = 20 rpm'], ...
    ['\omega = 25 rpm'],'Location','southeast');
xlabel("1/\itC_E");
ylabel('$$\alpha$$','Interpreter','latex');
set(gca,"FontSize",18)
set(gca,"FontName","times new roman")
set(gca,"Ytick",[2,3,4,5,6,7])
% set(gca,"Xtick",[0.00001:0.00005:0.0003])
set(gca,"LineWidth",2)
box on;

%% figure 2 plot DAOR flowdepth
color=colormap(turbo(8));
[num_DAOR,~,~]=xlsread(['C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\parametric_study.xlsx'],11);
[num_flowdepth,~,~]=xlsread(['C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\parametric_study.xlsx'],12);
figure(2);hold on;
plot(nan,nan,'ks','MarkerSize',8,LineWidth=1.5);
plot(nan,nan,'ko','MarkerSize',8,LineWidth=1.5);
plot(nan,nan,'k^','MarkerSize',8,LineWidth=1.5);
for i=1:3
    if i==1
        marker="s";
    elseif i==2
        marker="o";
    elseif i==3
        marker="^";
    end
    for ploti=1:6
        %% plot without errorbar
%         plot(num_DAOR(ploti,i*2-1),num_flowdepth(ploti,i*2-1), ...
%             marker, ...
%             'color',color(ploti+1,:), ...
%             'LineWidth',1.5, ...
%             'MarkerSize',8);
        errb=errorbar(num_DAOR(ploti,i*3-2),num_flowdepth(ploti,i*3-2),num_DAOR(ploti,i*3-1),num_DAOR(ploti,i*3-1),num_flowdepth(ploti,i*3-1),num_flowdepth(ploti,i*3-1), ...
            marker,'MarkerEdgeColor',color(ploti+1,:));
    %     errb.Marker = marker;
        errb.MarkerSize = 8;
        errb.Color = color(ploti+1,:);
        errb.CapSize = 8;
        errb.LineWidth =1.5;


    end

end
%% figure 2 settings
xlim([38 60]);
ylim([10 20]);
% set(gca, 'XScale', 'log');
% set(gca,'XTick',[0.1,1,10,100]);
legend('\omega = 15 rpm', ...
    '\omega = 20 rpm', ...
    '\omega = 25 rpm','Location','northwest');

xlabel("DAOR");
ylabel("{\itD}_{max}");
set(gca,"FontSize",18)
set(gca,"FontName","times new roman")
set(gca,"LineWidth",2);
% set(gca,"Ytick",[2,3,4,5,6,7])
box on;

%% figure 3 plot We flowdepth
[num_DAOR,~,~]=xlsread(['C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\parametric_study.xlsx'],11);
[num_flowdepth,~,~]=xlsread(['C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\parametric_study.xlsx'],12);
figure(3);hold on;
plot(nan,nan,'ks','MarkerSize',8,LineWidth=1.5);
plot(nan,nan,'ko','MarkerSize',8,LineWidth=1.5);
plot(nan,nan,'k^','MarkerSize',8,LineWidth=1.5);
CE=[];
Dmax=[];
for i=1:3
    if i==1
        marker="ks";
        linetype="k-.";
    elseif i==2
        marker="ko";
        linetype="k--";
    elseif i==3
        marker="k^";
        linetype="k:";
    end
    plot(BOplFR(:,i).^1,num_flowdepth(:,i*3-2),linetype,"LineWidth",1.5);
    for j=1:6
    errb=errorbar(BOplFR(j,i)^1,num_flowdepth(j,i*3-2),num_flowdepth(j,i*3-1),marker);
%     errb.Marker = marker;
    errb.MarkerSize = 8;
    errb.Color = color(j+1,:);
    errb.CapSize = 8;
    errb.LineWidth =1.5;
    Dmax=[Dmax;num_flowdepth(j,i*3-2)];
    CE=[CE;1./BOplFR(j,i)];
    end
end
%% guide of eye fit
[xData, yData] = prepareCurveData( CE, Dmax );
% Set up fittype and options.
ft = fittype( 'power2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [13.2713891861167 -0.19788955298435 -0.013691674029998];

% Fit model to data.
[fitresult, ~] = fit( xData, yData, ft, opts );
fitx=0:0.01:2.5;
fity=fitresult.a.*fitx.^fitresult.b+fitresult.c;
h = plot(fitx,fity,"k--");
h.LineWidth=2;

%% figure 3 settings
figure(3)
% xlim([0 2.5]);
% ylim([11 18]);
% xlim([-5 45]);% 1/C_E
% ylim([-1 6]);% 1/C_E
xlim([1 3.5]);% C_E
ylim([11 18]);
set(gca, 'XScale', 'linear');
% set(gca,'XTick',[0,5000,10000,15000,20000,25000]);
% set(gca,"Xtick",[0:0.5:2.5])
% legend(['\omega = 15 rpm'], ...
%     ['\omega = 20 rpm'], ...
%     ['\omega = 25 rpm'],'Location','northeast');
xlabel("\itC_E");
% ylabel("({\ith}_{max,wet} - {\ith}_{max,dry}) / \itd_p");
ylabel("{\ith}_{max} / \itd_p");
set(gca,"FontSize",18)
set(gca,"FontName","times new roman");
set(gca,"LineWidth",2);

box on;

%% figure 4 plot We DAOR
color=colormap(turbo(8));
[num_DAOR,~,~]=xlsread('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\parametric_study.xlsx',11);
% [num_flowdepth,~,~]=xlsread(['C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\parametric_study.xlsx'],12);
% [num_Jarray,~,~]=xlsread(['C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\parametric_study.xlsx'],15);
figure(4);hold on;
plot(nan,nan,'ks','MarkerSize',8,LineWidth=1.5);
plot(nan,nan,'ko','MarkerSize',8,LineWidth=1.5);
plot(nan,nan,'k^','MarkerSize',8,LineWidth=1.5);
CE=[];
DAOR=[];
Bo=[];
Fr=[];
for i=1:3
    if i==1
        marker="s";
    elseif i==2
        marker="o";
    elseif i==3
        marker="^";
    end
    for j=1:6
    errb=errorbar(BOplFR(j,i)^1,num_DAOR(j,i*3-2)-22.3,num_DAOR(j,i*3-1),marker,'LineWidth',1);
    errb.MarkerSize = 10;
    errb.Color = color(j+1,:);
    errb.CapSize = 10;
    errb.LineWidth =1.2;
    DAOR=[DAOR;num_DAOR(j,i*3-2)-22.3];
    CE=[CE;BOplFR(j,i)];
    Bo=[Bo;BO(j,i)];
    Fr=[Fr;FR(:,i)'];

    end
end

% plot(num_Jarray(33:53,2),num_Jarray(33:53,3),'k*'); % Jarray results
% We(modeli,column)=2460*dp/2*(omega/180*pi*0.1)^2/surften;
for i=1:3
    if i==1
        marker="s";
    elseif i==2
        marker="o";
    elseif i==3
        marker="^";
    end
    index=[1;4];
    for j=1:2
    errb=errorbar(BOplFR_exp(index(j,1),i)^1,num_DAOR(i+17,index(j,1))-22.3,num_DAOR(i+17,index(j,1)+1),marker,'LineWidth',1,...
    'MarkerEdgeColor',color(index(j,1)+1,:),'MarkerFaceColor',color(index(j,1)+1,:));
    errb.MarkerSize = 8;
    errb.Color = color(index(j,1)+1,:);
    errb.CapSize = 10;
    errb.LineWidth =1.5;
    DAOR=[DAOR;num_DAOR(i+17,index(j,1))-22.3];
    CE=[CE;BOplFR(index(j,1),i)];
    Bo=[Bo;BO(index(j,1),i)];
    Fr=[Fr;FR(1,i)];
    end

end
%% guide of eye fit
[xData, yData] = prepareCurveData( CE, DAOR );

% Set up fittype and options.
% ft = fittype( 'power2' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
% opts.StartPoint = [46.3047222040371 -0.202865003729663 -0.0321392907141581];

ft = fittype( 'a*x+b', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.585267750979777 0.223811939491137];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
fitx=0:0.01:10;
% fity=fitresult.a.*fitx.^fitresult.b+fitresult.c;
fity=fitresult.a.*fitx+fitresult.b;
h = plot(fitx,fity,"k--");
h.LineWidth=2;


%% figure 4 settings
% xlim([0 2.5]);
% ylim([35 60]);
% xlim([-5 45]);% 1/C_E
% ylim([-1 10]);% 1/C_E
xlim([1 3.5]);% C_E
% ylim([-2 11])
% set(gca, 'XScale', 'log');
% set(gca,'XTick',[1:0.5:4.5]);%beta 20 1/
% set(gca,'XTick',[0:0.5:2.5]);%beta 10 1/
% set(gca,"Ytick",[35,40,45,50,55,60])
legend('\omega = 15 rpm', ...
    '\omega = 20 rpm', ...
    '\omega = 25 rpm', ...
    'Location','southeast');
xlabel("\itC_E");
% xlabel("1/({\itWe}");
% ylabel("DAOR (\circ)");
% ylabel("\theta_{wet} - \theta_{dry}, (\circ)");
ylabel("\theta_{dyna} - \theta_{\itd} (\circ)");
set(gca,"FontSize",18)
set(gca,"FontName","times new roman")
set(gca,'LineWidth',2)
box on;

%% figure 5 plot alpha flowdepth
figure(5)
color=colormap(turbo(8));
[num_DAOR,~,~]=xlsread(['C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\parametric_study.xlsx'],11);
[num_flowdepth,~,~]=xlsread(['C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\parametric_study.xlsx'],12);
[num_alpha,~,~]=xlsread(['C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\parametric_study.xlsx'],10);
[num_Nmean,~,~]=xlsread(['C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\parametric_study.xlsx'],14);
figure(5);hold on;
plot(nan,nan,'ks','MarkerSize',8,LineWidth=1.5);
plot(nan,nan,'ko','MarkerSize',8,LineWidth=1.5);
plot(nan,nan,'k^','MarkerSize',8,LineWidth=1.5);
for i=1:3
    if i==1
        marker="s";
    elseif i==2
        marker="o";
    elseif i==3
        marker="^";
    end
    for ploti=1:6
        %% plot without errorbar
%         plot(num_DAOR(ploti,i*2-1),num_flowdepth(ploti,i*2-1), ...
%             marker, ...
%             'color',color(ploti+1,:), ...
%             'LineWidth',1.5, ...
%             'MarkerSize',8);
        errb=errorbar(-num_alpha(ploti,i*3-2),num_flowdepth(ploti,i*3-2),num_flowdepth(ploti,i*3-1),num_flowdepth(ploti,i*3-1),num_alpha(ploti,i*3-1),num_alpha(ploti,i*3-1), ...
            marker,'MarkerEdgeColor',color(ploti+1,:));
    %     errb.Marker = marker;
        errb.MarkerSize = 8;
        errb.Color = color(ploti+1,:);
        errb.CapSize = 8;
        errb.LineWidth =1.5;


    end

end
%% figure 5 settings
xlim([2 7]);
ylim([11 18]);
% set(gca, 'XScale', 'log');
% set(gca,'XTick',[0.1,1,10,100]);
legend('\omega = 15 rpm', ...
    '\omega = 20 rpm', ...
    '\omega = 25 rpm','Location','southwest');

xlabel("\alpha");
ylabel("{\ith}_{max}/\itd_p");
set(gca,"FontSize",18)
set(gca,"FontName","times new roman")
set(gca,"LineWidth",2);
% set(gca,"Ytick",[2,3,4,5,6,7])
box on;