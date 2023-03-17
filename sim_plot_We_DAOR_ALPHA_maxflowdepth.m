color=colormap(turbo(8));
%% figure 1 plot WE alpha
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
    errb=errorbar(num_alpha(j,i*3),-num_alpha(j,i*3-2),num_alpha(j,i*3-1),marker);
%     errb.Marker = marker;
    errb.MarkerSize = 8;
    errb.Color = color(j+1,:);
    errb.CapSize = 8;
    errb.LineWidth =1.5;
    end
end
%% figure 1 settings
xlim([10 50000]);
ylim([2 7]);
set(gca, 'XScale', 'log');
set(gca,'XTick',[0.1,1,10,100]);
legend(['\omega = 15 rpm'], ...
    ['\omega = 20 rpm'], ...
    ['\omega = 25 rpm'],'Location','southeast');
xlabel("\itWe");
ylabel('$$\alpha$$','Interpreter','latex');
set(gca,"FontSize",18)
set(gca,"FontName","times new roman")
set(gca,"Ytick",[2,3,4,5,6,7])
set(gca,"Xtick",[10,100,1000,10000])
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
for i=1:3
    if i==1
        marker="ks";
    elseif i==2
        marker="ko";
    elseif i==3
        marker="k^";
    end
    for j=1:6
    errb=errorbar(num_flowdepth(j,i*3),num_flowdepth(j,i*3-2),num_flowdepth(j,i*3-1),marker);
%     errb.Marker = marker;
    errb.MarkerSize = 8;
    errb.Color = color(j+1,:);
    errb.CapSize = 8;
    errb.LineWidth =1.5;
    end

end
%% figure 3 settings
xlim([10 50000]);
ylim([11 18]);
set(gca, 'XScale', 'log');
set(gca,'XTick',[10,100,1000,10000]);
legend(['\omega = 15 rpm'], ...
    ['\omega = 20 rpm'], ...
    ['\omega = 25 rpm'],'Location','northeast');
xlabel("\itWe");
ylabel("{\ith}_{max}/\itd_p");
set(gca,"FontSize",18)
set(gca,"FontName","times new roman");
set(gca,"LineWidth",2);
% set(gca,"Ytick",[2,3,4,5,6,7])
box on;

%% figure 4 plot We DAOR
[num_DAOR,~,~]=xlsread(['C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\parametric_study.xlsx'],11);
% [num_flowdepth,~,~]=xlsread(['C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\parametric_study.xlsx'],12);
% [num_Jarray,~,~]=xlsread(['C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\parametric_study.xlsx'],14);
figure(4);hold on;
plot(nan,nan,'ks','MarkerSize',8,LineWidth=1.5);
plot(nan,nan,'ko','MarkerSize',8,LineWidth=1.5);
plot(nan,nan,'k^','MarkerSize',8,LineWidth=1.5);
for i=1:3
    if i==1
        marker="ks";
    elseif i==2
        marker="ko";
    elseif i==3
        marker="k^";
    end
    for j=1:6
    errb=errorbar(1./num_DAOR(j,i*3),num_DAOR(j,i*3-2),num_DAOR(j,i*3-1),marker,'LineWidth',1);
    errb.MarkerSize = 10;
    errb.Color = color(j+1,:);
    errb.CapSize = 10;
    errb.LineWidth =1.2;
    end

end

% plot(num_Jarray(33:53,2),num_Jarray(33:53,3),'k*'); % Jarray results
% We(modeli,column)=2460*dp/2*(omega/180*pi*0.1)^2/surften;
for i=1:3
    if i==1
        marker="ks";
    elseif i==2
        marker="ko";
    elseif i==3
        marker="k^";
    end
    
    index=[1;4];
    for j=1:2
    errb=errorbar(1./num_DAOR(index(j,1),i*3),num_DAOR(i+17,index(j,1)),num_DAOR(i+17,index(j,1)+1),marker,'LineWidth',1,...
    'MarkerEdgeColor',color(j*j+j,:),'MarkerFaceColor',color(j*j+j,:));
    errb.MarkerSize = 8;
    errb.Color = color(j*j+j,:);
    errb.CapSize = 10;
    errb.LineWidth =1.5;
    end

end

%% figure 4 settings
xlim([-0.003 0.032]);
ylim([35 60]);
set(gca, 'XScale', 'linear');
% set(gca,'XTick',[0.1,1,10,100]);
% legend('\omega = 15 rpm', ...
%     '\omega = 20 rpm', ...
%     '\omega = 25 rpm', ...
%     'Location','southeast');
legend('', ...
    '', ...
    '', ...
    '\gamma=0', ...
    '\gamma=0.007\itN/m', ...
    '\gamma=0.018\itN/m', ...
    '\gamma=0.037\itN/m', ...
    '\gamma=0.073\itN/m', ...
    '\gamma=0.146\itN/m', ...
    'Location','southeast','Orientation','vertical','Numcolumns',2);
% xlabel("1/({\itBo}^{-1}+\beta{\itFr}^2)");
xlabel("1/{\itWe}");
ylabel("\theta_{dyna, \gamma, \omega} (\circ)");
set(gca,"FontSize",18)
set(gca,"FontName","times new roman")
set(gca,"Ytick",[35,40,45,50,55,60])
set(gca,'LineWidth',2)
box on;