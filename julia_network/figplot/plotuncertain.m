IAcasedata=readtable('IAweekly_cases_ensemble1.csv');
IAbetadata=readtable('IAdaily_beta_ensemble1.csv');
IAbasecsae=table2array(IAcasedata(:,2));
IApredcasse=table2array(IAcasedata(:,4));
IAcaseup=table2array(IAcasedata(:,7));
IAcasedown=table2array(IAcasedata(:,6));
IAbeta=table2array(IAbetadata(:,3));
IAbetaup=table2array(IAbetadata(:,6));
IAbetadown=table2array(IAbetadata(:,5));
ILcasedata=readtable('ILweekly_cases_ensemble1.csv');
ILbetadata=readtable('ILdaily_beta_ensemble1.csv');
ILbasecsae=table2array(ILcasedata(:,2));
ILpredcasse=table2array(ILcasedata(:,4));
ILcaseup=table2array(ILcasedata(:,7));
ILcasedown=table2array(ILcasedata(:,6));
ILbeta=table2array(ILbetadata(:,3));
ILbetaup=table2array(ILbetadata(:,6));
ILbetadown=table2array(ILbetadata(:,5));
INcasedata=readtable('INweekly_cases_ensemble2.csv');
INbetadata=readtable('INdaily_beta_ensemble2.csv');
INbasecsae=table2array(INcasedata(:,2));
INpredcasse=table2array(INcasedata(:,4));
INcaseup=table2array(INcasedata(:,7));
INcasedown=table2array(INcasedata(:,6));
INbeta=table2array(INbetadata(:,3));
INbetaup=table2array(INbetadata(:,6));
INbetadown=table2array(INbetadata(:,5));







week1 = 1:312;
startYear = 2013;
startWeek = 40;
startDate = datetime(startYear, 1, 1) + calweeks(startWeek - 1);
date1 = startDate + calweeks(week1 - 1);
x=1:312;
figure(1);
subplot(2,3,4)
hold on;
% 绘制均值线
plot(date1, IAbeta(1:7:end), '-', 'LineWidth', 2,'Color','#6badd7');
% 绘制置信区间阴影
fill([date1'; flipud(date1')], [IAbetadown(1:7:end); flipud(IAbetaup(1:7:end))], [110/255, 169/255, 207/255], ...
    'FaceAlpha', 0.5, 'EdgeColor', 'none');
xtickformat("yyyy")
xlabel("year", 'FontSize', 12)
ylabel("Transmission rate", 'FontSize', 12)
title("(d)", 'FontSize', 22)
% 添加图例和标签
grid on;
hold off;
subplot(2,3,5)
hold on;
% 绘制均值线
plot(date1, ILbeta(1:7:end), '-', 'LineWidth', 2,'Color','#6badd7');
% 绘制置信区间阴影
fill([date1'; flipud(date1')], [ILbetadown(1:7:end); flipud(ILbetaup(1:7:end))],[110/255, 169/255, 207/255], ...
    'FaceAlpha', 0.5, 'EdgeColor', 'none');
xtickformat("yyyy")
xlabel("year", 'FontSize', 12)
ylabel("Transmission rate", 'FontSize', 12)
title("(e)", 'FontSize', 22)
% 添加图例和标签
grid on;
hold off;
subplot(2,3,6)
hold on;
% 绘制均值线
plot(date1, INbeta(1:7:end), '-', 'LineWidth', 2,'Color','#6badd7');
% 绘制置信区间阴影
fill([date1'; flipud(date1')], [INbetadown(1:7:end); flipud(INbetaup(1:7:end))],[110/255, 169/255, 207/255], ...
    'FaceAlpha', 0.5, 'EdgeColor', 'none');
xtickformat("yyyy")
xlabel("year", 'FontSize', 12)
ylabel("Transmission rate", 'FontSize', 12)
title("(f)", 'FontSize', 22)
% 添加图例和标签
grid on;
hold off;




x1=date1;

subplot(2,3,1)
plot(x1, IAbasecsae, '*', 'LineWidth', 2,'Color','#f3683c');
hold on;

% 绘制均值线
plot(x1, IApredcasse, '-', 'LineWidth', 2,'Color','#206fb6');

% 绘制置信区间阴影
fill([x1'; flipud(x1')], [max(IAcasedown,0); flipud(IAcaseup)], [110/255, 169/255, 207/255], ...
    'FaceAlpha', 0.6, 'EdgeColor', 'none');
xtickformat("yyyy")
xlabel("year", 'FontSize', 12)
ylabel("Influenza cases", 'FontSize', 12)
legend("Real data","Fitting curve", 'FontSize', 10)
title("(a)", 'FontSize', 22)
% 添加图例和标签
grid on;
hold off;
subplot(2,3,2)
plot(x1, ILbasecsae, '*', 'LineWidth', 2,'Color','#f3683c');
hold on;

% 绘制均值线
plot(x1, ILpredcasse, '-', 'LineWidth', 2,'Color','#206fb6');

% 绘制置信区间阴影
fill([x1'; flipud(x1')], [max(ILcasedown,0); flipud(ILcaseup)], [110/255, 169/255, 207/255], ...
    'FaceAlpha', 0.6, 'EdgeColor', 'none');
xtickformat("yyyy")
xlabel("year", 'FontSize', 12)
ylabel("Influenza cases", 'FontSize', 12)
legend("Real data","Fitting curve", 'FontSize', 10)
title("(b)", 'FontSize', 22)
% 添加图例和标签
grid on;
hold off;
subplot(2,3,3)
plot(x1, INbasecsae, '*', 'LineWidth', 2,'Color','#f3683c');
hold on;

% 绘制均值线
plot(x1, INpredcasse, '-', 'LineWidth', 2,'Color','#206fb6');

% 绘制置信区间阴影
fill([x1'; flipud(x1')], [max(INcasedown,0); flipud(INcaseup)],[110/255, 169/255, 207/255], ...
    'FaceAlpha', 0.6, 'EdgeColor', 'none');
xtickformat("yyyy")
xlabel("year", 'FontSize', 12)
ylabel("Influenza cases", 'FontSize', 12)
legend("Real data","Fitting curve", 'FontSize', 10)
title("(c)", 'FontSize', 22)
% 添加图例和标签
grid on;
hold off;