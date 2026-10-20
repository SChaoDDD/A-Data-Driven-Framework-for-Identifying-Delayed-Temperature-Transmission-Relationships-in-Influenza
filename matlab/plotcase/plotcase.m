clear;
clc;
close all;

%% =========================================================
% Read data
% ==========================================================

IA = readtable('IAEMA.csv');
IAo = table2array(IA(:,1)) .* 2.74;

IL = readtable('ILEMA.csv');
ILo = table2array(IL(:,1)) .* 2.74;

IN = readtable('INEMA.csv');
INo = table2array(IN(:,1)) .* 2.74;

IAfit = load('IAfit.mat');
ILfit = load('ILfit.mat');
INfit = load('INfit.mat');

%% Model reconstruction

IApred = IAfit.newcase;
ILpred = ILfit.newcase;
INpred = INfit.newcase;

% Convert to column vectors
IApred = IApred(:);
ILpred = ILpred(:);
INpred = INpred(:);

IAo = IAo(:);
ILo = ILo(:);
INo = INo(:);

%% =========================================================
% Time axis
% ==========================================================

startDate = datetime(2013,9,30);

% newcase normally contains 311 weekly intervals
nIA = min(length(IAo), length(IApred));
nIL = min(length(ILo), length(ILpred));
nIN = min(length(INo), length(INpred));

dateIA = startDate + calweeks(0:nIA-1);
dateIL = startDate + calweeks(0:nIL-1);
dateIN = startDate + calweeks(0:nIN-1);

yearTicks = datetime(2014:2019,1,1);

%% =========================================================
% Plotting parameters
% ==========================================================

fontName       = 'Times New Roman';
fontSizeTick   = 9;
fontSizeLabel  = 10;
fontSizePanel  = 13;
fontSizeLegend = 9;

axisLineWidth = 0.8;
fitLineWidth  = 1.8;

% Solid observed points
markerSize = 4.0;

% Same colors as OH / SD
obsColor = [210 85 50] / 255;       % muted dark orange-red
fitColor = [0.16 0.45 0.68];        % blue

% Figure size
figWidth  = 22.0;     % cm
figHeight = 6.5;      % cm

%% =========================================================
% Create figure
% ==========================================================

fig = figure( ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[2 2 figWidth figHeight]);

tl = tiledlayout(fig,1,3, ...
    'TileSpacing','compact', ...
    'Padding','loose');

%% =========================================================
% (a) IA
% ==========================================================

ax1 = nexttile(tl,1);
hold(ax1,'on');

hObs = plot(ax1, ...
    dateIA, IAo(1:nIA), ...
    'o', ...
    'LineStyle','none', ...
    'MarkerSize',markerSize, ...
    'MarkerEdgeColor',obsColor, ...
    'MarkerFaceColor',obsColor);

hFit = plot(ax1, ...
    dateIA, IApred(1:nIA), ...
    '-', ...
    'Color',fitColor, ...
    'LineWidth',fitLineWidth);

xlabel(ax1,'Year', ...
    'FontName',fontName, ...
    'FontSize',fontSizeLabel);

ylabel(ax1,'Estimated influenza cases', ...
    'FontName',fontName, ...
    'FontSize',fontSizeLabel);

text(ax1,0.03,0.95,'(a) IA', ...
    'Units','normalized', ...
    'FontName',fontName, ...
    'FontSize',fontSizePanel, ...
    'FontWeight','bold', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','top');
ylim([0 300])
lgd1 = legend(ax1, ...
    {'Real data','Fitting curve'}, ...
    'FontName',fontName, ...
    'FontSize',fontSizeLegend, ...
    'Box','off');

lgd1.Units = 'normalized';
lgd1.Position = [0.22 0.82 0.14 0.10];
hold(ax1,'off');

%% =========================================================
% (b) IL
% ==========================================================

ax2 = nexttile(tl,2);
hold(ax2,'on');

plot(ax2, ...
    dateIL, ILo(1:nIL), ...
    'o', ...
    'LineStyle','none', ...
    'MarkerSize',markerSize, ...
    'MarkerEdgeColor',obsColor, ...
    'MarkerFaceColor',obsColor);

plot(ax2, ...
    dateIL, ILpred(1:nIL), ...
    '-', ...
    'Color',fitColor, ...
    'LineWidth',fitLineWidth);

xlabel(ax2,'Year', ...
    'FontName',fontName, ...
    'FontSize',fontSizeLabel);

ylabel(ax2,'Estimated influenza cases', ...
    'FontName',fontName, ...
    'FontSize',fontSizeLabel);
ylim([0 2500])
text(ax2,0.03,0.95,'(b) IL', ...
    'Units','normalized', ...
    'FontName',fontName, ...
    'FontSize',fontSizePanel, ...
    'FontWeight','bold', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','top');

% Put one shared legend in panel (b)
lgd = legend(ax2, ...
    {'Real data','Fitting curve'}, ...
    'FontName',fontName, ...
    'FontSize',fontSizeLegend, ...
    'Location','northeast', ...
    'Box','off');
lgd.Units = 'normalized';
lgd.Position = [0.5 0.82 0.14 0.10];
hold(ax2,'off');

%% =========================================================
% (c) IN
% ==========================================================

ax3 = nexttile(tl,3);
hold(ax3,'on');

plot(ax3, ...
    dateIN, INo(1:nIN), ...
    'o', ...
    'LineStyle','none', ...
    'MarkerSize',markerSize, ...
    'MarkerEdgeColor',obsColor, ...
    'MarkerFaceColor',obsColor);

plot(ax3, ...
    dateIN, INpred(1:nIN), ...
    '-', ...
    'Color',fitColor, ...
    'LineWidth',fitLineWidth);

xlabel(ax3,'Year', ...
    'FontName',fontName, ...
    'FontSize',fontSizeLabel);

ylabel(ax3,'Estimated influenza cases', ...
    'FontName',fontName, ...
    'FontSize',fontSizeLabel);

text(ax3,0.03,0.95,'(c) IN', ...
    'Units','normalized', ...
    'FontName',fontName, ...
    'FontSize',fontSizePanel, ...
    'FontWeight','bold', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','top');
ylim([0 500])
lgd3 = legend(ax3, ...
    {'Real data','Fitting curve'}, ...
    'FontName',fontName, ...
    'FontSize',fontSizeLegend, ...
    'Box','off');

lgd3.Units = 'normalized';
lgd3.Position = [0.79 0.82 0.14 0.10];
hold(ax3,'off');

%% =========================================================
% Common formatting
% ==========================================================

axAll = [ax1 ax2 ax3];

for i = 1:3

    axAll(i).FontName = fontName;
    axAll(i).FontSize = fontSizeTick;

    axAll(i).LineWidth = axisLineWidth;
    axAll(i).TickDir = 'out';
    axAll(i).Box = 'off';

    axAll(i).XGrid = 'off';
    axAll(i).YGrid = 'on';

    axAll(i).GridAlpha = 0.10;
    axAll(i).GridLineStyle = '-';

    xticks(axAll(i),yearTicks);
    xtickformat(axAll(i),'yyyy');

end

%% Keep time axes aligned

linkaxes(axAll,'x');

%% =========================================================
% Export
% ==========================================================

drawnow;

exportgraphics(fig, ...
    'Figure8.eps', ...
    'ContentType','vector', ...
    'BackgroundColor','white');