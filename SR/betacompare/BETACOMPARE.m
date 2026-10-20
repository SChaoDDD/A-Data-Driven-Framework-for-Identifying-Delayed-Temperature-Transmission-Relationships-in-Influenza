clear;
clc;
close all;

%% =========================================================
% Comparison of transmission-rate reconstruction methods
% IA / IL / IN
% ==========================================================

%% 1. Read data

IA = readtable("IAbeta.csv");
IL = readtable("ILbeta.csv");
IN = readtable("INbeta.csv");

% Weekly UDE-inferred beta
IAbaseline = table2array(IA(1:7:end,3));
ILbaseline = table2array(IL(1:7:end,3));
INbaseline = table2array(IN(1:7:end,3));

%% HDLNM-SR

tmp = load("IAbetafit.mat");
IAdlnm = tmp.beta_pred;

tmp = load("ILbetafit.mat");
ILdlnm = tmp.beta_pred;

tmp = load("INbetafit.mat");
INdlnm = tmp.beta_pred;

% Preserve original construction
IA0 = [IAbaseline(1:30); IAdlnm];
IL0 = [ILbaseline(1:30); ILdlnm];
IN0 = [INbaseline(1:30); INdlnm];

% Weekly HDLNM-SR values
IAHDLNM = IA0(1:7:end);
ILHDLNM = IL0(1:7:end);
INHDLNM = IN0(1:7:end);

%% Previous method 1

tmp = load("IAmethod1.mat");
IAmethod1 = tmp.Beta1;

tmp = load("ILmethod1.mat");
ILmethod1 = tmp.Beta1;

tmp = load("INmethod1.mat");
INmethod1 = tmp.Beta1;

%% Previous method 2

tmp = load("IAmethod2.mat");
IAmethod2 = tmp.Beta2;

tmp = load("ILmethod2.mat");
ILmethod2 = tmp.Beta2;

tmp = load("INmethod2.mat");
INmethod2 = tmp.Beta2;

%% Ensemble uncertainty

IAlow = table2array(IA(1:7:end,5));
IAup  = table2array(IA(1:7:end,6));

ILlow = table2array(IL(1:7:end,5));
ILup  = table2array(IL(1:7:end,6));

INlow = table2array(IN(1:7:end,5));
INup  = table2array(IN(1:7:end,6));


%% =========================================================
% 2. Time axis
% ==========================================================

startDate = datetime(2013,9,30);

% Maximum weekly time axis
maxWeeks = 312;

dateWeekly = startDate + calweeks(0:maxWeeks-1);

% Year ticks
yearTicks = datetime(2014:2019,1,1);


%% =========================================================
% 3. Plotting parameters
% ==========================================================

fontName       = 'Times New Roman';
fontSizeTick   = 9;
fontSizeLabel  = 10;
fontSizePanel  = 13;
fontSizeLegend = 9;

axisLineWidth = 0.8;

% Main methods
udeLineWidth = 2.5;
srLineWidth  = 1.3;

% Previous methods
previousLineWidth = 1.2;

% Preserve original colors
udeColor     = [0.10 0.50 0.60];
srColor      = '#e17000';
method1Color = '#4cd09a';
method2Color = '#ffb3dd';

% Uncertainty band
bandColor = [0.882 0.435 0];
bandAlpha = 0.18;

% Figure size
figWidth  = 20.0;
figHeight = 7.2;


%% =========================================================
% 4. Individual y-axis limits
% Modify these yourself
% ==========================================================

ylimIA = [0.3 2.5];
ylimIL = [0.3 2.5];
ylimIN = [0 3.5];

% Store together for loop
yLimits = {
    ylimIA
    ylimIL
    ylimIN
};


%% =========================================================
% 5. Organize data
% ==========================================================

baseline = {
    IAbaseline
    ILbaseline
    INbaseline
};

HDLNM = {
    IAHDLNM
    ILHDLNM
    INHDLNM
};

method1 = {
    IAmethod1
    ILmethod1
    INmethod1
};

method2 = {
    IAmethod2
    ILmethod2
    INmethod2
};

lowerBand = {
    IAlow
    ILlow
    INlow
};

upperBand = {
    IAup
    ILup
    INup
};

regionNames = {'IA','IL','IN'};
panelLabels = {'(a)','(b)','(c)'};


%% =========================================================
% 6. Create figure
% ==========================================================

fig = figure( ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[2 2 figWidth figHeight]);

tl = tiledlayout(fig,1,3, ...
    'TileSpacing','compact', ...
    'Padding','compact');

ax = gobjects(3,1);

legendHandles = cell(3,1);


%% 7. Plot
% ==========================================================

for i = 1:3

    ax(i) = nexttile(tl,i);

    hold(ax(i),'on');

    %% Determine common data length

    n = min([ ...
        length(baseline{i}), ...
        length(HDLNM{i}), ...
        length(method1{i}), ...
        length(method2{i}), ...
        length(lowerBand{i}), ...
        length(upperBand{i}), ...
        length(dateWeekly) ...
        ]);

    datePlot = dateWeekly(1:n);


    %% Ensemble uncertainty first

    fill(ax(i), ...
        [datePlot, fliplr(datePlot)], ...
        [lowerBand{i}(1:n)', ...
         fliplr(upperBand{i}(1:n)')], ...
        bandColor, ...
        'FaceAlpha',bandAlpha, ...
        'EdgeColor','none', ...
        'HandleVisibility','off');


    %% UDE

    hUDE = plot(ax(i), ...
        datePlot, ...
        baseline{i}(1:n), ...
        '-', ...
        'LineWidth',udeLineWidth, ...
        'Color',udeColor);


    %% Equation [52]

    hM1 = plot(ax(i), ...
        datePlot, ...
        method1{i}(1:n), ...
        '-', ...
        'LineWidth',previousLineWidth, ...
        'Color',method1Color);


    %% Equation [54]

    hM2 = plot(ax(i), ...
        datePlot, ...
        method2{i}(1:n), ...
        '-', ...
        'LineWidth',previousLineWidth, ...
        'Color',method2Color);


    %% HDLNM-SR (最后绘制)

    hSR = plot(ax(i), ...
        datePlot, ...
        HDLNM{i}(1:n), ...
        '-', ...
        'LineWidth',srLineWidth, ...
        'Color',srColor);


    legendHandles{i} = [hUDE hSR hM1 hM2];


    %% Axis labels

    xlabel(ax(i), ...
        'Year', ...
        'FontName',fontName, ...
        'FontSize',fontSizeLabel);

    ylabel(ax(i), ...
        'Transmission rate', ...
        'FontName',fontName, ...
        'FontSize',fontSizeLabel);


    %% X axis

    xlim(ax(i), ...
        [datePlot(1) datePlot(end)]);

    xticks(ax(i),yearTicks);

    xtickformat(ax(i),'yyyy');


    %% Individual Y axis

    ylim(ax(i),yLimits{i});


    %% Panel labels

    text(ax(i), ...
        0.03,0.99, ...
        sprintf('%s %s', ...
        panelLabels{i}, ...
        regionNames{i}), ...
        'Units','normalized', ...
        'FontName',fontName, ...
        'FontSize',fontSizePanel, ...
        'FontWeight','bold', ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','top', ...
        'Color','k');

    hold(ax(i),'off');

end


%% =========================================================
% 8. Common formatting
% ==========================================================

for i = 1:3

    ax(i).FontName = fontName;
    ax(i).FontSize = fontSizeTick;

    ax(i).LineWidth = axisLineWidth;

    ax(i).TickDir = 'out';
    ax(i).Box = 'off';

    ax(i).XGrid = 'off';
    ax(i).YGrid = 'on';

    ax(i).GridAlpha = 0.08;
    ax(i).GridLineStyle = '-';

end


%% =========================================================
% 9. Legends
% Each legend can be moved independently
% ==========================================================

%% IA

lgdIA = legend(ax(1), ...
    legendHandles{1}, ...
    {'UDE', ...
     'Symbolic regression', ...
     'Formulation [42]', ...
     'Formulation [44]'}, ...
    'FontName',fontName, ...
    'FontSize',fontSizeLegend, ...
    'Box','off');

lgdIA.Units = 'normalized';

% Modify this yourself
lgdIA.Position = [0.2 0.79 0.12 0.17];


%% IL

lgdIL = legend(ax(2), ...
    legendHandles{2}, ...
    {'UDE', ...
     'Symbolic regression', ...
     'Formulation [42]', ...
     'Formulation [44]'}, ...
    'FontName',fontName, ...
    'FontSize',fontSizeLegend, ...
    'Box','off');

lgdIL.Units = 'normalized';

% Modify this yourself
lgdIL.Position = [0.5 0.79 0.12 0.17];


%% IN

lgdIN = legend(ax(3), ...
    legendHandles{3}, ...
    {'UDE', ...
     'Symbolic regression', ...
     'Formulation [42]', ...
     'Formulation [44]'}, ...
    'FontName',fontName, ...
    'FontSize',fontSizeLegend, ...
    'Box','off');

lgdIN.Units = 'normalized';

% Modify this yourself
lgdIN.Position = [0.82 0.79 0.12 0.17];


%% =========================================================
% 10. Export
% ==========================================================

drawnow;

exportgraphics(fig, ...
    'Figure10.eps', ...
    'ContentType','vector', ...
    'BackgroundColor','white');
% 
% observed = IAbaseline;
% predicted = IAmethod2';

%% =======================
%  方法一：手动计算
% ========================

% % 计算残差
% residuals = observed - predicted;
% 
% % 计算 RMSE
% n = length(observed);
% rmse_manual = sqrt(sum(residuals.^2) / n);
% 
% % 计算 R²
% SS_res = sum(residuals.^2);                    % 残差平方和
% SS_tot = sum((observed - mean(observed)).^2);  % 总平方和
% R2_manual = 1 - SS_res / SS_tot;
% 
% fprintf('=== 手动计算结果 ===\n');
% fprintf('RMSE = %.4f\n', rmse_manual);
% fprintf('R²   = %.4f\n', R2_manual);