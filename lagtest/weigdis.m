clear; clc; close all;

%% =========================================================
% Figure: Bayesian lag analysis
% (a) Cumulative lag weights under different Lmax
% (b) Posterior mean lag weights for L = 7
% ==========================================================

%% 1. Plotting parameters

fontName       = 'Times New Roman';
fontSizeTick   = 10;
fontSizeLabel  = 10;
fontSizePanel  = 13;
fontSizeLegend = 9;

lineWidth      = 1.8;
markerSize     = 5.5;
axisLineWidth  = 0.8;

% Bar color for panel (b)
barColor = [0.32 0.52 0.68];

% Figure size
figWidth  = 17.5;   % cm
figHeight = 7.5;    % cm


%% 2. File names

files = { ...
    'weight7.csv', ...
    'weight10.csv', ...
    'weight14.csv', ...
    'weight21.csv', ...
    'weight28.csv' ...
    };

Lmax_list = [7, 10, 14, 21, 28];
nGroup = numel(files);

CW  = cell(nGroup,1);
Lag = cell(nGroup,1);

C7  = nan(nGroup,1);
L99 = nan(nGroup,1);


%% 3. Read data and calculate cumulative weights

for i = 1:nGroup

    T = readtable(files{i}, ...
        'VariableNamingRule','preserve');

    w = T.mean(:);

    CW{i}  = cumsum(w);
    Lag{i} = (0:length(w)-1)';

    if length(w) >= 8
        C7(i) = sum(w(1:8));      % lag 0--7
    else
        C7(i) = sum(w);
    end

    idx = find(CW{i} >= 0.99, 1, 'first');

    if ~isempty(idx)
        L99(i) = Lag{i}(idx);
    end

    fprintf( ...
        'Lmax = %2d | L99 = %2g | C(7) = %.5f (%.3f%%)\n', ...
        Lmax_list(i), L99(i), C7(i), 100*C7(i));

end


%% 4. Read posterior mean lag weights for L = 7

we = readtable('weight7.csv', ...
    'VariableNamingRule','preserve');

lag = (0:7)';

% Preserve your original choice: second column
weight = table2array(we(:,2));

% Safety: use only lag 0--7
weight = weight(1:min(8,length(weight)));
lag    = lag(1:length(weight));


%% 5. Create figure

fig = figure( ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[2 2 figWidth figHeight]);

tl = tiledlayout(fig,1,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');


%% =========================================================
% Panel (a): cumulative lag weights
% ==========================================================

ax1 = nexttile(tl,1);
hold(ax1,'on');

for i = 1:nGroup

    % Do not manually specify Color:
    % preserves MATLAB's current default color order
    plot(ax1, ...
        Lag{i}, CW{i}, ...
        '-o', ...
        'LineWidth',lineWidth, ...
        'MarkerSize',markerSize, ...
        'MarkerFaceColor','auto', ...
        'DisplayName', ...
        sprintf('L_{max} = %d',Lmax_list(i)));

end

% Reference lines
yline(ax1,0.99,'--', ...
    '99% threshold', ...
    'LineWidth',1.1, ...
    'Color',[0.35 0.35 0.35], ...
    'FontName',fontName, ...
    'FontSize',fontSizeTick, ...
    'LabelHorizontalAlignment','left', ...
    'HandleVisibility','off');

xline(ax1,7,'--', ...
    'Lag 7', ...
    'LineWidth',1.1, ...
    'Color',[0.35 0.35 0.35], ...
    'FontName',fontName, ...
    'FontSize',fontSizeTick, ...
    'LabelVerticalAlignment','bottom', ...
    'HandleVisibility','off');

xlabel(ax1, ...
    'Lag (days)', ...
    'FontName',fontName, ...
    'FontSize',fontSizeLabel);

ylabel(ax1, ...
    'Cumulative lag weight, C(k)', ...
    'FontName',fontName, ...
    'FontSize',fontSizeLabel);

xlim(ax1,[0 12]);
xticks(ax1,0:2:12);

ylim(ax1,[0.90 1.002]);
yticks(ax1,[0.90 0.92 0.94 0.96 0.98 0.99 1.00]);

% Panel label
text(ax1, ...
    0.03,1.04, ...
    '(a)', ...
    'Units','normalized', ...
    'FontName',fontName, ...
    'FontSize',fontSizePanel, ...
    'FontWeight','bold', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','top');

legend(ax1, ...
    'Location','east', ...
    'FontName',fontName, ...
    'FontSize',fontSizeLegend, ...
    'Box','off');

hold(ax1,'off');


%% =========================================================
% Panel (b): posterior mean lag weights for L = 7
% ==========================================================

ax2 = nexttile(tl,2);
hold(ax2,'on');

bar(ax2, ...
    lag, weight, ...
    0.68, ...
    'FaceColor',barColor, ...
    'EdgeColor','none');

xlabel(ax2, ...
    'Lag (days)', ...
    'FontName',fontName, ...
    'FontSize',fontSizeLabel);

ylabel(ax2, ...
    'Posterior mean lag weight', ...
    'FontName',fontName, ...
    'FontSize',fontSizeLabel);

xlim(ax2,[-0.6 7.6]);
xticks(ax2,0:7);

% Make bar panel start from zero
ylim(ax2,[0, max(weight)*1.15]);

% Panel label
text(ax2, ...
    0.03,1.04, ...
    '(b)', ...
    'Units','normalized', ...
    'FontName',fontName, ...
    'FontSize',fontSizePanel, ...
    'FontWeight','bold', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','top');

hold(ax2,'off');


%% 6. Common axis formatting

axAll = [ax1 ax2];

for i = 1:2

    axAll(i).FontName = fontName;
    axAll(i).FontSize = fontSizeTick;

    axAll(i).LineWidth = axisLineWidth;
    axAll(i).TickDir = 'out';
    axAll(i).Box = 'off';

    axAll(i).XGrid = 'off';
    axAll(i).YGrid = 'on';

    axAll(i).GridAlpha = 0.10;
    axAll(i).GridLineStyle = '-';

end


%% 7. Export

drawnow;

exportgraphics(fig, ...
    'Figure6.eps', ...
    'ContentType','vector', ...
    'BackgroundColor','white');