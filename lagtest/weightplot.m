clear; clc; close all;

%% =========================================
% File names - two categories
% =========================================
files_3 = {'weight7_3.csv', 'weight10_3.csv', 'weight14_3.csv', 'weight21_3.csv', 'weight28_3.csv'};
files_2 = {'weight7_2.csv', 'weight10_2.csv', 'weight14_2.csv', 'weight21_2.csv', 'weight28_2.csv'};

Lmax_list = [7, 10, 14, 21, 28];
nGroup = numel(files_3);

%% =========================================
% Function to process one group of files
% =========================================
function process_and_plot(files, Lmax_list, ax, label_str)
    nGroup = numel(files);
    CW  = cell(nGroup,1);
    Lag = cell(nGroup,1);
    C7  = nan(nGroup,1);
    L99 = nan(nGroup,1);

    for i = 1:nGroup
        T = readtable(files{i}, 'VariableNamingRule', 'preserve');
        w = T.mean(:);

        CW{i}  = cumsum(w);
        Lag{i} = (0:length(w)-1)';

        if length(w) >= 8
            C7(i) = sum(w(1:8));   % lag 0--7
        else
            C7(i) = sum(w);
        end

        idx = find(CW{i} >= 0.99, 1, 'first');
        L99(i) = Lag{i}(idx);

        fprintf('%s: Lmax = %2d | L99 = %2d | C(7) = %.5f (%.3f%%)\n', ...
            label_str, Lmax_list(i), L99(i), C7(i), 100*C7(i));
    end

    %% Plot on given axes
    hold(ax, 'on');

    for i = 1:nGroup
        plot(ax, Lag{i}, CW{i}, '-o', ...
            'LineWidth', 2.0, ...
            'MarkerSize', 6, ...
            'MarkerFaceColor', 'auto', ...
            'DisplayName', sprintf('L_{max} = %d', Lmax_list(i)));
    end

    %% Reference lines
    yline(ax, 0.99, '--', '99% threshold', ...
        'LineWidth', 1.4, ...
        'FontSize', 11, ...
        'LabelHorizontalAlignment', 'left',...
        'HandleVisibility','off');

    xline(ax, 7, '--', 'Lag 7', ...
        'LineWidth', 1.4, ...
        'FontSize', 11, ...
        'LabelVerticalAlignment', 'bottom',...
        'HandleVisibility','off');

    %% Axis labels
    xlabel(ax, 'Lag (days)', 'FontSize', 13);
    ylabel(ax, 'Cumulative lag weight, C(k)', 'FontSize', 13);

    %% Axis limits
    xlim(ax, [0 12]);
    xticks(ax, 0:2:12);

    ylim(ax, [0.85 1.01]);
    yticks(ax, [0.86 0.88 0.90 0.92 0.94 0.96 0.98 0.99 1.00]);

    %% Axis style
    box(ax, 'off');
    set(ax, ...
        'FontName', 'Arial', ...
        'FontSize', 12, ...
        'LineWidth', 1.2, ...
        'TickDir', 'out', ...
        'TickLength', [0.015 0.015]);

    %% Legend
    lgd = legend(ax, 'Location', 'southeast', 'Box', 'off');
    lgd.FontSize = 11;

    %% Add (a) or (b) label in top-left corner
    text(ax, 0.02, 0.97, label_str, ...
        'Units', 'normalized', ...
        'FontSize', 14, ...
        'FontWeight', 'bold', ...
        'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'left',...
        Position=[0.1,1.04]);
end

%% =========================================
% Create figure with two subplots
% =========================================
figure('Units', 'centimeters', ...
       'Position', [2 2 22 9], ...
       'Color', 'w');

% Left subplot - weightX_3
ax1 = subplot(1, 2, 1);
process_and_plot(files_3, Lmax_list, ax1, '(a)');

% Right subplot - weightX_2
ax2 = subplot(1, 2, 2);
process_and_plot(files_2, Lmax_list, ax2, '(b)');