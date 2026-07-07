%% Parameter recovery visualization
% MATLAB version of param_recovery_5_visualization.ipynb

clear;
clc;
close all;

%% Folders and input data
outputFolderName = 'param_recovery_5_visualization';
inputFolderName  = 'param_recovery_4_param_recovery';

if ~exist(outputFolderName, 'dir')
    mkdir(outputFolderName);
end

filePath = fullfile(inputFolderName, 'alpha_comparison.csv');
df = readtable(filePath);

fitAlphaPlus  = double(df.fit_alpha_plus);
fitAlphaMinus = double(df.fit_alpha_minus);
simAlphaPlus  = double(df.sim_alpha_plus);
simAlphaMinus = double(df.sim_alpha_minus);

%% Fitted-versus-simulated parameter-recovery plots
fig1 = figure( ...
    'Color', 'w', ...
    'Units', 'inches', ...
    'Position', [1, 1, 12, 5]);

layout = tiledlayout(fig1, 1, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

ax1 = nexttile(layout, 1);
plotWithRecovery(ax1, fitAlphaPlus, simAlphaPlus, '$\alpha^{+}$');

ax2 = nexttile(layout, 2);
plotWithRecovery(ax2, fitAlphaMinus, simAlphaMinus, '$\alpha^{-}$');

savePath = fullfile(outputFolderName, 'alpha_plus_minus_scatter.pdf');
exportgraphics(fig1, savePath, 'Resolution', 300);

%% Alpha asymmetry
% (alpha_minus - alpha_plus) / (alpha_minus + alpha_plus)
% > 0: risk aversion
% < 0: risk seeking

denom = fitAlphaMinus + fitAlphaPlus;
alphaAsymmetry = nan(size(denom));

valid = denom ~= 0;
alphaAsymmetry(valid) = ...
    (fitAlphaMinus(valid) - fitAlphaPlus(valid)) ./ denom(valid);

fig2 = figure( ...
    'Color', 'w', ...
    'Units', 'inches', ...
    'Position', [1, 1, 10, 5]);

% Python participant indices begin at zero, so preserve that convention.
participantIndex = 0:(numel(alphaAsymmetry) - 1);

plot(participantIndex, alphaAsymmetry, 'o-', ...
    'LineWidth', 1);
hold on;
yline(0, '--', 'LineWidth', 1);
hold off;

xlabel('participant');
ylabel('$(\alpha^- - \alpha^+)/(\alpha^- + \alpha^+)$', ...
    'Interpreter', 'latex');
title('niv');

xticks(participantIndex);
box off;

savePath = fullfile(outputFolderName, ...
    'alpha_asymmetry_by_participant.pdf');
exportgraphics(fig2, savePath, 'Resolution', 300);

%% Local function
function plotWithRecovery(ax, x, y, titleText)
    x = double(x(:));
    y = double(y(:));

    % Remove NaNs, matching the Python notebook.
    validMask = ~isnan(x) & ~isnan(y);
    x = x(validMask);
    y = y(validMask);

    if numel(x) < 2
        text(ax, 0.5, 0.5, 'Not enough valid points', ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'Units', 'normalized');
        title(ax, titleText, 'Interpreter', 'latex');
        return;
    end

    hold(ax, 'on');

    % Scatter plot. This color matches Matplotlib's default blue.
    scatter(ax, x, y, 30, ...
        'MarkerFaceColor', [0.1216, 0.4667, 0.7059], ...
        'MarkerEdgeColor', 'none', ...
        'MarkerFaceAlpha', 0.5);

    % Shared axis limits.
    minVal = min([x; y]);
    maxVal = max([x; y]);
    padding = 0.02;

    % Identity line.
    identityLimits = [minVal - padding, maxVal + padding];
    plot(ax, identityLimits, identityLimits, '--k', ...
        'LineWidth', 0.8);

    % Linear regression: simulated ~ fitted.
    regressionCoefficients = polyfit(x, y, 1);
    slope = regressionCoefficients(1);
    intercept = regressionCoefficients(2); 

    % Pearson correlation.
    correlationMatrix = corrcoef(x, y, 'Rows', 'complete');
    r = correlationMatrix(1, 2);

    % Two-sided p-value for Pearson correlation.
    n = numel(x);
    if n > 2 && abs(r) < 1
        tStatistic = r * sqrt((n - 2) / (1 - r^2));
        p = betainc((n - 2) / ((n - 2) + tStatistic^2), ...
            (n - 2) / 2, 0.5); 
    else
        p = 0; 
    end

    % RMSE relative to the identity line.
    rmseIdentity = sqrt(mean((y - x).^2)); 

    title(ax, ...
        {titleText, sprintf('$r = %.2f,\; \mathrm{slope} = %.2f$', r, slope)}, ...
        'Interpreter', 'latex');

    xlabel(ax, 'Fitted');
    ylabel(ax, 'Simulated');

    xlim(ax, identityLimits);
    ylim(ax, identityLimits);
    pbaspect(ax, [1, 1, 1]);

    box(ax, 'off');
    ax.TickDir = 'out';
    ax.XAxis.TickLabelFormat = '%.1f';
    ax.YAxis.TickLabelFormat = '%.1f';

    hold(ax, 'off');
end
