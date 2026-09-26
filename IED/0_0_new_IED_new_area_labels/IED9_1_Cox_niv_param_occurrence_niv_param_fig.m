clear;
clc;
close all;

nivFile = ...
    'D:\Nill\code\BART\IED\0_0_new_IED_new_area_labels\context_modeling\param_recovery_4_param_recovery\alpha_comparison.csv';

outputFolder = ...
    'D:\Nill\code\BART\IED\0_0_new_IED_new_area_labels\IED9_1_Cox_niv_param_occurrence\';

expectedParticipants = 71;

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

%% Read all participant rows without filtering or deduplication
raw = readtable(nivFile, 'TextType', 'string', ...
    'VariableNamingRule', 'preserve');

requiredVariables = ["ptID", "fit_alpha_plus", "fit_alpha_minus"];
missingVariables = requiredVariables(~ismember(requiredVariables, ...
    string(raw.Properties.VariableNames)));

if ~isempty(missingVariables)
    error('Missing required CSV column(s): %s', ...
        strjoin(missingVariables, ', '));
end

numberOfParticipants = height(raw);

if numberOfParticipants ~= expectedParticipants
    error(['Expected %d participant rows in alpha_comparison.csv, ' ...
        'but found %d. No rows were removed by this script.'], ...
        expectedParticipants, numberOfParticipants);
end

patientID = strip(string(raw.("ptID")));
alphaPlus = str2double(string(raw.("fit_alpha_plus")));
alphaMinus = str2double(string(raw.("fit_alpha_minus")));

blankID = strlength(patientID) == 0 | ismissing(patientID);
patientID(blankID) = "Missing ID (row " + find(blankID) + ")";

denominator = alphaMinus + alphaPlus;
validAlpha = isfinite(alphaPlus) & isfinite(alphaMinus);
validNiv = validAlpha & denominator > 0;
invalidNiv = ~validNiv;

nivRaw = NaN(numberOfParticipants, 1);
nivRaw(validNiv) = ...
    (alphaMinus(validNiv) - alphaPlus(validNiv)) ./ ...
    denominator(validNiv);

fprintf('Participant rows retained for visualization: %d\n', ...
    numberOfParticipants);
fprintf('Participants with valid Niv values: %d\n', sum(validNiv));
fprintf('Participants shown as N/A: %d\n', sum(invalidNiv));

%% Colors
maximumAbsoluteNiv = max(abs(nivRaw(validNiv)), [], 'omitnan');

if isempty(maximumAbsoluteNiv) || ~isfinite(maximumAbsoluteNiv) || ...
        maximumAbsoluteNiv == 0
    maximumAbsoluteNiv = 1;
end

colorLimits = [-maximumAbsoluteNiv maximumAbsoluteNiv];

blue = [0.12 0.35 0.75];
white = [0.97 0.97 0.97];
red = [0.78 0.18 0.18];
gray = [0.45 0.45 0.45];

colorMap = [ ...
    linspace(blue(1), white(1), 128)', ...
    linspace(blue(2), white(2), 128)', ...
    linspace(blue(3), white(3), 128)'; ...
    linspace(white(1), red(1), 128)', ...
    linspace(white(2), red(2), 128)', ...
    linspace(white(3), red(3), 128)' ...
    ];

%% Wide figure and layout
fig = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [40 40 2600 900]);

set(fig, 'Colormap', colorMap);

layout = tiledlayout(fig, 1, 3, ...
    'Padding', 'compact', ...
    'TileSpacing', 'compact');

%% Left subplot: alpha+ versus alpha-
ax1 = nexttile(layout, 1);
hold(ax1, 'on');

finiteAlphaValues = [alphaPlus(validAlpha); alphaMinus(validAlpha)];
axisLimit = max(finiteAlphaValues, [], 'omitnan');

if isempty(axisLimit) || ~isfinite(axisLimit) || axisLimit <= 0
    axisLimit = 1;
else
    axisLimit = 1.08 * axisLimit;
end

plot(ax1, [0 axisLimit], [0 axisLimit], '--', ...
    'Color', [0.35 0.35 0.35], ...
    'LineWidth', 1.4);

scatter(ax1, alphaPlus(validNiv), alphaMinus(validNiv), ...
    92, nivRaw(validNiv), 'filled', ...
    'MarkerFaceAlpha', 0.75, ...
    'MarkerEdgeColor', 'none');

plottableInvalid = validAlpha & invalidNiv;

if any(plottableInvalid)
    scatter(ax1, alphaPlus(plottableInvalid), ...
        alphaMinus(plottableInvalid), ...
        92, gray, 'x', 'LineWidth', 1.5);
end

set(ax1, 'CLim', colorLimits);

xlim(ax1, [0 axisLimit]);
ylim(ax1, [0 axisLimit]);
axis(ax1, 'square');

xlabel(ax1, 'Fitted positive learning rate, alpha+');
ylabel(ax1, 'Fitted negative learning rate, alpha-');

title(ax1, 'alpha+ versus alpha- by participant', ...
    'FontWeight', 'bold');

colorbar(ax1);
box(ax1, 'off');
hold(ax1, 'off');

%% Right subplot: wide rectangle with participants on x-axis
ax2 = nexttile(layout, 2, [1 2]);
hold(ax2, 'on');

validRows = find(validNiv);
[nivSorted, validSortOrder] = sort(nivRaw(validRows));
sortedValidRows = validRows(validSortOrder);

invalidRows = find(invalidNiv);

numberOfValidParticipants = numel(sortedValidRows);
validPositions = (1:numberOfValidParticipants)';
invalidPositions = ...
    (numberOfValidParticipants + 1:numberOfParticipants)';

scatter(ax2, validPositions, nivSorted, ...
    62, nivSorted, 'filled', ...
    'MarkerFaceAlpha', 0.85, ...
    'MarkerEdgeColor', 'none');

if ~isempty(invalidPositions)
    scatter(ax2, invalidPositions, ...
        zeros(size(invalidPositions)), ...
        62, gray, 'x', 'LineWidth', 1.5);

    text(ax2, invalidPositions, ...
        0.04 * maximumAbsoluteNiv * ...
        ones(size(invalidPositions)), ...
        repmat("N/A", size(invalidPositions)), ...
        'Color', gray, ...
        'FontSize', 8, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom');
end

yline(ax2, 0, '--', ...
    'Color', [0.25 0.25 0.25]);

set(ax2, 'CLim', colorLimits);

xlim(ax2, [0.5 numberOfParticipants + 0.5]);
ylim(ax2, 1.08 * colorLimits);

% Explicitly make the right plot box a wide rectangle:
% x-axis length is three times the y-axis length.
pbaspect(ax2, [1.5 1 1]);

participantTicks = unique([1, 5:5:numberOfParticipants, ...
    numberOfParticipants]);
xticks(ax2, participantTicks);

xlabel(ax2, sprintf( ...
    'Participants sorted by Niv (n = %d)', ...
    numberOfParticipants));

ylabel(ax2, 'Niv learning-rate asymmetry');

title(ax2, sprintf( ...
    'Participant-level Niv asymmetry (%d valid, %d N/A)', ...
    sum(validNiv), sum(invalidNiv)), ...
    'FontWeight', 'bold');

box(ax2, 'off');
hold(ax2, 'off');

title(layout, 'Fitted learning rates and Niv asymmetry', ...
    'FontWeight', 'bold', ...
    'FontSize', 17);

%% Export visualization
outputFile = fullfile(outputFolder, ...
    'Niv_alpha_plus_alpha_minus_participant_visualization.pdf');

exportgraphics(fig, outputFile, 'ContentType', 'vector');
close(fig);

fprintf('Saved visualization: %s\n', outputFile);