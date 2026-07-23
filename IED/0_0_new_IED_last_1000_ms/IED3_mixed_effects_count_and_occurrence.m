% Test within-patient log-transformed IED count and IED occurrence
% using pooled mixed-effects models

% Author: Nill

clear;
clc;
close all;

inputFolder = 'D:\Nill\data\BART\0_0_new_IED_last_1000_ms\IED1_find_number_of_IEDs\';

outputFolder = ...
    'D:\Nill\code\BART\IED\0_0_new_IED_last_1000_ms\IED3_mixed_effects_count_and_occurrence\';

mkdir(outputFolder);

files = dir(fullfile(inputFolder, '*.LFPIED.mat'));
allTrials = table();

% Put all participants and trials together

for pt = 1:length(files)

    fileName = files(pt).name;
    parts = strsplit(fileName, '.');
    ptID = string(parts{1});

    disp("Processing patient: " + ptID)

    data = load(fullfile(inputFolder, fileName));
    LFPIED = data.LFPIED;

    nTrials = LFPIED.nTrials;

    RT = LFPIED.RTs(1:nTrials);
    IT = LFPIED.ITs(1:nTrials);
    BR = LFPIED.BankedTrials(1:nTrials);
    control = LFPIED.isControl(1:nTrials);

    RT = RT(:);
    IT = IT(:);
    BR = BR(:);
    control = control(:);

    % Trial number for each IED

    RTtrials = round(LFPIED.IED_occurance_RT(:, 1));
    ITtrials = round(LFPIED.IED_occurance_IT(:, 1));

    RTtrials = RTtrials(isfinite(RTtrials) & ...
        RTtrials >= 1 & RTtrials <= nTrials);

    ITtrials = ITtrials(isfinite(ITtrials) & ...
        ITtrials >= 1 & ITtrials <= nTrials);

    % Count IEDs in each trial

    nIED_RT = accumarray(RTtrials, 1, [nTrials 1]);
    nIED_IT = accumarray(ITtrials, 1, [nTrials 1]);

    % Log-transform IED counts

    logIED_RT = log(nIED_RT + 1);
    logIED_IT = log(nIED_IT + 1);

    % IED occurrence: zero or one

    IEDoccurred_RT = double(nIED_RT > 0);
    IEDoccurred_IT = double(nIED_IT > 0);

    patientID = repmat(ptID, nTrials, 1);

    patientTrials = table(patientID, RT, IT, BR, control, ...
        nIED_RT, nIED_IT, logIED_RT, logIED_IT, ...
        IEDoccurred_RT, IEDoccurred_IT);

    allTrials = [allTrials; patientTrials];
end

allTrials.patientID = categorical(allTrials.patientID);

% Run the six within-patient mixed-effects models

analysisType = ["Log IED count"; ...
                "Log IED count"; ...
                "Log IED count"; ...
                "IED occurrence"; ...
                "IED occurrence"; ...
                "IED occurrence"];

outcome = ["RT"; "IT"; "BR"; "RT"; "IT"; "BR"];

results = table();

for a = 1:6

    % Select predictor and outcome

    if a == 1

        % Log IED count during RT predicting RT

        x = allTrials.logIED_RT;
        y = allTrials.RT;

    elseif a == 2

        % Log IED count during IT predicting IT

        x = allTrials.logIED_IT;
        y = allTrials.IT;

    elseif a == 3

        % Log IED count during IT predicting banking

        x = allTrials.logIED_IT;
        y = allTrials.BR;

    elseif a == 4

        % IED occurrence during RT predicting RT

        x = allTrials.IEDoccurred_RT;
        y = allTrials.RT;

    elseif a == 5

        % IED occurrence during IT predicting IT

        x = allTrials.IEDoccurred_IT;
        y = allTrials.IT;

    else

        % IED occurrence during IT predicting banking

        x = allTrials.IEDoccurred_IT;
        y = allTrials.BR;
    end

    % Keep valid non-control trials

    keep = allTrials.control == 0 & ...
        isfinite(allTrials.RT) & ...
        allTrials.RT > 0 & ...
        allTrials.RT <= 30 & ...
        isfinite(x) & ...
        isfinite(y);

    if outcome(a) == "IT"

        keep = keep & allTrials.IT > 0;

    elseif outcome(a) == "BR"

        keep = keep & ...
            allTrials.IT > 0 & ...
            (allTrials.BR == 0 | allTrials.BR == 1);
    end

    T = table(y(keep), x(keep), allTrials.patientID(keep), ...
        'VariableNames', {'y', 'x', 'patientID'});

    % Calculate each patient's mean predictor

    groups = findgroups(T.patientID);

    patientMean = splitapply(@mean, T.x, groups);

    T.xMean = patientMean(groups);

    % Center predictor around each patient's own mean

    T.xWithin = T.x - T.xMean;

    disp("Running within-patient " + analysisType(a) + ...
        " model for " + outcome(a))

    % Fit mixed-effects model

    if outcome(a) == "BR"

        % Logistic mixed-effects model for banking

        model = fitglme(T, ...
            ['y ~ xWithin + xMean + ' ...
             '(1 | patientID) + ' ...
             '(xWithin - 1 | patientID)'], ...
            'Distribution', 'Binomial', ...
            'Link', 'logit');

    else

        % Linear mixed-effects model for RT and IT

        model = fitlme(T, ...
            ['y ~ xWithin + xMean + ' ...
             '(1 | patientID) + ' ...
             '(xWithin - 1 | patientID)']);
    end

    % Get the within-patient coefficient

    coefficients = model.Coefficients;
    CI = coefCI(model);

    slopeRow = find(strcmp(coefficients.Name, 'xWithin'));

    beta = coefficients.Estimate(slopeRow);
    SE = coefficients.SE(slopeRow);
    testStatistic = coefficients.tStat(slopeRow);
    pValue = coefficients.pValue(slopeRow);

    CILow = CI(slopeRow, 1);
    CIHigh = CI(slopeRow, 2);

    % Calculate odds ratio for banking

    oddsRatio = NaN;
    oddsLow = NaN;
    oddsHigh = NaN;

    if outcome(a) == "BR"

        oddsRatio = exp(beta);
        oddsLow = exp(CILow);
        oddsHigh = exp(CIHigh);
    end

    % Save within-patient results

    nPatients = length(unique(T.patientID));
    nTrialsUsed = height(T);
    effectType = "Within-patient";

    newRow = table(analysisType(a), outcome(a), ...
        effectType, nPatients, nTrialsUsed, ...
        beta, SE, testStatistic, pValue, ...
        CILow, CIHigh, oddsRatio, oddsLow, oddsHigh, ...
        'VariableNames', ...
        {'analysisType', 'outcome', 'effectType', ...
         'nPatients', 'nTrials', ...
         'beta', 'SE', ...
         'testStatistic', 'pValue', ...
         'CILow', 'CIHigh', ...
         'oddsRatio', ...
         'oddsRatioCILow', ...
         'oddsRatioCIHigh'});

    results = [results; newRow];
end

% FDR correction across all six models

p = results.pValue;

[sortedP, order] = sort(p);

fdr = sortedP .* length(p) ./ (1:length(p))';

for a = length(fdr)-1:-1:1

    fdr(a) = min(fdr(a), fdr(a + 1));
end

fdr(fdr > 1) = 1;

correctedP = zeros(length(p), 1);
correctedP(order) = fdr;

results.pValue_FDR = correctedP;
results.significant_FDR = correctedP < 0.05;

% Save results

writetable(results, ...
    fullfile(outputFolder, ...
    'within_patient_mixed_effects_count_and_occurrence.csv'));


colors = [0.204 0.459 0.702; ...
          0.847 0.333 0.153; ...
          0.250 0.600 0.250];

figure('Color', 'w', 'Position', [100 100 750 750]);

for panel = 1:2

    subplot(2, 1, panel)
    hold on

    if panel == 1

        rows = 1:3;
        panelTitle = 'Within-patient log(IED count + 1)';
        yLabelText = 'Within-patient coefficient';

    else

        rows = 4:6;
        panelTitle = 'Within-patient IED occurrence';
        yLabelText = 'Within-patient coefficient';
    end

    beta = results.beta(rows);
    lowError = beta - results.CILow(rows);
    highError = results.CIHigh(rows) - beta;

    % Draw coefficients

    for k = 1:3

        errorbar(k, beta(k), ...
            lowError(k), highError(k), 'o', ...
            'Color', colors(k, :), ...
            'MarkerFaceColor', colors(k, :), ...
            'MarkerSize', 7, ...
            'LineWidth', 1.5, ...
            'CapSize', 10);
    end

    yline(0, '--', ...
        'Color', [0.45 0.45 0.45], ...
        'LineWidth', 1);

    % Put one star on FDR-significant coefficients

    yRange = max(results.CIHigh(rows)) - ...
        min(results.CILow(rows));

    if yRange == 0

        yRange = 1;
    end

    for k = 1:3

        if results.pValue_FDR(rows(k)) < 0.05

            text(k, ...
                results.CIHigh(rows(k)) + 0.03*yRange, ...
                '*', ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom', ...
                'FontSize', 16, ...
                'FontWeight', 'bold');
        end
    end

    xlim([0.5 3.5]);
    xticks(1:3);
    xticklabels({'RT', 'IT', 'BR'});

    ylabel(yLabelText, ...
        'FontSize', 10, ...
        'FontWeight', 'bold');

    title(panelTitle, ...
        'FontSize', 13, ...
        'FontWeight', 'bold');

    set(gca, ...
        'FontSize', 11, ...
        'LineWidth', 1.1, ...
        'Box', 'off', ...
        'TickDir', 'out');
end

% Save figure

exportgraphics(gcf, ...
    fullfile(outputFolder, ...
    'within_patient_mixed_effects_count_and_occurrence.pdf'), ...
    'ContentType', 'vector');

close(gcf)

