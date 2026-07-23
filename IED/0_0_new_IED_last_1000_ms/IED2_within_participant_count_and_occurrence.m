% Test IED count and IED occurrence separately for each participant
% Top subplot: IED count
% Bottom subplot: IED occurrence
% RT and IT stay in raw seconds
% BR coefficients are logistic regression coefficients
% Significance is based on permutation p < 0.05
% Author: Nill

clear;
clc;
close all;

inputFolder = 'D:\Nill\data\BART\0_0_new_IED_last_1000_ms\IED1_find_number_of_IEDs\';

outputFolder = ...
    'D:\Nill\code\BART\IED\0_0_new_IED_last_1000_ms\IED2_within_participant_count_and_occurrence\';

mkdir(outputFolder);

nPerm = 1000;

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

    % IED occurrence: zero or one

    IEDoccurred_RT = double(nIED_RT > 0);
    IEDoccurred_IT = double(nIED_IT > 0);

    patientID = repmat(ptID, nTrials, 1);

    patientTrials = table(patientID, RT, IT, BR, control, ...
        nIED_RT, nIED_IT, IEDoccurred_RT, IEDoccurred_IT);

    allTrials = [allTrials; patientTrials];
end

analysisType = ["IED count"; ...
                "IED count"; ...
                "IED count"; ...
                "IED occurrence"; ...
                "IED occurrence"; ...
                "IED occurrence"];

outcome = ["RT"; "IT"; "BR"; "RT"; "IT"; "BR"];

patientList = unique(allTrials.patientID);
results = table();

rng(1);

% Run six models separately for every participant

for pt = 1:length(patientList)

    ptID = patientList(pt);

    disp(" ")
    disp("Running models for patient: " + ptID)

    patientRows = allTrials.patientID == ptID;

    for a = 1:6

        if a == 1

            % IED count during RT predicting RT

            x = allTrials.nIED_RT;
            y = allTrials.RT;

        elseif a == 2

            % IED count during IT predicting IT

            x = allTrials.nIED_IT;
            y = allTrials.IT;

        elseif a == 3

            % IED count during IT predicting banking

            x = allTrials.nIED_IT;
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

        % Keep valid non-control trials for this participant

        keep = patientRows & ...
            allTrials.control == 0 & ...
            isfinite(allTrials.RT) & ...
            allTrials.RT > 0 & ...
            allTrials.RT <= 15 & ...
            isfinite(x) & ...
            isfinite(y);

        if outcome(a) == "IT"

            keep = keep & allTrials.IT > 0;

        elseif outcome(a) == "BR"

            keep = keep & ...
                allTrials.IT > 0 & ...
                (allTrials.BR == 0 | allTrials.BR == 1);
        end

        xPatient = x(keep);
        yPatient = y(keep);

        % A slope cannot be calculated when the predictor does not vary

        if length(xPatient) < 3 || length(unique(xPatient)) < 2

            disp("Skipping " + analysisType(a) + " " + ...
                outcome(a) + ": predictor does not vary")

            continue
        end

        % Logistic regression also needs both banking outcomes

        if outcome(a) == "BR" && length(unique(yPatient)) < 2

            disp("Skipping " + analysisType(a) + ...
                " BR: banking does not vary")

            continue
        end

        T = table(yPatient, xPatient, ...
            'VariableNames', {'y', 'x'});

        disp("Running " + analysisType(a) + ...
            " model for " + outcome(a))

        if outcome(a) == "BR"

            % Logistic model for banking

            model = fitglm(T, ...
                'y ~ x', ...
                'Distribution', 'Binomial', ...
                'Link', 'logit');

        else

            % Linear model for RT and IT

            model = fitlm(T, 'y ~ x');
        end

        coefficients = model.Coefficients;
        CI = coefCI(model);

        slopeRow = find(strcmp( ...
            coefficients.Properties.RowNames, 'x'));

        beta = coefficients.Estimate(slopeRow);
        SE = coefficients.SE(slopeRow);
        testStatistic = coefficients.tStat(slopeRow);
        modelPValue = coefficients.pValue(slopeRow);

        CILow = CI(slopeRow, 1);
        CIHigh = CI(slopeRow, 2);

        oddsRatio = NaN;
        oddsLow = NaN;
        oddsHigh = NaN;

        if outcome(a) == "BR"

            oddsRatio = exp(beta);
            oddsLow = exp(CILow);
            oddsHigh = exp(CIHigh);
        end

        % Permutation test

        disp("Running " + nPerm + " permutations...")

        permutedBeta = NaN(nPerm, 1);

        for perm = 1:nPerm

            % Shuffle IED values across this participant's trials

            shuffledOrder = randperm(length(xPatient));
            xShuffled = xPatient(shuffledOrder);
            xShuffled = xShuffled(:);

            shuffledTable = table(yPatient, xShuffled, ...
                'VariableNames', {'y', 'x'});

            if outcome(a) == "BR"

                shuffledModel = fitglm(shuffledTable, ...
                    'y ~ x', ...
                    'Distribution', 'Binomial', ...
                    'Link', 'logit');

            else

                shuffledModel = fitlm(shuffledTable, 'y ~ x');
            end

            shuffledCoefficients = shuffledModel.Coefficients;

            shuffledSlopeRow = find(strcmp( ...
                shuffledCoefficients.Properties.RowNames, 'x'));

            permutedBeta(perm) = ...
                shuffledCoefficients.Estimate(shuffledSlopeRow);

            if mod(perm, 100) == 0

                disp("Permutation " + perm + ...
                    " of " + nPerm)

            end
        end

        % Two-sided permutation p-value

        validPermutations = isfinite(permutedBeta);

        permutationPValue = ...
            (sum(abs(permutedBeta(validPermutations)) >= ...
            abs(beta)) + 1) / ...
            (sum(validPermutations) + 1);

        significantPermutation = permutationPValue < 0.05;

        nTrialsUsed = height(T);

        newRow = table(ptID, analysisType(a), outcome(a), ...
            nTrialsUsed, beta, SE, testStatistic, ...
            modelPValue, permutationPValue, ...
            significantPermutation, ...
            CILow, CIHigh, oddsRatio, oddsLow, oddsHigh, ...
            'VariableNames', ...
            {'patientID', 'analysisType', 'outcome', ...
             'nTrials', 'beta', 'SE', ...
             'testStatistic', 'modelPValue', ...
             'permutationPValue', ...
             'significantPermutation', ...
             'CILow', 'CIHigh', ...
             'oddsRatio', ...
             'oddsRatioCILow', ...
             'oddsRatioCIHigh'});

        results = [results; newRow];

        disp("Real coefficient: " + beta)
        disp("Permutation p-value: " + permutationPValue)
    end
end

% Save participant coefficients

writetable(results, ...
    fullfile(outputFolder, ...
    'per_participant_count_and_occurrence_permutation.csv'));

disp(" ")
disp("All participant results:")
disp(results)

disp(" ")
disp("Significant coefficients using permutation p < 0.05:")
disp(results(results.significantPermutation, :))

% Create scatter boxplots

colors = [0.204 0.459 0.702; ...
          0.847 0.333 0.153; ...
          0.250 0.600 0.250];

outcomeNames = ["RT", "IT", "BR"];

rng(1);

figure('Color', 'w', 'Position', [100 100 750 750]);

for panel = 1:2

    subplot(2, 1, panel)
    hold on

    if panel == 1

        wantedAnalysis = "IED count";
        panelTitle = 'IED count';

    else

        wantedAnalysis = "IED occurrence";
        panelTitle = 'IED occurrence';
    end

    for k = 1:3

        rows = results.analysisType == wantedAnalysis & ...
            results.outcome == outcomeNames(k);

        coefficients = results.beta(rows);
        significant = results.significantPermutation(rows);

        validRows = isfinite(coefficients);
        coefficients = coefficients(validRows);
        significant = significant(validRows);

        xBox = k * ones(length(coefficients), 1);

        boxchart(xBox, coefficients, ...
            'BoxFaceColor', colors(k, :), ...
            'BoxWidth', 0.45, ...
            'MarkerStyle', 'none', ...
            'LineWidth', 1.2);

        xScatter = k + ...
            (rand(length(coefficients), 1) - 0.5) * 0.24;

        % Nonsignificant participant coefficients

        scatter(xScatter(~significant), ...
            coefficients(~significant), 25, ...
            colors(k, :), ...
            'filled', ...
            'MarkerFaceAlpha', 0.65, ...
            'MarkerEdgeColor', 'none');

        % Significant participant coefficients

        scatter(xScatter(significant), ...
            coefficients(significant), 38, ...
            colors(k, :), ...
            'filled', ...
            'MarkerEdgeColor', 'k', ...
            'LineWidth', 1);
    end

    yline(0, '--', ...
        'Color', [0.45 0.45 0.45], ...
        'LineWidth', 1);

    xlim([0.5 3.5]);
    xticks(1:3);
    xticklabels({'RT', 'IT', 'BR'});

    ylabel('Participant coefficient', ...
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

exportgraphics(gcf, ...
    fullfile(outputFolder, ...
    'per_participant_count_and_occurrence_permutation.pdf'), ...
    'ContentType', 'vector');

close(gcf)

% Save significant participants based on permutation p-value

significantResults = ...
    results(results.significantPermutation, :);

writetable(significantResults, ...
    fullfile(outputFolder, ...
    'significant_participant_coefficients_permutation.csv'));

% Find the two smallest BR coefficients

BRresults = results(results.outcome == "BR", :);
BRresults = sortrows(BRresults, 'beta', 'ascend');

disp(" ")
disp("Two smallest BR coefficients:")
disp(BRresults(1:min(2, height(BRresults)), :))