%% Sort regions by decreasing number of participants

regionSummary = sortrows( ...
    regionSummary, ...
    {'TotalParticipantsWithElectrodes', 'TotalIEDs'}, ...
    {'descend', 'descend'});

%% Calculate mean and standard deviation of nTrials

validTrialRows = ...
    isfinite(participantTrials.NTrials) & ...
    participantTrials.NTrials > 0;

validNTrials = participantTrials.NTrials(validTrialRows);

if isempty(validNTrials)

    meanNTrials = NaN;
    stdNTrials = NaN;
    minimumNTrials = NaN;
    maximumNTrials = NaN;
    nParticipantsWithTrialData = 0;

else

    meanNTrials = mean(validNTrials);
    stdNTrials = std(validNTrials, 0);
    minimumNTrials = min(validNTrials);
    maximumNTrials = max(validNTrials);
    nParticipantsWithTrialData = length(validNTrials);

end

trialSummary = table( ...
    nParticipantsWithTrialData, ...
    meanNTrials, ...
    stdNTrials, ...
    minimumNTrials, ...
    maximumNTrials, ...
    'VariableNames', { ...
        'NParticipants', ...
        'MeanNTrials', ...
        'SDNTrials', ...
        'MinimumNTrials', ...
        'MaximumNTrials'});

%% Display results

fprintf('\n============================================================\n');
fprintf('REGIONAL SUMMARY\n');
fprintf('============================================================\n');

disp(regionSummary);

fprintf('\n============================================================\n');
fprintf('NUMBER OF TRIALS\n');
fprintf('============================================================\n');

disp(trialSummary);

fprintf('nTrials = %.2f +/- %.2f across participants\n', ...
    meanNTrials, stdNTrials);

%% Save only the two requested CSV files

regionalSummaryFile = fullfile( ...
    outputFolderName, ...
    'regional_summary_table.csv');

nTrialsSummaryFile = fullfile( ...
    outputFolderName, ...
    'nTrials_summary.csv');

writetable(regionSummary, regionalSummaryFile);
writetable(trialSummary, nTrialsSummaryFile);

fprintf('\nSaved:\n');
fprintf('%s\n', regionalSummaryFile);
fprintf('%s\n', nTrialsSummaryFile);