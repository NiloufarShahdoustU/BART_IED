function p = chiSquaredTest(vec1, vec2)
    % Counts
    n11 = sum(vec1 == 1); % 1s in vec1
    n10 = sum(vec1 == 0); % 0s in vec1
    n21 = sum(vec2 == 1); % 1s in vec2
    n20 = sum(vec2 == 0); % 0s in vec2
    
    % Construct contingency table
    contingencyTable = [n10 n20; n11 n21];
    
    % Use MATLAB's chi2 test function assuming it's available in your version,
    % or adapt this section according to the statistical function you have.
    [h, p, chi2stat] = fishertest(contingencyTable); % Adapt this line as needed
    

end
