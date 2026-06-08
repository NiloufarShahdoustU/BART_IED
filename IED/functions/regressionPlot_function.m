function regressionPlot_function(x1, x2)
    % Ensure inputs are column vectors
    x1 = x1(:);
    x2 = x2(:);
    
    % Create a table for the regression model
    tbl = table(x2, x1, 'VariableNames', {'x2', 'x1'});
    
    % Fit the linear regression model
    mdl = fitlm(tbl, 'x1 ~ x2');
    
    % Plot the regression
    h = plot(mdl);
    h(3).Visible = 'off'; % Turn off the lower confidence bound
    h(4).Visible = 'off'; % Turn off the upper confidence bound
    
    % Remove confidence bounds from the legend
    legend('off'); % Update legend to include only data and fit lines
    
    % Adjust plot appearance
    h(1).Marker = 'o'; % Set marker to circle
    h(1).MarkerFaceColor = 'k'; % Fill marker with black color
    h(1).MarkerEdgeColor = 'k'; % Set marker edge color to black
    h(1).LineStyle = 'none'; % Remove lines connecting markers
    
    h(2).LineWidth = 2; % Set regression line to be thicker
    
    xlim([0 max(x2)]+1);
    ylim([0 max(x1)]+1);
    xlabel(' ');
    ylabel(' ');
    title(' ');
    axis square;
    set(gca, 'box', 'off');
    set(gca, 'TickDir', 'out');
end
