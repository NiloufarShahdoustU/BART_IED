function annotate_significance(x1, x2, p_val, y, threshold)
    line([x1 x2], [y y], 'Color', 'k', 'LineWidth', 1.5);
    if p_val < threshold
        text(mean([x1 x2]), y * 1.05, '*', 'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
    else
        text(mean([x1 x2]), y * 1.05, 'ns', 'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
    end
end

