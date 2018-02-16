function [h] = plotData(x,y,z, xName, yName, titleName)
    h = figure;
    errorbar(x,y,z,'o','Color',[0.64, 0.08, 0.16],'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0.64, 0.08, 0.16]);
    
    xlabel(xName);
    ylabel(yName);
    title(titleName);
    
    plot_template(1)
end