function [] = plotDeltaXc(cyc, xc1, xc1_err, xc3, xc3_err)

errCalc = sqrt(xc1_err.^2 + xc3_err.^2);
deltaXc = abs(xc1 - xc3);

figure
errorbar(cyc, deltaXc, errCalc, 'ow')

xlabel('Applied AC Cycles')
ylabel('\Delta \phi (degrees)')

plot_template()
end

