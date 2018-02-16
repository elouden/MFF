function [] = plotXC(cm_han, cyc, xc1, xc1_err, xc3, xc3_err)
figure(cm_han)
hold on

errorbar(cyc, xc1, xc1_err, 'ow')
errorbar(cyc, xc3, xc3_err, 'ow')
end