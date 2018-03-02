function [fitobject ci fig_h chi2] = matlabFit(phi, int, interr, problemVarNames, problemVarValues, start, lower, upper, opt)
% v 9.2 3/1/2018 E R Louden    

% function [fitobject ci fig_h chi2] = matlabFit(phi, int, interr, problemVarNames, problemVarValues, start, lower, upper, opt)
% matlabFit: calls the matlab fitting function, fits several different varieties of gaussians
%   opt = 0: default; standard 3-peak Guassian
%   opt = 1: restricted to a 2-peak Guassian, rewrites i01 and i03 in terms of 
%               itot - the total intensity (itot = i01 + i03)
%                 i1 - the fraction of intensiy in the first peak (i1 = i01 / itot) 
%   opt = 2:  completely unrestrained 3-peak fit
%   opt = 3:  completely unrestrained 2-peak fit
%
%   fitboject       -   smoothed/binned x values, NUM
%   ci              -   fonfidence intervals for the fit
%   fig_h           -   handle for the plot with data / fit
%   chi2            -   calculated chi2 for the fit
%
%   phi             -   dependent data
%   int             -   independent data
%   interr          -   error on the independent data
%   probVarNames    -   parameters that are 'fixed' in a cycle 
%   probVarValues   -   the value to use for 'fixed' parameters
%   start           -   inital guess values
%   lower           -   lower bound for fit parameters
%   upper           -   upper bound for fit parameters
%   opt             -   specifies which fit type to use


%% Set the fit-type

% if opt is not specified, use the default
if (nargin<9)
    opt = 0;
end

% set the generic fit options
fo = fitoptions('Method','NonlinearLeastSquares','StartPoint',start,'Lower',lower,'Upper',upper);

if(opt == 0)
    % default fit type

    % ORDER
    % fwhm, i01,  xc1, i02, xc2, i03, xc3, y0
    ft = fittype('y0 + (i01/(fwhm * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(x-xc1)^2/(fwhm^2/log(4))) + (i02/(fwhm * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(x-xc2)^2/(fwhm^2/log(4))) + (i03/(fwhm * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(x-xc3)^2/(fwhm^2/log(4)))',...
                 'options',fo,'problem',problemVarNames);    

elseif(opt == 1)
   % fit type to better determine total intensity
   % NOTE: only for two-peak Gaussian fits

    % ORDER
    % i1, itot,  xc1, xc3, y0
   ft = fittype('y0 + itot*((i1/(fwhm * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(x-xc1)^2/(fwhm^2/log(4))) + ((1-i1)/(fwhm * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(x-xc3)^2/(fwhm^2/log(4))))',...
                'options',fo,'problem',problemVarNames);    

elseif(opt == 2)
    %open fwhm

    % ORDER
    % fwhm1, i01, xc1, fwhm2, i02, xc2, fwhm3, i03, xc3, y0
    ft = fittype('y0 + (i01/(fwhm1 * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(x-xc1)^2/(fwhm1^2/log(4))) + (i02/(fwhm2 * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(x-xc2)^2/(fwhm2^2/log(4))) + (i03/(fwhm3 * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(x-xc3)^2/(fwhm3^2/log(4)))',...
                 'options',fo,'problem',problemVarNames);   

elseif(opt == 3)
   % open fwhm
   % NOTE: only for two-peak Gaussian fits

   % ORDER
   % fwhm1, itot, i1 fwhm3, xc1, xc3, y0
   ft = fittype('y0 + itot*((i1/(fwhm1 * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(x-xc1)^2/(fwhm1^2/log(4))) + ((1-i1)/(fwhm3 * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(x-xc3)^2/(fwhm3^2/log(4))))',...
                'options',fo,'problem',problemVarNames);    

else
    disp('error')
end


%% Perform the fit

[fitobject, gof] = fit(phi', int', ft, 'problem', problemVarValues)
ci = confint(fitobject);

 
%% Plot and format the fit
 
fig_h = plotData(phi, int, interr, '\phi - \phi_0', 'Intensity (arb. units)', 'Fit Cycle')
hold on
plot(fitobject)
set(findobj('Type','Line','Color','r'),'Color',[0.08, 0.17, 0.55])
xlabel('\phi - \phi_0')
ylabel('Intensity (arb. units)')
plot_template(2)
hold off


%% Calculate the chi^2

if(opt == 0)
    y =fitobject.y0 + (fitobject.i01/(fitobject.fwhm * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(phi-fitobject.xc1).^2/(fitobject.fwhm^2/log(4))) + (fitobject.i02/(fitobject.fwhm * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(phi-fitobject.xc2).^2/(fitobject.fwhm^2/log(4))) + (fitobject.i03/(fitobject.fwhm * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(phi-fitobject.xc3).^2/(fitobject.fwhm^2/log(4)));

elseif(opt == 1)
    y = fitobject.y0 + fitobject.itot*((fitobject.i1/(fitobject.fwhm * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(phi-fitobject.xc1).^2/(fitobject.fwhm^2/log(4))) + ((1-fitobject.i1)/(fitobject.fwhm * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(phi-fitobject.xc3).^2/(fitobject.fwhm^2/log(4))));

elseif(opt == 2)
    y = fitobject.y0 + (fitobject.i01/(fitobject.fwhm1 * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(phi-fitobject.xc1).^2/(fitobject.fwhm1^2/log(4))) + (fitobject.i02/(fitobject.fwhm2 * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(phi-fitobject.xc2).^2/(fitobject.fwhm2^2/log(4))) + (fitobject.i03/(fitobject.fwhm3 * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(phi-fitobject.xc3).^2/(fitobject.fwhm3^2/log(4)));

elseif(opt == 3)
    y = fitobject.y0 + fitobject.itot*((fitobject.i1/(fitobject.fwhm1 * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(phi-fitobject.xc1).^2/(fitobject.fwhm1^2/log(4))) + ((1-fitobject.i1)/(fitobject.fwhm3 * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(phi-fitobject.xc3).^2/(fitobject.fwhm3^2/log(4))));

else
    disp('error')

end

chi2_l = (y - int).^2 ./ int;
chi2 = sum(chi2_l);

end