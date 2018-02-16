function [fitobject ci fig_h] = matlabFit(phi, int, interr, problemVarNames, problemVarValues, start, lower, upper, opt)
% Code to call the matlab fitting function
% Offers two different fit types
%   opt = 0: default; standard 3-peak Guassian
%   opt = 1: restricted to a 2-peak Guassian, rewrites i01 and i03 in terms of 
%               itot - the total intensity (itot = i01 + i03)
%                 i1 - the fraction of intensiy in the first peak (i1 = i01 / itot) 
%   opt = 2:  completely unrestrained

% v 9.1
% 5/13/2015 MFF 
% last updated 2/15/2018 - E R Louden

        % if opt is not specified, use the default
        if (nargin<9)
            opt = 0;
        end

        fo = fitoptions('Method','NonlinearLeastSquares','StartPoint',start,'Lower',lower,'Upper',upper);

        if(opt == 0)
            % defualt fit type
            
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
        
        
         [fitobject, gof] = fit(phi', int', ft, 'problem', problemVarValues)
         ci = confint(fitobject);
         fig_h = plotData(phi, int, interr, '\phi - \phi_0', 'Intensity (arb. units)', 'Fit Cycle')
         hold on
         plot(fitobject)
         set(findobj('Type','Line','Color','r'),'Color',[0.08, 0.17, 0.55])
         xlabel('\phi - \phi_0')
         ylabel('Intensity (arb. units)')
         plot_template(1)
         hold off
   
end