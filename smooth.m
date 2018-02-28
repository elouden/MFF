function [phiS, intS, intErrS] = smooth( nphi, nInt, nIntErr, fn, phisym, symtrz, phirng, gwidth) 
% v 9.2 2/21/2018 E R Louden    

% function [phiS, intS, intErrS] = smooth( nphi, nInt, nIntErr, fn, phisym, symtrz, phirng, gwidth) 
% smooth:  performs binning and smoothing on raw neutron counts
%           exported from GRASP with an angle binning of 0.1 (see function getData
%           for the processof exporting).  
%           The original verion of this function was drafted by M R Eskildsen.
%           Data can either by a single vector of ints, or a matrix of int vs some control parameter
%   
%   phiS        -    smoothed/binned x values, NUM
%   intS        -    smoothed/binned y values, NUM
%   ntErrS      -    error on y, NUM
%
%   nphi        -    initial x values, ideally exported with AB = 0.1, NUM 
%   nInt        -    initial y values, NUM
%   nIntErr     -    initial error on y, NUM
%   fn          -    the number of passed y data sets if giving a matrix of   y, otherwise 1, INT
%   physym      -    the value to center data at, NUM
%   symtrz      -    if data should be symmetrized, input 1 to reflect data about physym, INT     
%   phirng      -    the final x values to plot against, [NUM (start), NUM (end), NUM (step size)] 
%   gwidth      -    the guassian fwhm determining smoothing amount, NUM


N = numel(fn);
phi = phirng(1):phirng(3):phirng(2);
Int = zeros([N length(phi)]);
IntErr = Int;
      
    for n = 1:N
        
        % Center data about phisym

        % GRASP will output data above 360 as being at 1 - check for this case
        % if it did, the sorted phi will not match input nphi
        nphiSORT = sort(nphi);
        if(nphiSORT ~= nphi)
            nphi(find(nphi < nphi(1))) = nphi(find(nphi < nphi(1))) + 360;
        end

        nphi = nphi - phisym;

        
        % Symmetrize
        if(symtrz);
            nphi = [nphi-phisym(n); -(nphi-phisym(n))];
            nInt = [nInt; nInt];
            nIntErr = [nIntErr; nIntErr];
        end

        
         % Smooth

         WhtInt = Int(1,:); % weighted intensity
         WhtIntErrSq = WhtInt; % weighted error
         Wht = WhtInt;  % total weight for normalizing

         % loop over the desired bin locations
         for i = 1:length(phi)

             % for each bin, sum up the intensity based on weighted average for distance from nominal bin position 
             for j = 1:length(nphi)
                 WhtInt(j) = exp(-(phi(i)-nphi(j))^2/gwidth^2)*nInt(j);
                 WhtIntErrSq(j) = (exp(-(phi(i)-nphi(j))^2/gwidth^2)*nIntErr(j))^2;
                 Wht(j) = exp(-(phi(i)-nphi(j))^2/gwidth^2);
             end

             % The final intensity at this phi is the sum of the weighted intensity divided by the total weight
             Int(n,i) = sum(WhtInt)/sum(Wht);
             IntErr(n,i) = sqrt(sum(WhtIntErrSq))/sum(Wht);
         end

         % once smoothing is done, store the final variables into the output variables
         phiS = phi;
         intS = Int;
         intErrS = IntErr;

    end
    
end