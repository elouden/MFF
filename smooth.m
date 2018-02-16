function [phiS, intS, intErrS] = smooth( nphi, nInt, nIntErr, fn, phisym, symtrz, phirng, gwidth) %, intrng, savfig, clsfig, hl, fig_h, cm, cm_x, path)   
%% Allocating memory for processed data
N = numel(fn);
phi = phirng(1):phirng(3):phirng(2);
Int = zeros([N length(phi)]);
IntErr = Int;
      
for n = 1:N
% center phi range around zero (subtract 360 from values > 180)
%     data = readtable([pathname char(fn(n))],'HeaderLines',hl(n),'Delimiter','\t'); 
%     nphi = data{:,{'Az_Angle'}};
%     nInt = data{:,{'I'}};
%     nIntErr = data{:,{'Err_I'}};
%     
%     %center about I-phase position (stored in phisym)

%      nphi = nphi - phisym
%       
     phiabove = find(nphi > 180);
%     phisym
     if(phisym(n) > 180)
         j = 1;
        if(length(phisym) == 1) pnum = 1; else pnum = n;   end
         nphi = nphi - phisym(pnum)
     else
         nphi(phiabove) = nphi(phiabove) - 360;
%     
    %center about I-phase position (stored in phisym)
     j = 0;
% %     if(phisym > 180)
% %         j = 1;
% %     end
         if(length(phisym) == 1) pnum = 1; else pnum = n;   end
%         phisym(pnum)
         max(nphi)
         min(nphi)
         nphi = nphi - phisym(pnum) + j*360;
     end
%     
    % Symmetrize
    if(symtrz);
        nphi = [nphi-phisym(n); -(nphi-phisym(n))];
        nInt = [nInt; nInt];
        nIntErr = [nIntErr; nIntErr];
    end
    
     % Smooth
     WhtInt = Int(1,:);
     WhtIntErrSq = WhtInt;
     Wht = WhtInt;
%     if(phisym > 180)
%         disp('break')
%     end
     for i = 1:length(phi)
         for j = 1:length(nphi)
             WhtInt(j) = exp(-(phi(i)-nphi(j))^2/gwidth^2)*nInt(j);
             WhtIntErrSq(j) = (exp(-(phi(i)-nphi(j))^2/gwidth^2)*nIntErr(j))^2;
             Wht(j) = exp(-(phi(i)-nphi(j))^2/gwidth^2);
         end
         length(fn)
         phisym
         if(phisym(n) > 180 & length(fn)>1)
             disp('break')
         end
         Int(n,i) = sum(WhtInt)/sum(Wht);
         IntErr(n,i) = sqrt(sum(WhtIntErrSq))/sum(Wht);
     end
%     
       phiS = phi;
       intS = Int;
       intErrS = IntErr;
%     %Store in data structure
%     name = ['Field_' char(fh(n))];
%     ILL_Dat.(name).phi = phi;
%     ILL_Dat.(name).Int = Int(n,:);
%     ILL_Dat.(name).IntErr = IntErr(n,:);
%     
end
end