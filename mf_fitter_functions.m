%% General Functions
% [phi, Int, Int_err] = getData(img_num)
% [h] = plotData(x,y,z, xName, yName, titleName)
% [] = plot_template()
% [] = plotXC(cm_han, cyc, xc1, xc1_err, xc3, xc3_err)
% [h] = CM(phi, Int, cm_x)
% [fitobject ci fig_h] = matlabFit(phi, int, interr, problemVarNames, problemVarValues, start, lower, upper)
% [phiS, intS, intErrS] = smooth( nphi, nInt, nIntErr, fn, phisym, symtrz, phirng, gwidth)
 
 disp('test')

function [phi, Int, Int_err] = getData(img_num)
        global status_flags
        global mf_fitter
        global plot_info
        global grasp_handles
    
        curAB = get(grasp_handles.window_modules.radial_average.azimuth_bin,'String');
       
        % load desired depth file and update main grasp GUI
        % the +1 is necessary as file 1 represents a sum
        status_flags.selector.fd = img_num+1;
        set(grasp_handles.window_modules.radial_average.azimuth_bin,'String',0.1);
        radial_average_callbacks2('azimuth_bin');
        grasp_update;

        % plot I vs xi and create the current figure, store handle
        radial_average_callbacks('averaging_control','azimuthal');
        mf_fitter.handles.plot_handle = gcf;
        close(gcf)
        grasp_update;
            
        % save data into appropriate variables
        phi = plot_info.export_data(:,1);
        Int = plot_info.export_data(:,2);
        Int_err = plot_info.export_data(:,3);
        
        % return to original selected angle binning
        set(grasp_handles.window_modules.radial_average.azimuth_bin,'String',curAB)
        radial_average_callbacks2('azimuth_bin');
        grasp_update;

end

function [h] = plotData(x,y,z, xName, yName, titleName)
    h = figure
    errorbar(x,y,z,'ok')
    
    xlabel(xName)
    ylabel(yName)
    title(titleName)
    
    plot_template
end

function [] = plot_template()
    scale = 1.5;
    fig_size = [7.3 5.7];
    fig_position = [1.25 1.25];
    position = scale*[fig_position fig_size];

    set(gca,...
        'Color',[1 1 1],...
        'FontName','Arial','FontSize',scale*8,...
        'Xcolor',[0 0 0],'Ycolor',[0 0 0],...
        'XMinorTick','on','YMinorTick','on',...
        'Box','on');
        %'XLim',xlim,'XMinorTick', 'on',...
        %'YLim',ylim,'YMinorTick', 'on',...
        %'Units','centimeters','Position', position,...

hold off

end

function [fitobject ci fig_h] = matlabFit(phi, int, interr, problemVarNames, problemVarValues, start, lower, upper)
        % ORDER
        % y0, fwhm, i01,  xc1, i02, xc2, i03, xc3

% 
         fo = fitoptions('Method','NonlinearLeastSquares','StartPoint',start,'Lower',lower,'Upper',upper);
% 
         ft = fittype('y0 + (i01/(fwhm * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(x-xc1)^2/(fwhm^2/log(4))) + (i02/(fwhm * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(x-xc2)^2/(fwhm^2/log(4))) + (i03/(fwhm * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(x-xc3)^2/(fwhm^2/log(4)))',...
                         'options',fo,'problem',problemVarNames);        

%        
         [fitobject, gof] = fit(phi', int', ft, 'problem', problemVarValues)
         ci = confint(fitobject);
         fig_h = plotData(phi, int, interr, '\phi - \phi_0', 'Intensity (arb. units)', 'Fit Cycle 1')
         hold on
         plot(fitobject)
         hold off
%   
end

function [] = plotXC(cm_han, cyc, xc1, xc1_err, xc3, xc3_err)
figure(cm_han)
hold on

errorbar(cyc, xc1, xc1_err, 'ow')
errorbar(cyc, xc3, xc3_err, 'ow')
end

function [h] = CM(phi, Int, cm_x)
global mf_fitter

  indicesToDelete = (find(abs(phi) > 31)); % change to reference phi ring later
  sizeInt = size(Int)
  for i = 1:length(indicesToDelete)
      for j = 1:sizeInt(1)
      Int(j,indicesToDelete(i)) = min(Int(j,:))
      end
  end
    
   % data already smoothed and stored in 2D array above
   % just swap so that x corresponds to rows
   Int = Int';
   
   % offset phi to account for the way pcolor positions pixels
    %global phi_shift
    phi_shift = phi;
    for i=1:(length(phi)-1)
        delta_phi(i) = (phi(i+1) - phi(i))/2;
        phi_shift(i) = phi(i) - delta_phi(i);
    end  
    
    % perform the background subtraction and intensity integration (=1) 
    ang = numel(phi);
    
    N = length(mf_fitter.fit_data.names);
    
    for i=1:N
        % Find the average background 
            %bg =  sum(Int(1:3,i) + Int((ang-2):ang,i)) / 6;
            %backsub_int(:,i) = Int(:,i) - bg;
            
            %new BG sub method
            % grab the (25% of total) minimum values and average them for BG
            temp_int = sort(Int(:,i));
            BG = temp_int(1:0.25*numel(temp_int));
            BG = sum(BG)/numel(BG);
            
            backsub_int(:,i) = Int(:,i) - BG;

        % Make the integrated intensity equal to one across every column
            %s = sum(backsub_int(:,i));
            %integrated_int(:,i) = backsub_int(:,i)/s;
            % check
            %sum(integrated_int(:,i))
        
        % Normalize to the maximum in each column
            m = max(backsub_int(:,i));
            norm_int(:,i) = backsub_int(:,i)/m;
            
    end

    
    % For Logarithmic intensity scale    
        %log_i = log(integrated_int-min(min(integrated_int)));
        log_i = log(norm_int-min(min(norm_int)));
        %log_i = log(norm_int)
        
        
    % Add junk x and phi value (based on the way matlab command pcolor works)
        cm_x(N+1) = 2*cm_x(N);
        phi_shift(ang+1) = phi_shift(ang) + delta_phi(ang-1);
        phi(ang+1) = phi(ang) + delta_phi(ang-1);
        
        
   % plot the colormap
   intName = {'Int','Smoothed Intensity','log_i','Logarithmic Intensity Scale','norm_int','Normalized to Max Per Field'}
   for i = 2:2
   %for i = 3
       h = figure
       
       % get intensity data
       cData = eval(char(intName(2*i-1)));
       
       % Add junk Int column, row
       %size(cData)
       %ang
       %length(phi)
       %length(cm_x)
       cData(2:(ang+1),:) = cData;
       %cData(1,:) = 0; % shows up as red
       %cData(ang,:) = 0; % shows up as red
       cData(ang+1,:) = 0;
       %cData(:,2) = 0; % shows up as red
       cData(:,N+1) = 0;
       %cData(:,N) = 0; % shows up as red
       
       % make the colormap and label axes
       pcolor(cm_x,phi_shift,cData)
       title(intName(2*i))
       xlabel('H(T)')
       ylabel('\phi  - \phi_0  (degrees)')
       
       CM_lim = get(gca,'CLim')
       sort_int = sort(cData(:))
       l = ceil(length(sort_int(:))/10)
       CM_lim(1) = mean(sort_int(2:l)) % always one value that is - inf, want to ignore this one
       %CM_lim = [min(min())]
       CM_lim(2) = max(max(cData))
       set(gca,'XScale','log','CLim',CM_lim)
       
       plot_template
       
   end
   
end


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
