function [h] = CM(phi, Int, cm_x)
global mf_fitter
global status_flags
   
 if(get(mf_fitter.handles.fit_dropdown,'Value') == 3 )
%     offset = status_flags.analysis_modules.sectors.theta;
     xc2 = mf_fitter.fit_data.center2(:,1);
 end
% else
%     offset = 0;
% end 
    % Code for smoothing out the data fluctuations for the 
    % initialize
        phi = mf_fitter.SmoothedData.phi;
        intensity = mf_fitter.SmoothedData.Int';
        xc1 = mf_fitter.fit_data.center1(:,1);% - offset; %mf_fitter.center_plot.x1;
        %xc2 = mf_fitter.fit_data.center2; % mf_fitter.center_plot.x2;
        xc3 = mf_fitter.fit_data.center3(:,1);% - offset; %mf_fitter.center_plot.x3;
        %xgs = mf_fitter.center_plot.xgs;
        %xms = mf_fitter.center_plot.xms;
        CYC = mf_fitter.fit_data.cycles;
        cycles = CYC
        %climits = [mf_fitter.climits(1), mf_fitter.climits(2)]
     
        num_angles = length(phi);
        cyc = length(cycles);
        bg = 0;
        xcycles = cycles;
        intensity(num_angles+1,:) = 0;  % gives junk column/row not plotted by pcolor
        intensity(:,cyc+1) = 0;
        backsub_int = zeros(size(intensity));  
        integrated_intensity = zeros(size(intensity));
        %smoothbox_I = zeros(size(intensity));
        %smoothgaus_I = zeros(size(intensity));
        

    % offset phi to account for the way pcolor positions pixels
        global phi_shift
        phi_shift = phi;
        for i=1:(length(phi)-1)
            delta_phi(i) = (phi(i+1) - phi(i))/2;
            phi_shift(i) = phi(i) - delta_phi(i);
        end  

    % perform the background subtraction and intensity integration (=1)    
        for i=1:cyc
            % Find the average background 
                bg =  sum(intensity(1:5,i) + intensity((num_angles-4):num_angles,i)) / 8;
                backsub_int(:,i) = intensity(:,i) - bg;
   
            % Make the integrated intensity equal to one across every column
                s = sum(backsub_int(:,i));
                integrated_intensity(:,i) = backsub_int(:,i)/s;
        end
        
%     % Smoothing
% 
%     % Box Smoothing
%     % Take the mean of the pixel and it's two nearest neighbors
%         for i = 1:cyc  
%             for j=1:num_angles
%                 if(j==1) smoothbox_I(j,i) = mean(integrated_intensity((j):(j+1),i));
%                 elseif (j==num_angles) smoothbox_I(j,i) = mean(integrated_intensity((j-1):(j),i));
%                 else smoothbox_I(j,i) = mean(integrated_intensity((j-1):(j+1),i));
%                 end
%             end
%         end
% 
%     % Gaussian
%     % Take the mean again, this time weighted by the gaussian with specified fwhm 
%         fwhm = 0.5;
%         for i = 1:cyc
%             total_I = [];
%             sum_weights = [];
%             for j=1:num_angles
%                 for pix=1:num_angles
%                     total_I(pix) = exp((-(j-pix)^2)/fwhm)*integrated_intensity(pix,i);
%                     sum_weights(pix) = exp((-(j-pix)^2)/fwhm);
%                 end
%                 smoothgaus_I(j,i) = sum(total_I)/sum(sum_weights);
%             end
%         end

    % For Logarithmic intensity scale    
        %log_i = log(smoothgaus_I-min(min(smoothgaus_I)));
        log_i = log(integrated_intensity - min(min(integrated_intensity)));
        
    % Add junk cycle and phi value (based on the way matlab command pcolor works)
        cycles(cyc+1) = 2*cycles(cyc);
        xcycles(length(xcycles)+1) = 2*xcycles(length(xcycles));
        %xms(length(xms)+1) = 2*xms(length(xms));
        phi_shift(num_angles+1) = phi_shift(num_angles) + delta_phi(num_angles-1);
        phi(num_angles+1) = phi(num_angles) + delta_phi(num_angles-1);
        
    % Define values for asprep column    
        asprep_log_i =  log_i(:,1:2);
        asprep_cycles = [0,1];
%         
         for i=1:(length(xcycles)-1)
             xcycles(i) = sqrt(xcycles(i+1)*xcycles(i));
         end
         
         xcycles = xcycles(1:(length(xcycles)-1))
%         end
%         for i=1:(length(xms)-1)
%             xmscycles(i) = sqrt(xms(i+1)*xms(i));
%         end
%         xgscycles(length(xgs)) = sqrt(xgs(length(xgs))*2*xgs(length(xgs)));
%         xmscycles(length(xms)) = sqrt(xms(length(xms))*2*xms(length(xms)));
%         
    % Figure Formatting (gives the correct aspect ratio)
        f = 2;
        fig_size = f*[5.7 4.45];
        offset = f*[0.1 0.1];
        fig_start = [0 0];
        graph_start = f*[1.25 1.25];
        graph_sep = 0.15;
        asprep_size = 0.3;

        position1 = [graph_start asprep_size fig_size(2)-graph_start(2)];
        position2 = [graph_start(1)+asprep_size+graph_sep graph_start(2) fig_size(1)-(asprep_size+graph_sep+graph_start(1)) fig_size(2)-graph_start(2)];

    % Make the figure window    
        figure
        set(gcf,'Units','centimeters','Position',[fig_start fig_size+offset]);

    % Make the as prepped column
        g = subplot(1,2,1);
            g = pcolor(asprep_cycles,phi_shift,asprep_log_i); 
            %s = sprintf('Azimuthal Angle (%c)', char(176));
           % ylabel(s);
            ylabel('\phi - \phi_0 (degrees)')
      
    % Make the rest of the colormap figure        
        h = subplot(1,2,2);
            h = pcolor(cycles,phi_shift,log_i); 
            %str = {'Logarithmic Integrated Intensity', ['\fontsize{10}', sprintf(' with %d fwhm Gaussian Smoothing', fwhm)], '\fontsize{10} 1.29 Angle Binning'}
            xlabel('Applied || AC Cycles');
            hold on
            % Add the fitted data centers
                %plot(xgscycles,xc1,'.','MarkerSize',10,'MarkerEdgeColor',[1 1 1]);
                %plot(xmscycles,xc2,'.','MarkerSize',10,'MarkerEdgeColor',[1 1 1]);
                %plot(xgscycles,xc3,'.','MarkerSize',10,'MarkerEdgeColor',[1 1 1]);
                plot(xcycles',xc1,'o','MarkerSize',4,'MarkerEdgeColor',[1 1 1]);
                plot(xcycles',xc3,'o','MarkerSize',4,'MarkerEdgeColor',[1 1 1]);
           
                if(get(mf_fitter.handles.fit_dropdown,'Value') == 3 )
                    plot(xcycles',xc2,'o','MarkerSize',4,'MarkerEdgeColor',[1 1 1]);  
                end
           hold off
  
    % Format the figures
        CM_lim = get(gca,'CLim')
        sort_int = sort(log_i(:))
        l = ceil(length(sort_int(:))/10)
        CM_lim(1) = mean(sort_int(2:l)) % always one value that is - inf, want to ignore this one
        %CM_lim = [min(min())]
        CM_lim(2) = max(max(log_i))
        Climits = [CM_lim(1) CM_lim(2)]
    
        subplot(1,2,1); 
        set(gca,...
            'Color',[1 1 1],...
            'FontName','Arial','FontSize',f*8,...
            'Xcolor',[0 0 0],'Ycolor',[0 0 0],...
            'Clim',Climits,...
            'XTickLabelMode','Manual','XTickLabel',{' ','0',' '},...
            'Ylim',[-22 22],...
            'Units','centimeters','Position', position1);
            ax.Box = 'on';
            
        h_ylabel = get(gca,'YLabel');
        set(h_ylabel,'FontName','Arial','FontSize',f*8,'Color',[0 0 0]);
 
        subplot(1,2,2);
        set(gca,...
            'Color',[1 1 1],...
            'FontName','Arial','FontSize',f*8,...
            'Xcolor',[0 0 0],'Ycolor',[0 0 0],...
            'YTickLabelMode','Manual','YTickLabel',{},...
            'XTickLabelMode','Auto','XTick',[10,100,1000,10000],...
            'XScale','Log','Clim',Climits,...
            'Ylim',[-22 22],...
            'Units','centimeters','Position', position2);
         ax.Box = 'on';
         

        %set(gca,'XScale','log','CLim',CM_lim)
%        
       % plot_template
       
             
        h_xlabel = get(gca,'XLabel');
        set(h_xlabel,'FontName','Arial','FontSize',f*8,'Color',[0 0 0]);
        
        titlename = [mf_fitter.folder ' Colormap - Smoothing FWHM ' num2str(mf_fitter.smoothing.fwhm) ' Step Size ' num2str(mf_fitter.smoothing.step)];
        title(titlename,'FontSize',f*12,'Fontname','Arial','Color','black');
       
        h = gcf;
        mf_fitter_save('save', 'Colormap', h)

%   indicesToDelete = (find(abs(phi) > 31)); % change to reference phi ring later
%   sizeInt = size(Int)
%   for i = 1:length(indicesToDelete)
%       for j = 1:sizeInt(1)
%       Int(j,indicesToDelete(i)) = min(Int(j,:))
%       end
%   end
%     
%    % data already smoothed and stored in 2D array above
%    % just swap so that x corresponds to rows
%    Int = Int';
%    
%    % offset phi to account for the way pcolor positions pixels
%     %global phi_shift
% %     phi_shift = phi;
% %     for i=1:(length(phi)-1)
% %         delta_phi(i) = (phi(i+1) - phi(i))/2
% %         phi_shift(i) = phi(i) - delta_phi(i)
% %     end  
% %  
%     
%     % perform the background subtraction and intensity integration (=1) 
%     ang = numel(phi);    
%     N = length(mf_fitter.fit_data.names);
%     
%     for i=1:N
%         % Find the average background 
%             %bg =  sum(Int(1:3,i) + Int((ang-2):ang,i)) / 6;
%             %backsub_int(:,i) = Int(:,i) - bg;
%             
%             %new BG sub method
%             % grab the (25% of total) minimum values and average them for BG
%             temp_int = sort(Int(:,i));
%             BG = temp_int(1:0.25*numel(temp_int));
%             BG = sum(BG)/numel(BG);
%             
%             backsub_int(:,i) = Int(:,i) - BG;
% 
%         % Make the integrated intensity equal to one across every column
%             %s = sum(backsub_int(:,i));
%             %integrated_int(:,i) = backsub_int(:,i)/s;
%             % check
%             %sum(integrated_int(:,i))
%         
%         % Normalize to the maximum in each column
%             m = max(backsub_int(:,i));
%             norm_int(:,i) = backsub_int(:,i)/m;
%             
%     end
% 
%     
%     % For Logarithmic intensity scale    
%         %log_i = log(integrated_int-min(min(integrated_int)));
%         log_i = log(norm_int-min(min(norm_int)));
%         %log_i = log(norm_int)
%         
%         
% %     % Add junk x and phi value (based on the way matlab command pcolor works)
% %         cm_x(N+1) = 2*cm_x(N);
% %         phi_shift(ang+1) = phi_shift(ang) + delta_phi(ang-1);
% %         phi(ang+1) = phi(ang) + delta_phi(ang-1);
% %         
%         
%    % plot the colormap
%    intName = {'Int','Smoothed Intensity','log_i','Logarithmic Intensity Scale','norm_int','Normalized to Max Per Field'}
%    for i = 2:2
%    %for i = 3
%        h = figure
%        
%        % get intensity data
%        cData = eval(char(intName(2*i-1)));
%        
%        % Add junk Int column, row
%        %size(cData)
%        %ang
%        %length(phi)
%        %length(cm_x)
%        cData(2:(ang+1),:) = cData;
%        %cData(1,:) = 0; % shows up as red
%        %cData(ang,:) = 0; % shows up as red
%        %cData(ang+1,:) = 0;
%        %cData(:,2) = 0; % shows up as red
%        %cData(:,N+1) = 0;
%        %cData(:,N) = 0; % shows up as red
%        
%        % make the colormap and label axes
%        phi 
%        %phi_shift
%        %pcolor(cm_x,phi,cData)
%        figure
%        imagesc(cData)
%        title(intName(2*i))
%        xlabel('H(T)')
%        ylabel('\phi  - \phi_0  (degrees)')
%        
%        figure
%        size(cm_x)
%        size(phi)
%        size(cData)
%        imagesc(cm_x,phi,cData)
%        
%        CM_lim = get(gca,'CLim')
%        sort_int = sort(cData(:))
%        l = ceil(length(sort_int(:))/10)
%        CM_lim(1) = mean(sort_int(2:l)) % always one value that is - inf, want to ignore this one
%        %CM_lim = [min(min())]
%        CM_lim(2) = max(max(cData))
%        set(gca,'XScale','log','CLim',CM_lim)
%        
%        plot_template
       
%   end
   
end