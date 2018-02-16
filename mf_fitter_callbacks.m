function [ output_args ] = mf_fitter_callbacks( to_do, num_peaks, img_num, variable)
%MF_FITTER_CALLBACKS VERSION 8.0
%   to_do:  STRING, switches between the various callback tasks with cases
%                initialize, set_params, set_(free, fwhm, centers, both), fit, data storage, discriminate, refit, fraction, update
%   num_peaks:  INTEGER, sets the number of gaussian peaks for fitting
%   img_num:  INTEGER, selects the index of the file to manipulate
%   variable:  STRING, used to average various data arrays

% v 9.1
% 5/8/2017 MFF Liz

% Global Variables
global mf_fitter;
global status_flags;
global grasp_handles;
global plot_info;


% To_Do Switch
switch to_do
    
    case 'initialize'
        % case run upon first hitting 'Go!'
        % clears out any remaining data from previous runs
        
        % determine number of files ('depth') for use elsewhere
        mf_fitter.depth = status_flags.selector.fdpth_max - 1;
        
        % create inner and outer peaks vector
        % initizlize it to the 3 peak option
        for j = 1:mf_fitter.depth
            mf_fitter.inner_peak(j,1) = 1;
            mf_fitter.outer_peaks(j,1) = 1;
        end
        
        % initialize fit_data, and curve fit window parameters to zero
        var_fit_names = {'background','fwhm','center1','center2','center3','intensity1','intensity2','intensity3','chi2','fms','fgs'};
        var_data_names = {'intensity','intensity_err','angle'};
        var_cplot_names = {'x1','x1_err','x2','x2_err','x3','x3_err'};

        for i=1:length(var_fit_names)
            variable = char(var_fit_names(i));
            mf_fitter.fit_data.(variable) = [];
            mf_fitter.averages.(variable) = [];
        end
      
        for i=1:length(var_data_names)
           variable = char(var_data_names(i));
           mf_fitter.data.(variable) = []; 
        end
        
        for i=1:length(var_cplot_names)
           variable = char(var_cplot_names(i));
           mf_fitter.center_plot.(variable) = [];
        end
              
        mf_fitter.fix = [];
        mf_fitter.SmoothedData = [];
        mf_fitter.fit_data.I_tot = [];
        mf_fitter.fit_data.I1 = [];
        
        % set the curve fit general parameters that are always the same
        status_flags.fitter.fn1d = 1;
        status_flags.fitter.function_info_1d.name = 'Gaussian';
        status_flags.fitter.function_info_1d.math_code = {'y = y0 + (i0 / (fwhm*sqrt(pi/2) / sqrt(log(4)))) .* exp((-2*(x-xc).^2) / (fwhm.^2 / log(4)));'};
        status_flags.fitter.function_info_1d.auto_guess_code = {'[amp i_amp] = max(y);    %Peak Intensity','xc = x(i_amp);         %Centre Position','y0 = min(y);         %Background','x_diff = diff(x); %Differences between adjacent x''s.','ytemp = y(1:length(y)-1);','i0 = abs((sum((ytemp-y0).*x_diff)));','sigma = i0/(amp*sqrt(2*pi));','fwhm = 2*sqrt(2*log(2))*sigma; %The sqrt is because this is a 2D gaussian','guess_values = [y0, i0, xc, fwhm];'};
        status_flags.fitter.function_info_1d.point_click_code = {'text_handle = grasp_message([''Click on Background''],1,''sub''); %Grasp Message','[x y0]=ginput(1); %Mouse input','delete(text_handle);pause(0.1);','text_handle = grasp_message([''Click on Peak''],1,''sub'');  %Grasp Message','[xc amp]=ginput(1); %Mouse input, Centre and Peak intensity','delete(text_handle);pause(0.1);','text_handle = grasp_message([''Click on Peak Width''],1,''sub'');  %Grasp Message','[fwhm y]=ginput(1); %Mouse input','delete(text_handle);pause(0.1);','fwhm=abs(fwhm-xc);','amp=amp-y0;','i0 = amp*fwhm*sqrt(pi)/(sqrt(2)*sqrt(2*log(2))); %Convert Peak Intensity to IntInt','guess_values = [y0, i0, xc, fwhm];'};
        status_flags.fitter.function_info_1d.no_parameters = 4;
        
        mf_fitter.exp = mf_fitter.folder;
  
        
    case 'set_params'
        % case run before every fit
        % called by all of the various set_ cases
        % these are parameters strictly determined by num_peaks
        status_flags.fitter.number1d = num_peaks;
        status_flags.fitter.function_info_1d.no_parameters = 4;
        status_flags.fitter.function_info_1d.variable_names = repmat({'y0', 'i0', 'xc', 'fwhm'},1,num_peaks);
        status_flags.fitter.function_info_1d.long_names = repmat({'Background', 'Integrated Intensity', 'Center', 'FWHM'},1,num_peaks);
        
        % clear out old data/options
        status_flags_names = {'group','fix','values','err_values'}
        for i=1:length(status_flags_names)
            variable = char(status_flags_names(i));
           status_flags.fitter.function_info_1d.(variable) = []; 
        end      
       
        
    case 'fit'
        % load desired depth file and update main grasp GUI
        % the +1 is necessary as file 1 represents a sum
        status_flags.selector.fd = img_num+1;
        grasp_update;

        % plot I vs xi and create the current figure, store handle
        radial_average_callbacks('averaging_control','azimuthal');
        mf_fitter.handles.plot_handle = gcf;
        grasp_update;

  
        angles = plot_info.export_data(:,1);
        intensity = plot_info.export_data(:,2);
        intensity_err = plot_info.export_data(:,3);

        mf_fitter.data.angle = angles;
        mf_fitter.data.intensity(:,img_num) = intensity';
        mf_fitter.data.intensity_err(:,img_num) = intensity_err';
        
        %Perform auto guess (only for 1-peak fit) and fit functions
        if( num_peaks==1 ), grasp_plot_fit_callbacks_2('auto_guess','auto_guess');
        else grasp_plot_fit_callbacks_2('auto_guess','off');
        end  
        
        % perform the fit
        set(grasp_handles.window_modules.curve_fit1d.curve_number,'value',1)
        grasp_plot_fit_callbacks_2('fit_it');
        
        % we don't always want all data stored
        % % store the temporary fit parameters to a permanent structure
        % mf_fitter_callbacks('data_storage',num_peaks,img_num)
        
    
    case 'data_storage'
        % set the handle based on whether grabbing the values or the errors
        for j=1:2
            if(j==1) h=status_flags.fitter.function_info_1d.values; 
            else h=status_flags.fitter.function_info_1d.err_values;
            end
       
            mf_fitter.fit_data.background(img_num,j) = h(1);
            mf_fitter.fit_data.fwhm(img_num,j) = h(4);

            if(num_peaks==1 || num_peaks==3) %grab inner peak
                if(num_peaks==1) index = 3;
                else index = 7;
                end
                mf_fitter.fit_data.intensity2(img_num,j) = h(index-1);
                mf_fitter.fit_data.center2(img_num,j) = h(index);
                mf_fitter.fit_data.intensity1(img_num,j) = 0;
                mf_fitter.fit_data.center1(img_num,j) = 0;
                mf_fitter.fit_data.intensity3(img_num,j) = 0;
                mf_fitter.fit_data.center3(img_num,j) = 0;
            end
            if(num_peaks==2 || num_peaks==3) %grab outer peaks
                mf_fitter.fit_data.intensity1(img_num,j) = h(2);
                mf_fitter.fit_data.center1(img_num,j) = h(3);
                mf_fitter.fit_data.intensity3(img_num,j) = h(4*num_peaks-2);
                mf_fitter.fit_data.center3(img_num,j) = h(4*num_peaks-1);
                if(num_peaks==2)
                    mf_fitter.fit_data.intensity2(img_num,j) = 0;
                    mf_fitter.fit_data.center2(img_num,j) = 0;
                end
            end
        end
        
        
    case 'avg'
        Sum = 0;
        Num = 0;
        for j=1:length(mf_fitter.fit_data.(variable))
            if(mf_fitter.fit_data.(variable)(j,1)~=0)
                Sum = Sum + mf_fitter.fit_data.(variable)(j,1);
                Num = Num + 1;
            end
        end
        mf_fitter.averages.(variable) = Sum / Num;

        
    case 'adv_avg'
        med = median(mf_fitter.fit_data.(variable));
        stdev = std(mf_fitter.fit_data.(variable));
        Sum = 0;
        Num = 0;
        for j=1:length(mf_fitter.fit_data.(variable))
            if(mf_fitter.fit_data.(variable)(j,1)<=med(1,1)+stdev(1,1)  && mf_fitter.fit_data.(variable)(j,1)>=med(1,1)-stdev(1,1))
                Sum = Sum + mf_fitter.fit_data.(variable)(j,1);
                Num = Num + 1;
            end
        end
        mf_fitter.averages.(variable) = Sum / Num;


    case 'discriminate'
        cutoff = mf_fitter.cutoff;

        j=1;
        TotalIntensity = 0;
        Fraction1=0;
        Fraction2=0;
        Fraction3=0;

        while (j<= mf_fitter.depth)
            TotalIntensity = abs(mf_fitter.fit_data.intensity1(j,1) +mf_fitter.fit_data.intensity2(j,1)+ mf_fitter.fit_data.intensity3(j,1));
            Fraction1= mf_fitter.fit_data.intensity1(j,1)/TotalIntensity;
            Fraction2= mf_fitter.fit_data.intensity2(j,1)/TotalIntensity;
            Fraction3= mf_fitter.fit_data.intensity3(j,1)/TotalIntensity;

            if (Fraction1<cutoff && Fraction3<cutoff)
                mf_fitter.outer_peaks(j,1) = 0;
            else 
                mf_fitter.outer_peaks(j,1) = 1;
            end

            if (Fraction2<cutoff)
                mf_fitter.inner_peak(j,1) = 0;
            else 
                mf_fitter.inner_peak(j,1) = 1;
            end
            j=j+1;
        end

    case 'center_store'
        mf_fitter.center_plot.x1 = mf_fitter.fit_data.center1(:,1);
        mf_fitter.center_plot.x1_err = mf_fitter.fit_data.center1(:,2);
        mf_fitter.center_plot.x2 = mf_fitter.fit_data.center2(:,1);
        mf_fitter.center_plot.x2_err = mf_fitter.fit_data.center2(:,2);
        mf_fitter.center_plot.x3 = mf_fitter.fit_data.center3(:,1);
        mf_fitter.center_plot.x3_err = mf_fitter.fit_data.center3(:,2);

    case 'center_plot'
        mf_fitter_callbacks('adv_avg',0,0,'center2');
        xc = mf_fitter.averages.center2;
        mf_fitter.handles.cplot = figure
        axes1 = axes('Parent',mf_fitter.handles.cplot,...
            'ZColor',[0 0 0],'YColor',[0 0 0],'XColor',[0 0 0],...
            'YLim',[xc-20, xc+20],...
            'Color',[1 1 1]); 
         hold(axes1,'all');
        
        if( isfield(mf_fitter.fit_data, 'cycles'))
            x = mf_fitter.fit_data.cycles';
            xaxis = 'AC Cycles';
            set(gca,'XScale','log');
        else
            x = 1:mf_fitter.depth;
            xaxis = 'Index';
        end
        
         titlename = ['\fontsize{24}\color{black}\fontname{Arial} ' mf_fitter.exp ' Peak Centers AB: ' get(grasp_handles.window_modules.radial_average.azimuth_bin,'String')]
         title(titlename);
         xlabel(['\fontsize{16}\color{black}\fontname{Arial} ' xaxis]);
         ylabel(['\fontsize{16}\color{black}\fontname{Arial} Azimuthal Angle (Degrees)']);
        
        
        % find centers that are > than intensity cutoff
        m=1;
        g=1;
        ic = 0.05;
        for i=1:length(mf_fitter.fit_data.intensity1)
            igs1 = mf_fitter.fit_data.intensity1(i);
            ims = mf_fitter.fit_data.intensity2(i);
            igs2 = mf_fitter.fit_data.intensity3(i);
            
            if ( (ims / (igs1+igs2+ims)) > (ic) )
                x2(m) = mf_fitter.center_plot.x2(i);
                x2_err(m) = mf_fitter.center_plot.x2_err(i);
                xms(m) = x(i);
                m = m+1;
            end
            if ( (igs1 / (igs1+igs2+ims)) > (ic) && (igs2 / (igs1+igs2+ims)) > (ic)  )
                x1(g) = mf_fitter.center_plot.x1(i);
                x1_err(g) = mf_fitter.center_plot.x1_err(i);
                x3(g) = mf_fitter.center_plot.x3(i);
                x3_err(g) = mf_fitter.center_plot.x3_err(i);
                xgs(g) = x(i);
                g = g+1;
            end
        end
        
        % plot the center data with error bars 
        % errorbar(x,y,L,U)
        hold on
        errorbar(x,mf_fitter.fit_data.center1(:,1),mf_fitter.fit_data.center1(:,2),'k*','MarkerSize',7);
        errorbar(x,mf_fitter.fit_data.center2(:,1),mf_fitter.fit_data.center2(:,2),'k*','MarkerSize',7);
        errorbar(x,mf_fitter.fit_data.center3(:,1),mf_fitter.fit_data.center3(:,2),'k*','MarkerSize',7);
        errorbar(xgs,x1,x1_err,'bo','MarkerSize',7);
        errorbar(xms,x2,x2_err,'ro','MarkerSize',7);
        errorbar(xgs,x3,x3_err,'bo','MarkerSize',7);
        hold off;
        
        % store the values of the centers/cycles that are above the cutoff for later use
        mf_fitter.center_plot.xgs = xgs;
        mf_fitter.center_plot.xms = xms;
        mf_fitter.center_plot.x1 = x1;
        mf_fitter.center_plot.x2 = x2;
        mf_fitter.center_plot.x3 = x3;
        mf_fitter.center_plot.x1_err = x1_err;
        mf_fitter.center_plot.x2_err = x2_err;
        mf_fitter.cetner_plot.x3_err = x3_err;

       
    case 'refit'
        %Parameters that vary based on number of peaks (num_peaks)
        try
            close(grasp_handles.window_modules.curve_fit1d.window);
        end
        grasp_plot_fit_window()
        %mf_fitter_callbacks('set_params',num_peaks,img_num);
        
        status_flags.fitter.number1d = num_peaks;
        status_flags.fitter.function_info_1d.no_parameters = 4;
        status_flags.fitter.function_info_1d.variable_names = repmat({'y0', 'i0', 'xc', 'fwhm'},1,num_peaks);
        status_flags.fitter.function_info_1d.long_names = repmat({'Background', 'Integrated Intensity', 'Center', 'FWHM'},1,num_peaks);
        status_flags.fitter.function_info_1d.group = repmat([1,0,0,1],1,num_peaks);
        status_flags.fitter.function_info_1d.err_values = repmat([0,0,0,0],1,num_peaks);
        status_flags.fitter.function_info_1d.fix = repmat([1],1,4*num_peaks)
      
        %First give guess values to be overwritten if should be fixed
        if(num_peaks~=1)
           status_flags.fitter.function_info_1d.values(1) = mf_fitter.fit_data.background(img_num,1);
           status_flags.fitter.function_info_1d.values(4) = mf_fitter.fit_data.fwhm(img_num,1);
           status_flags.fitter.function_info_1d.values(3) = mf_fitter.fit_data.center1(img_num,1);
           status_flags.fitter.function_info_1d.values(7) = mf_fitter.fit_data.center3(img_num,1);
           status_flags.fitter.function_info_1d.values(2) = mf_fitter.fit_data.intensity1(img_num,1);
           status_flags.fitter.function_info_1d.values(6) = mf_fitter.fit_data.intensity3(img_num,1);
            if(num_peaks==3)
                status_flags.fitter.function_info_1d.values(7) = mf_fitter.fit_data.center2(img_num,1);
                status_flags.fitter.function_info_1d.values(11) = mf_fitter.fit_data.center3(img_num,1);
                status_flags.fitter.function_info_1d.values(6) = mf_fitter.fit_data.intensity2(img_num,1);
                status_flags.fitter.function_info_1d.values(10) = mf_fitter.fit_data.intensity3(img_num,1);
            end
        end
        
         grasp_update();
         grasp_plot_fit_callbacks('update_curve_fit_window');
         grasp_plot_fit_window()
         mf_fitter_callbacks('fit',3,img_num)
         
    case 'refit-Matlab'
        % get data
        phi = mf_fitter.SmoothedData.phi;
        Int = mf_fitter.SmoothedData.Int(img_num,:);
        Int_err = mf_fitter.SmoothedData.Int_err(img_num,:);
        
        figure
        errorbar(phi,Int,Int_err,'ok','MarkerFaceColor',[1 1 1])
        hold on
        plot_template
        hold on
        
        % get fit
        y0 = mf_fitter.fit_data.background(img_num,1);
        itot = mf_fitter.fit_data.I_tot(img_num,1);
        i1 = mf_fitter.fit_data.I1(img_num,1);
        fwhm = mf_fitter.fit_data.fwhm(img_num,1);
        xc1 = mf_fitter.fit_data.center1(img_num,1);
        xc3 = mf_fitter.fit_data.center3(img_num,1);
        
        % plot fit
        fplot(@(x) y0 + itot*((i1/(fwhm * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(x-xc1)^2/(fwhm^2/log(4))) + ((1-i1)/(fwhm * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(x-xc3)^2/(fwhm^2/log(4)))), [min(phi) max(phi)], 'r')
        hold off              
        
    
    case 'fraction'
         %create actual plot
         %turn hold on so it will add each data point as calculated
         mf_fitter.handles.fig = axes('Parent',figure,...
            'ZColor',[0 0 0],'YColor',[0 0 0],'XColor',[0 0 0],...
            'YLim',[0, 1.2],...
            'Color',[1 1 1]); 
         hold(mf_fitter.handles.fig,'all');
         
         titlename = [mf_fitter.exp ' Fractional Volume in One Peak State AB: ' get(grasp_handles.window_modules.radial_average.azimuth_bin,'String')];
         title(titlename,'FontSize',16,'Fontname','Arial','Color','black');
         xlabel(['\fontsize{16}\color{black}\fontname{Arial} Index']);
         ylabel(['\fontsize{16}\color{black}\fontname{Arial} Fraction of Vortices']);
         hold on;

         %loop through fit_data structure and add relative intensity column
        for j=1:mf_fitter.depth
           mf_fitter.fit_data.fms(j) = mf_fitter.fit_data.intensity2(j) / (mf_fitter.fit_data.intensity1(j) + mf_fitter.fit_data.intensity2(j) + mf_fitter.fit_data.intensity3(j));
           mf_fitter.fit_data.fgs(j) = 1 - mf_fitter.fit_data.fms(j);
        end
        
        indices = 1:mf_fitter.depth;
        % calculate the error
        mf_errorprop();
        
        if( isfield(mf_fitter.fit_data, 'cycles'))
            x = mf_fitter.fit_data.cycles';
            xlabelname = ['\fontsize{16}\color{black}\fontname{Arial} Number of Applied AC Cycles'];
            set(gca,'XScale','log');
        else
            x = indices;
            xlabelname = ['\fontsize{16}\color{black}\fontname{Arial} File Index'];
        end
        
        xlabel(['\fontsize{16}\color{black}\fontname{Arial}' xlabelname]);
        
        % plot the fms data with error bars 
        errorbar(x,mf_fitter.fit_data.fms,mf_fitter.error_prop.sig_frac,'bo','MarkerSize',7);
        set(gca,'YScale','log');
        
        hold off;
         
    
    case 'colormap'     
    % Code for smoothing out the data fluctuations for the 
    % initialize
        phi = mf_fitter.data.angle;
        intensity = mf_fitter.data.intensity;
        xc1 = mf_fitter.center_plot.x1;
        xc2 = mf_fitter.center_plot.x2;
        xc3 = mf_fitter.center_plot.x3;
        xgs = mf_fitter.center_plot.xgs;
        xms = mf_fitter.center_plot.xms;
        cycles = mf_fitter.fit_data.cycles;
        
        climits = [mf_fitter.climits(1), mf_fitter.climits(2)]
     
        num_angles = length(phi);
        cyc = length(cycles);
        bg = 0;
        xcycles = [];
        intensity(num_angles+1,:) = 0;  % gives junk column/row not plotted by pcolor
        intensity(:,cyc+1) = 0;
        backsub_int = zeros(size(intensity));  
        integrated_intensity = zeros(size(intensity));
        smoothbox_I = zeros(size(intensity));
        smoothgaus_I = zeros(size(intensity));
        

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
        
    % Smoothing

    % Box Smoothing
    % Take the mean of the pixel and it's two nearest neighbors
        for i = 1:cyc  
            for j=1:num_angles
                if(j==1) smoothbox_I(j,i) = mean(integrated_intensity((j):(j+1),i));
                elseif (j==num_angles) smoothbox_I(j,i) = mean(integrated_intensity((j-1):(j),i));
                else smoothbox_I(j,i) = mean(integrated_intensity((j-1):(j+1),i));
                end
            end
        end

    % Gaussian
    % Take the mean again, this time weighted by the gaussian with specified fwhm 
        fwhm = 0.5;
        for i = 1:cyc
            total_I = [];
            sum_weights = [];
            for j=1:num_angles
                for pix=1:num_angles
                    total_I(pix) = exp((-(j-pix)^2)/fwhm)*integrated_intensity(pix,i);
                    sum_weights(pix) = exp((-(j-pix)^2)/fwhm);
                end
                smoothgaus_I(j,i) = sum(total_I)/sum(sum_weights);
            end
        end

    % For Logarithmic intensity scale    
        log_i = log(smoothgaus_I-min(min(smoothgaus_I)));
        
    % Add junk cycle and phi value (based on the way matlab command pcolor works)
        cycles(cyc+1) = 2*cycles(cyc);
        %xgs(length(xgs)+1) = 2*xgs(length(xgs));
        %xms(length(xms)+1) = 2*xms(length(xms));
        phi_shift(num_angles+1) = phi_shift(num_angles) + delta_phi(num_angles-1);
        phi(num_angles+1) = phi(num_angles) + delta_phi(num_angles-1);
        
    % Define values for asprep column    
        asprep_log_i =  log_i(:,1:2);
        asprep_cycles = [0,1];
        
        for i=1:(length(xgs)-1)
            xgscycles(i) = sqrt(xgs(i+1)*xgs(i));
        end
        for i=1:(length(xms)-1)
            xmscycles(i) = sqrt(xms(i+1)*xms(i));
        end
        xgscycles(length(xgs)) = sqrt(xgs(length(xgs))*2*xgs(length(xgs)));
        xmscycles(length(xms)) = sqrt(xms(length(xms))*2*xms(length(xms)));
        
    % Figure Formatting (gives the correct aspect ratio)
        fig_size = 2*[7.3 5.7];
        offset = [0.1 0.1];
        fig_start = [0 0];
        graph_start = [1.25 1.25];
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
            s = sprintf('Azimuthal Angle (%c)', char(176));
            ylabel(s);
      
    % Make the rest of the colormap figure        
        h = subplot(1,2,2);
            h = pcolor(cycles,phi_shift,log_i); 
            %str = {'Logarithmic Integrated Intensity', ['\fontsize{10}', sprintf(' with %d fwhm Gaussian Smoothing', fwhm)], '\fontsize{10} 1.29 Angle Binning'}
            xlabel('Applied || AC Cycles');
            hold on
            % Add the fitted data centers
                plot(xgscycles,xc1,'.','MarkerSize',10,'MarkerEdgeColor',[1 1 1]);
                plot(xmscycles,xc2,'.','MarkerSize',10,'MarkerEdgeColor',[1 1 1]);
                plot(xgscycles,xc3,'.','MarkerSize',10,'MarkerEdgeColor',[1 1 1]);
            hold off
  
    % Format the figures
        subplot(1,2,1); 
        set(gca,...
            'Color',[1 1 1],...
            'FontName','Arial','FontSize',8,...
            'Xcolor',[0 0 0],'Ycolor',[0 0 0],...
            'XTickLabelMode','Manual','XTickLabel',{' ','0',' '},...
            'CLim', climits,...
            'Units','centimeters','Position', position1);
            ax.Box = 'on';
            
        h_ylabel = get(gca,'YLabel');
        set(h_ylabel,'FontName','Arial','FontSize',8,'Color',[0 0 0]);
 
        subplot(1,2,2);
        set(gca,...
            'Color',[1 1 1],...
            'FontName','Arial','FontSize',8,...
            'Xcolor',[0 0 0],'Ycolor',[0 0 0],...
            'YTickLabelMode','Manual','YTickLabel',{},...
            'XTickLabelMode','Auto','XTick',[10,100,1000,10000],...
            'XScale','Log','CLim', climits,...
            'Units','centimeters','Position', position2);
         ax.Box = 'on';
             
        h_xlabel = get(gca,'XLabel');
        set(h_xlabel,'FontName','Arial','FontSize',8,'Color',[0 0 0]);
        
        titlename = [mf_fitter.exp ' Colormap ' get(grasp_handles.window_modules.radial_average.azimuth_bin,'String')];
        title(titlename,'FontSize',12,'Fontname','Arial','Color','black');
       
        
        
    case 'peak_decay'
        disp('Peak Decay Plot')
        mf_fitter.handles.fig = figure('PaperSize',[8.3 11.7],...
            'Color',[0.80 0.80 0.80]);

        hold on
        mf_fitter_callbacks('avg',0,0,'center1');
        mf_fitter_callbacks('avg',0,0,'center2');
        mf_fitter_callbacks('avg',0,0,'center3');
        
        x1 = mf_fitter.averages.center1;
        x2 = mf_fitter.averages.center2;
        x3 = mf_fitter.averages.center3;
        xlow = x2-20;
        xhigh = x2+20;
        
        scale = 3;
        
         subplot2 = subplot(1,1,1,'Parent',mf_fitter.handles.fig,...
            'XLim', [xlow, xhigh],...
            'XDir','reverse',...
            'YLim', [-(scale*mf_fitter.depth), 10],...
            'Fontname','Timesnewroman',...
            'FontSize',12);
        
        set(gca,'Color',[1, 1, 1], 'XColor',[0, 0, 0], 'YColor',[0, 0, 0]);
        
        phi = mf_fitter.data.angle;
        for i = 1:mf_fitter.depth
           intensity2 = mf_fitter.data.intensity(:,i) - scale*i;
           errorbar(phi,intensity2,mf_fitter.data.intensity_err(:,i),...
                              'bo','LineWidth',1,...
                              'MarkerEdgeColor','k',...
                              'MarkerFaceColor','k',...
                              'MarkerSize',5);
                          
           y0 = mf_fitter.fit_data.background(i)-scale*i;
           fwhm = mf_fitter.fit_data.fwhm(i);
           i01 = mf_fitter.fit_data.intensity1(i);
           i02 = mf_fitter.fit_data.intensity2(i);
           i03 = mf_fitter.fit_data.intensity3(i);
           x1 = mf_fitter.fit_data.center1(i);
           x2 = mf_fitter.fit_data.center2(i);
           x3 = mf_fitter.fit_data.center3(i);
   
           fplot(@(x) y0 + (i01/(fwhm * sqrt(pi/2) / sqrt(log(4))))*exp(-2*((x-x1)^2/(fwhm^2/log(4)))) + (i02/(fwhm * sqrt(pi/2) / sqrt(log(4))))*exp(-2*((x-x2)^2/(fwhm^2/log(4)))) + (i03/(fwhm * sqrt(pi/2) / sqrt(log(4))))*exp(-2*((x-x3)^2/(fwhm^2/log(4))) ),...
                        [xlow xhigh], 2, 'r'); 
        end
        titlename = [mf_fitter.exp ' Peak Decay ' get(grasp_handles.window_modules.radial_average.azimuth_bin,'String')];
        title(titlename,'FontSize',16,'Fontname','Arial','Color','black');
        xlabel('Azimuthal Angle','FontSize',12,'Fontname','Arial','Color','black');
        ylabel('Relative Intensity','FontSize',12,'Fontname','Arial','Color','black');
       
        set(gca,'View',[-90,90]);
        set(gca,'YLim',[-scale*mf_fitter.depth-5,scale]);

        hold off
        
        
    case 'update'
        mf_fitter_callbacks(data_storage,num_peaks,img_num);
        mf_fitter_table();
   
        
%Close the switch   



%% For new way of fitting
    case 'matlab_fit1'
        % specify fit options and type
        % order of coeffecients for boundaries:
        %   a, i01, i02, i03, y0

        N = mf_fitter.depth;
        i = img_num;
        delta12 = mf_fitter.fit_data.center2(1) - mf_fitter.fit_data.center1(N);
        delta23 = mf_fitter.fit_data.center3(N) - mf_fitter.fit_data.center2(1);
        
        a = (delta12+delta23)/2;
        i01M = mf_fitter.fit_data.intensity1(N);
        i02M = mf_fitter.fit_data.intensity2(1);
        i03M = mf_fitter.fit_data.intensity3(N);
        
        xc = mf_fitter.averages.center2;
        fwhm = mf_fitter.averages.fwhm;
        y0 = mf_fitter.averages.background;

        start = [i/N * a, i/N*i01M, (N-i)/N*i02M, i/N*i03M];
        lower = [(1.5*fwhm)/2, 0, 0, 0];
        upper = [a, 1.2*i01M, 1.2*i02M, 1.2*i03M];

        fo = fitoptions('Method','NonlinearLeastSquares','StartPoint',start,'Lower',lower,'Upper',upper);

 %       ft = fittype('y0 + (i01/(fwhm * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(x-(xc-a))^2/(fwhm^2/log(4))) + (i02/(fwhm * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(x-xc)^2/(fwhm^2/log(4))) + (i03/(fwhm * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(x-(xc+a))^2/(fwhm^2/log(4)))',...
  %                      'options',fo,'problem',{'y0','fwhm','xc'});   
  
    
         ft = fittype('y0 + i0*((i1/(fwhm * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(x-(xc-a))^2/(fwhm^2/log(4))) + ((1-i1)/(fwhm * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(x-xc)^2/(fwhm^2/log(4))))',...
                        'options',fo,'problem',{'y0','fwhm','xc'});        


        % load desired depth file and update main grasp GUI
        % the +1 is necessary as file 1 represents a sum
        status_flags.selector.fd = img_num+1;
        grasp_update;

        % plot I vs xi and create the current figure, store handle
        radial_average_callbacks('averaging_control','azimuthal');
        mf_fitter.handles.plot_handle = gcf;
        grasp_update;
            
        angles = plot_info.export_data(:,1);
        intensity = plot_info.export_data(:,2);
        intensity_err = plot_info.export_data(:,3);
       
        [fitobject, gof] = fit(angles,intensity,ft,'problem',{y0,fwhm,xc})
        ci = confint(fitobject);
        plot(fitobject,angles,intensity);  
        
        mf_fitter.data.angle = angles;
        mf_fitter.data.intensity(:,img_num) = intensity;
        mf_fitter.data.intensity_err(:,img_num) = intensity_err;
        
        mf_fitter.fit_data.background(img_num,1) = y0;
        mf_fitter.fit_data.fwhm(img_num,1) = fwhm;
        %mf_fitter.fit_data.intensity1(img_num,1) = fitobject.i01;
        %mf_fitter.fit_data.intensity2(img_num,1) = fitobject.i02;
        %mf_fitter.fit_data.intensity3(img_num,1) = fitobject.i03;
        mf_fitter.fit_data.intensity1(img_num,1) = fitobject.i0;
        mf_fitter.fit_data.intensity2(img_num,1) = fitobject.i1;
        mf_fitter.fit_data.center1(img_num,1) = xc - fitobject.a;  
        mf_fitter.fit_data.center1(img_num,2) = ((fitobject.a - ci(1)) + (ci(2) - fitobject.a))/2;
        mf_fitter.fit_data.center2(img_num,1) = xc;
        mf_fitter.fit_data.center2(img_num,2) = ((fitobject.a-ci(1)) + (ci(2)-fitobject.a))/2; % + (((fitobject.b-ci(3))+(ci(4)-fitobject.b))/2)^2 );
        mf_fitter.fit_data.center3(img_num,1) = xc + fitobject.a; 
    
        
        case 'matlab_fit4'
         % specify fit options and type
        % order of coeffecients for boundaries:
        %   a, i01, i02, i03, y0

        N = mf_fitter.depth;
        i = img_num;

        i01 = mf_fitter.fit_data.intensity1(i);
        i02 = mf_fitter.fit_data.intensity2(i);
        i03 = mf_fitter.fit_data.intensity3(i);
        
        x1 = mf_fitter.fit_data.center1(i);
        x2 = mf_fitter.averages.center2;
        x3 = mf_fitter.fit_data.center3(i);
       
        fwhm = mf_fitter.averages.fwhm;
        y0 = mf_fitter.averages.background;

        xstd = 0.5;
        start = [x1, x3];
        lower = [ x1-xstd, x3-xstd];
        upper = [ x1+xstd, x3+xstd];
    
        fo = fitoptions('Method','NonlinearLeastSquares','StartPoint',start,'Lower',lower,'Upper',upper);

        ft = fittype('y0 + (i01/(fwhm * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(x-x1)^2/(fwhm^2/log(4))) + (i02/(fwhm * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(x-x2)^2/(fwhm^2/log(4))) + (i03/(fwhm * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(x-x3)^2/(fwhm^2/log(4)))',...
                        'options',fo,'problem',{'y0','fwhm','i01','i02','i03','x2'});        

        % load desired depth file and update main grasp GUI
        % the +1 is necessary as file 1 represents a sum
        status_flags.selector.fd = img_num+1;
        grasp_update;

        % plot I vs xi and create the current figure, store handle
        radial_average_callbacks('averaging_control','azimuthal');
        mf_fitter.handles.plot_handle = gcf;
        grasp_update;
            
        angles = plot_info.export_data(:,1);
        intensity = plot_info.export_data(:,2);
        intensity_err = plot_info.export_data(:,3);
       
        [fitobject, gof] = fit(angles,intensity,ft,'problem',{y0,fwhm, i01, i02, i03,x2});
        ci = confint(fitobject)
        plot(fitobject,angles,intensity);
        
        mf_fitter.data.angle = angles;
        mf_fitter.data.intensity(:,img_num) = intensity;
        mf_fitter.data.intensity_err(:,img_num) = intensity_err;
        
        mf_fitter.fit_data.background(img_num,1) = y0;
        mf_fitter.fit_data.fwhm(img_num,1) = fwhm;
        mf_fitter.fit_data.intensity1(img_num,1) = i01;
        mf_fitter.fit_data.intensity2(img_num,1) = i02;
        mf_fitter.fit_data.intensity3(img_num,1) = i03;
        mf_fitter.fit_data.center1(img_num,1) = fitobject.x1;
        mf_fitter.fit_data.center2(img_num,1) = x2;
        mf_fitter.fit_data.center3(img_num,1) = fitobject.x3;
        mf_fitter.fit_data.center1(img_num,2) = ((fitobject.x1 - ci(1)) + (ci(2) - fitobject.x1))/2;
        mf_fitter.fit_data.center2(img_num,2) = 0;
        mf_fitter.fit_data.center3(img_num,2) = ((fitobject.x3 - ci(3)) + (ci(4) - fitobject.x3))/2;

        
  case 'set_mfp1'
        % sets curve fit window
        % never fixes the background or intensities
        % cycles through fixing the GS centers, the fwhm, the MS center, or all of the above
        % groups background, fwhm

        i = img_num;
        mf_fitter_callbacks('set_params',num_peaks,i);
        
        status_flags.fitter.function_info_1d.group = repmat([1,0,0,1],1,num_peaks);
        status_flags.fitter.function_info_1d.err_values = repmat([0,0,0,0],1,num_peaks);         
        status_flags.fitter.function_info_1d.values = repmat([0],1,12);
        
        if(variable==1)
            status_flags.fitter.function_info_1d.fix = [0,0,1,1,1,0,1,1,1,0,1,1];
            status_flags.fitter.function_info_1d.values(7) = mf_fitter.averages.center2;
        elseif(variable==2)
            status_flags.fitter.function_info_1d.fix = [0,0,1,0,1,0,1,1,1,0,1,1];
            status_flags.fitter.function_info_1d.values(7) = mf_fitter.averages.center2;
        elseif(variable==3)
            status_flags.fitter.function_info_1d.fix = [0,0,1,1,1,0,0,1,1,0,1,1];
            status_flags.fitter.function_info_1d.values(7) = mf_fitter.fit_data.center2(i);
       elseif(variable==4)
            status_flags.fitter.function_info_1d.fix = [0,0,1,1,1,0,1,1,1,0,1,1];
            status_flags.fitter.function_info_1d.values(7) = mf_fitter.averages.center2;
        elseif(variable==5)
            status_flags.fitter.function_info_1d.fix = [0,0,1,1,1,0,1,1,1,0,1,1];
            status_flags.fitter.function_info_1d.values(7) = mf_fitter.fit_data.center2(i);
        end

        status_flags.fitter.function_info_1d.values(1) = mf_fitter.averages.background;
        status_flags.fitter.function_info_1d.values(2) = mf_fitter.fit_data.intensity1(i);
        status_flags.fitter.function_info_1d.values(3) = mf_fitter.fix.center1(i);
        status_flags.fitter.function_info_1d.values(4) = mf_fitter.averages.fwhm;
        status_flags.fitter.function_info_1d.values(6) = mf_fitter.fit_data.intensity2(i);
        status_flags.fitter.function_info_1d.values(10) = mf_fitter.fit_data.intensity3(i);
        status_flags.fitter.function_info_1d.values(11) = mf_fitter.fix.center3(i);

            
        case 'save'
            current_ab = str2num(get(grasp_handles.window_modules.radial_average.azimuth_bin,'String'));
            %mf_fitter.folder = ['AngleBinning_' num2str(current_ab) '/'];
            name = [mf_fitter.file ' AB ' num2str(current_ab)];
            dir = [mf_fitter.extension mf_fitter.folder];
            mkdir(dir);
            if(variable==2)
                figfile = [mf_fitter.extension mf_fitter.folder name '.fig'];
                saveas(mf_fitter.handles.fig,figfile);
                pdffile = [mf_fitter.extension mf_fitter.folder name '.pdf'];
                saveas(mf_fitter.handles.fig,pdffile);
            else
                jpgfile = [mf_fitter.extension mf_fitter.folder name '.jpg'];
                saveas(mf_fitter.handles.fig,jpgfile);
            end
            
            
    case 'plot_fit'
        %% Parameters
        i = img_num
        c1 = mf_fitter.fit_data.center1(i);
        c2 = mf_fitter.fit_data.center2(i);
        c3 = mf_fitter.fit_data.center3(i);
        i01 = mf_fitter.fit_data.intensity1(i);
        i02 = mf_fitter.fit_data.intensity2(i);
        i03 = mf_fitter.fit_data.intensity3(i);
        fwhm = mf_fitter.fit_data.fwhm(i);
        y0 = mf_fitter.fit_data.background(i);
        name = mf_fitter.fit_data.names(i);
        limits = [c2-20 c2+20];


        %% Fixed Peaks
        figure
        hold on
        label = ['Fit ' num2str(name)]
        title(label) % FontSize',16,'Fontname','Timesnewroman','Color','black')
        xlabel('Azimuthal Angle')%,'FontSize',12,'Fontname','Timesnewroman','Color','black')
        ylabel('Counts / standard monitor') %,'FontSize',12,'Fontname','Timesnewroman','Color','black')

        errorbar(mf_fitter.data.angle, mf_fitter.data.intensity(:,i), mf_fitter.data.intensity_err(:,i), 'bo')
        scatter(mf_fitter.data.angle, mf_fitter.data.intensity(:,i), 'w.')
        fplot(@(x) y0 + (i01 / (fwhm*sqrt(pi/2) / sqrt(log(4))) )*exp(-2*(x-c1)^2/(fwhm^2 / log(4))) + (i02 / (fwhm*sqrt(pi/2) / sqrt(log(4)) ))*exp(-2*(x-c2)^2/(fwhm^2 / log(4)))  + (i03 / (fwhm*sqrt(pi/2) / sqrt(log(4)) ))*exp(-2*(x-c3)^2/(fwhm^2 / log(4))), limits, '-g');

        plot_template()


    %% Currently unused callbacks
    case 'set_mfp_old'
        % sets curve fit window
        % fixes neither the fwhm nor the centers
        % groups background, fwhm
        % runs on the reference peaks (one or two peak fit)

        mf_fitter_callbacks('set_params',num_peaks,img_num);
        
        status_flags.fitter.function_info_1d.group = repmat([1,0,0,1],1,num_peaks);
        status_flags.fitter.function_info_1d.err_values = repmat([0,0,0,0],1,num_peaks);

        i01M = mf_fitter.fit_data.intensity1(mf_fitter.good_peaks.outers_num,1);
        i02M = mf_fitter.fit_data.intensity2(mf_fitter.good_peaks.inner_num,1)
        i03M = mf_fitter.fit_data.intensity3(mf_fitter.good_peaks.outers_num,1);
        x1M = mf_fitter.fit_data.center1(mf_fitter.good_peaks.outers_num,1);
        xc = mf_fitter.averages.center2;
        x3M = mf_fitter.fit_data.center3(mf_fitter.good_peaks.outers_num,1);
        fwhm = mf_fitter.averages.fwhm;
        i = img_num
        N = mf_fitter.depth
        
        i01M*(i/N)
        i
        N
        x1M
        xc
        fwhm
        (((i-1)/N)*(x1M)+(xc-(N-i-1)/N*(fwhm/2)))
        i02M*((N-i)/N)
        i03M*(i/N)
        (((i-1)/N)*(x3M)+(xc+(N-i-1)/N*(fwhm/2)))
                
            status_flags.fitter.function_info_1d.fix = [0,0,0,1,1,0,1,1,1,0,0,1];
            status_flags.fitter.function_info_1d.values = repmat([0],1,12);
            status_flags.fitter.function_info_1d.values(1) = mf_fitter.averages.background;
            status_flags.fitter.function_info_1d.values(2) = i01M*(i/N);
            status_flags.fitter.function_info_1d.values(3) = ((i/N)*(x1M)+(xc-(N-i)/N*(fwhm/2)))/2
            status_flags.fitter.function_info_1d.values(4) = fwhm;
            status_flags.fitter.function_info_1d.values(6) = i02M*(N-i)/N;
            status_flags.fitter.function_info_1d.values(7) = xc;
            status_flags.fitter.function_info_1d.values(10) = i03M*(i/N);
            status_flags.fitter.function_info_1d.values(11) = ((i/N)*(x3M)+(xc+(N-i)/N*(fwhm/2)))/2
            
            
             case 'matlab_fit3'
         % specify fit options and type
        % order of coeffecients for boundaries:
        %   a, i01, i02, i03, y0

        N = mf_fitter.depth;
        i = img_num;

        i01 = mf_fitter.fit_data.intensity1(i);
        i02 = mf_fitter.fit_data.intensity2(i);
        i03 = mf_fitter.fit_data.intensity3(i);
        
       x1 = mf_fitter.fit_data.center1(i);
       x2 = mf_fitter.fit_data.center2(i);
       x3 = mf_fitter.fit_data.center3(i);
       
        fwhm = mf_fitter.averages.fwhm;
        y0 = mf_fitter.averages.background;

        %x2std = std(mf_fitter.fit_data.center2);
        x2std = 0.5;
        xstd = std(mf_fitter.fit_data.center1);
        xstd = 0.5;
        start = [i01, i02, i03, x1, x2, x3];
        lower = [0.8*i01, 0.8*i02, 0.8*i03, x1-xstd, x2-x2std, x3-xstd];
        upper = [1.2*i01, 1.2*i02, 1.2*i03, x1+xstd, x2+x2std, x3+xstd];

        if( i01 < 0 )
            start(1) = 0;
            lower(1) = -5;
            upper(1) = 5;
        end
        if( i03 < 0 )
            start(3) = 0;
            lower(3) = -5;
            upper(3) = 5;
        end
    
        fo = fitoptions('Method','NonlinearLeastSquares','StartPoint',start,'Lower',lower,'Upper',upper);

        ft = fittype('y0 + (i01/(fwhm * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(x-x1)^2/(fwhm^2/log(4))) + (i02/(fwhm * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(x-x2)^2/(fwhm^2/log(4))) + (i03/(fwhm * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(x-x3)^2/(fwhm^2/log(4)))',...
                        'options',fo,'problem',{'y0','fwhm'});        

        % load desired depth file and update main grasp GUI
        % the +1 is necessary as file 1 represents a sum
        status_flags.selector.fd = img_num+1;
        grasp_update;

        % plot I vs xi and create the current figure, store handle
        radial_average_callbacks('averaging_control','azimuthal');
        mf_fitter.handles.plot_handle = gcf;
        grasp_update;
            
        angles = plot_info.export_data(:,1);
        intensity = plot_info.export_data(:,2);
        intensity_err = plot_info.export_data(:,3);
       
        [fitobject, gof] = fit(angles,intensity,ft,'problem',{y0,fwhm});
        ci = confint(fitobject)
        plot(fitobject,angles,intensity); 
        
        mf_fitter.data.angle = angles;
        mf_fitter.data.intensity(:,img_num) = intensity;
        mf_fitter.data.intensity_err(:,img_num) = intensity_err;
                
        mf_fitter.fit_data.background(img_num,1) = y0;
        mf_fitter.fit_data.fwhm(img_num,1) = fwhm;
        mf_fitter.fit_data.intensity1(img_num,1) = fitobject.i01;
        mf_fitter.fit_data.intensity2(img_num,1) = fitobject.i02;
        mf_fitter.fit_data.intensity3(img_num,1) = fitobject.i03;
        mf_fitter.fit_data.center1(img_num,1) = fitobject.x1;
        mf_fitter.fit_data.center2(img_num,1) = fitobject.x2;
        mf_fitter.fit_data.center3(img_num,1) = fitobject.x3;
        mf_fitter.fit_data.intensity1(img_num,2) = ((fitobject.i01 - ci(1)) + (ci(2) - fitobject.i01) )/2
        mf_fitter.fit_data.intensity2(img_num,2) = ((fitobject.i02 - ci(3)) + (ci(4) - fitobject.i02) )/2
        mf_fitter.fit_data.intensity1(img_num,2) = ((fitobject.i03 - ci(5)) + (ci(6) - fitobject.i03) )/2
        mf_fitter.fit_data.center1(img_num,2) = ((fitobject.x1 - ci(7)) + (ci(8) - fitobject.x1))/2;
        mf_fitter.fit_data.center2(img_num,2) = ((fitobject.x2 - ci(9)) + (ci(10) - fitobject.x2))/2;
        mf_fitter.fit_data.center3(img_num,2) = ((fitobject.x3 - ci(11)) + (ci(12) - fitobject.x3))/2;

case 'matlab_fit_cdet'
        % specify fit options and type
        % order of coeffecients for boundaries:
        %   a, b, i01, i02, i03, xc, y0

        %tol = str2num(get(mf_fitter.handles.tolerance,'String'));
        
        x1 = mf_fitter.fit_data.center1(img_num);
        x2 = mf_fitter.fit_data.center2(img_num);
        x3 = mf_fitter.fit_data.center3(img_num);
        
        fwhm = mf_fitter.averages.fwhm;
        i01 = mf_fitter.fit_data.intensity1(img_num);
        i02 = mf_fitter.fit_data.intensity2(img_num);
        i03 = mf_fitter.fit_data.intensity3(img_num);
        y0 = mf_fitter.fit_data.background(img_num);
        
        start = [x1, x2, x3];
        lower = start-0.75;
        upper = start+0.75;
% 
%         if( abs(0.3*i01) < 2*10^(-8) )
%             start = [a, b, abs(i01), abs(i02), abs(i03), xc, y0];
%             lower = [a-1, 0, -2, 0.85*i02, 0.85*i03, xc-0.75, 0];
%             upper = [a+1, tol, 5, 1.15*abs(i02), 1.15*abs(i03), xc+0.75, 2*y0]; 
%         elseif( abs(0.3*i03) < 2*10^(-8) )
%             start = [a, b, abs(i01), abs(i02), abs(i03), xc, y0];
%             lower = [a-1, 0, 0.85*i01, 0.85*i02, -2, xc-0.75, 0];
%             upper = [a+1, tol, 1.15*abs(i01), 1.15*abs(i02), 5, xc+0.75, 2*y0];
%         end

        fo = fitoptions('Method','NonlinearLeastSquares','StartPoint',start,'Lower',lower,'Upper',upper);

        ft = fittype('y0 + (i01/(fwhm * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(x-x1)^2/(fwhm^2/log(4))) + (i02/(fwhm * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(x-x2)^2/(fwhm^2/log(4))) + (i03/(fwhm * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(x-x3)^2/(fwhm^2/log(4)))',...
                        'options',fo,'problem',{'y0','i01','i02','i03','fwhm'});
                    
        % load desired depth file and update main grasp GUI
        % the +1 is necessary as file 1 represents a sum
        status_flags.selector.fd = img_num+1;
        grasp_update;

        % plot I vs xi and create the current figure, store handle
        radial_average_callbacks('averaging_control','azimuthal');
        mf_fitter.handles.plot_handle = gcf;
        grasp_update;
            
        angles = plot_info.export_data(:,1);
        intensity = plot_info.export_data(:,2);
        intensity_err = plot_info.export_data(:,3);
       
        [fitobject, gof] = fit(angles,intensity,ft,'problem',{y0,i01,i02,i03,fwhm})
        ci = confint(fitobject)
        plot(fitobject,angles,intensity);
        
        mf_fitter.data.angle = angles;
        mf_fitter.data.intensity(:,img_num) = intensity;
        mf_fitter.data.intensity_err(:,img_num) = intensity_err;
        
%         mf_fitter.fit_data.background(img_num,1) = fitobject.y0;
%         mf_fitter.fit_data.intensity1(img_num,1) = fitobject.i01;
%         mf_fitter.fit_data.intensity2(img_num,1) = fitobject.i02;
%         mf_fitter.fit_data.intensity3(img_num,1) = fitobject.i03;
        mf_fitter.fit_data.center1(img_num,1) = fitobject.x1
        mf_fitter.fit_data.center1(img_num,2) = ((fitobject.x1 - ci(1)) + (ci(2) - fitobject.x1))/2
        mf_fitter.fit_data.center1_err(img_num,1) = ci(1)
        mf_fitter.fit_data.center1_err(img_num,2) = ci(2)
        mf_fitter.fit_data.center2(img_num,1) = fitobject.x2;
        mf_fitter.fit_data.center2_err(img_num,1) = ci(3);
        mf_fitter.fit_data.center2_err(img_num,2) = ci(4);
        mf_fitter.fit_data.center2(img_num,2) = ((fitobject.x2 - ci(3)) + (ci(4) - fitobject.x2))/2
        mf_fitter.fit_data.center3(img_num,1) = fitobject.x3
        mf_fitter.fit_data.center3_err(img_num,1) = ci(5)
        mf_fitter.fit_data.center3_err(img_num,2) = ci(6) 
        mf_fitter.fit_data.center3(img_num,2) = ((fitobject.x3 - ci(5)) + (ci(6) - fitobject.x3))/2
     
        case 'matlab_fit_fixed'
        % specify fit options and type
        % order of coeffecients for boundaries:
        %   a, b, i01, i02, i03, xc, y0

        %tol = str2num(get(mf_fitter.handles.tolerance,'String'));
        delta12 = mf_fitter.fit_data.center2(img_num) - mf_fitter.fit_data.center1(img_num);
        delta23 = mf_fitter.fit_data.center3(img_num) - mf_fitter.fit_data.center2(img_num);
        
        a = (delta12+delta23)/2;
        %b = delta23-delta12;
        fwhm = 0.9*mf_fitter.averages.fwhm;
        i01 = mf_fitter.fit_data.intensity1(img_num);
        i02 = mf_fitter.fit_data.intensity2(img_num);
        i03 = mf_fitter.fit_data.intensity3(img_num);
        xc = mf_fitter.fit_data.center2(img_num);
        y0 = mf_fitter.fit_data.background(img_num);

%         start = [a, b, abs(i01), abs(i02), abs(i03), xc, y0];
%         lower = [a-1, 0, 0.85*i01, 0.85*i02, 0.85*i03, xc-0.75, 0];
%         upper = [a+1, tol, 1.15*abs(i01), 1.15*abs(i02), 1.15*abs(i03), xc+0.75, 4];

        start = [a, abs(i01), abs(i02), abs(i03), xc, y0];
        lower = [a-1, 0.85*i01, 0.85*i02, 0.85*i03, xc-0.75, 0];
        upper = [a+1, 1.15*abs(i01), 1.15*abs(i02), 1.15*abs(i03), xc+0.75, 4];

        if( abs(0.3*i01) < 2*10^(-8) )
            start = [a, abs(i01), abs(i02), abs(i03), xc, y0];
            lower = [a-1, -2, 0.85*i02, 0.85*i03, xc-0.75, 0];
            upper = [a+1, 5, 1.15*abs(i02), 1.15*abs(i03), xc+0.75, 2*y0]; 
        elseif( abs(0.3*i03) < 2*10^(-8) )
            start = [a, abs(i01), abs(i02), abs(i03), xc, y0];
            lower = [a-1, 0.85*i01, 0.85*i02, -2, xc-0.75, 0];
            upper = [a+1, 1.15*abs(i01), 1.15*abs(i02), 5, xc+0.75, 2*y0];
        end

        fo = fitoptions('Method','NonlinearLeastSquares','StartPoint',start,'Lower',lower,'Upper',upper);

        ft = fittype('y0 + (i01/(fwhm * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(x-(xc-a))^2/(fwhm^2/log(4))) + (i02/(fwhm * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(x-xc)^2/(fwhm^2/log(4))) + (i03/(fwhm * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(x-(xc+a))^2/(fwhm^2/log(4)))',...
                        'options',fo,'problem',{'fwhm'});        
        
%         ft = fittype('y0 + (i01/(fwhm * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(x-(xc-a))^2/(fwhm^2/log(4))) + (i02/(fwhm * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(x-xc)^2/(fwhm^2/log(4))) + (i03/(fwhm * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(x-(xc+a+b))^2/(fwhm^2/log(4)))',...
%                         'options',fo,'problem',{'fwhm'});
                    
        % load desired depth file and update main grasp GUI
        % the +1 is necessary as file 1 represents a sum
        status_flags.selector.fd = img_num+1;
        grasp_update;

        % plot I vs xi and create the current figure, store handle
        radial_average_callbacks('averaging_control','azimuthal');
        mf_fitter.handles.plot_handle = gcf;
        grasp_update;
            
        angles = plot_info.export_data(:,1);
        intensity = plot_info.export_data(:,2);
        intensity_err = plot_info.export_data(:,3);
       
        [fitobject, gof] = fit(angles,intensity,ft,'problem',fwhm);
        ci = confint(fitobject);
        plot(fitobject,angles,intensity);
        
        mf_fitter.data.angle = angles;
        mf_fitter.data.intensity(:,img_num) = intensity;
        mf_fitter.data.intensity_err(:,img_num) = intensity_err;
        
        
        % To put back in tolerance, remember to switch the ci(index)
        mf_fitter.fit_data.background(img_num,1) = fitobject.y0;
        mf_fitter.fit_data.intensity1(img_num,1) = fitobject.i01;
        mf_fitter.fit_data.intensity2(img_num,1) = fitobject.i02;
        mf_fitter.fit_data.intensity3(img_num,1) = fitobject.i03;
        mf_fitter.fit_data.center1(img_num,1) = fitobject.xc - fitobject.a;  
        mf_fitter.fit_data.center1(img_num,2) = sqrt( (((fitobject.xc - ci(9)) + (ci(10) - fitobject.xc))/2)^2 + (((fitobject.a - ci(1)) + (ci(2) - fitobject.a))/2)^2 );
        mf_fitter.fit_data.center1_err(img_num,1) = ci(9) + ci(1);
        mf_fitter.fit_data.center1_err(img_num,2) = ci(10) + ci(2);
        mf_fitter.fit_data.center2(img_num,1) = fitobject.xc;
        mf_fitter.fit_data.center2_err(img_num,1) = ci(9);
        mf_fitter.fit_data.center2_err(img_num,2) = ci(10);
        mf_fitter.fit_data.center2(img_num,2) = sqrt( (((fitobject.xc-ci(11)) + (ci(12)-fitobject.xc))/2)^2 + (((fitobject.a-ci(1)) + (ci(2)-fitobject.a))/2)^2); % + (((fitobject.b-ci(3))+(ci(4)-fitobject.b))/2)^2 );
        mf_fitter.fit_data.center3(img_num,1) = fitobject.xc + fitobject.a; % + fitobject.b;
        mf_fitter.fit_data.center3_err(img_num,1) = ci(9) + ci(1); % + ci(3);
        mf_fitter.fit_data.center3_err(img_num,2) = ci(10) + ci(2); % + ci(4);
        mf_fitter.fit_data.center3(img_num,2) = mf_fitter.fit_data.center1(img_num,2);

        
            case 'matlab_fit_free'
        % specify fit options and type
        % order of coeffecients for boundaries:
        %   a, b, fwhm, i01, i02, i03, xc, y0
        
        inner_num = mf_fitter.good_peaks.inner_num;
        outer_num = mf_fitter.good_peaks.outers_num;
        delta12 = mf_fitter.fit_data.center2(inner_num) - mf_fitter.fit_data.center1(outer_num);
        delta23 = mf_fitter.fit_data.center3(outer_num) - mf_fitter.fit_data.center2(inner_num);
        
        a = (delta12+delta23)/2  %(img_num  / mf_fitter.depth) * 
        %b = str2num(get(mf_fitter.handles.tolerance,'String'));
        fwhm = mf_fitter.averages.fwhm;
        i01 = mf_fitter.fit_data.intensity1(outer_num);
        i02 = mf_fitter.fit_data.intensity2(inner_num);
        i03 = mf_fitter.fit_data.intensity3(outer_num);
        xc = mf_fitter.fit_data.center2(mf_fitter.good_peaks.inner_num);
        y0 = mf_fitter.averages.background;
    
        % In matlab fit codes, commented out sections relate to tol.
        %start = [a, 0, fwhm, i01, i02, i03, xc, y0]
        %lower = [3, 0, fwhm-1, 0, 1, 0, xc-1, 0]
        %upper = [a+2, b, fwhm+1, 1.5*i01, 1.5*i02, 1.5*i03, xc+1, 2*y0]
            
        start = [a, fwhm, i01, i02, i03, xc, y0]
        lower = [3, fwhm-1, 0, 1, 0, xc-1, 0]
        upper = [a+2, fwhm+1, 1.5*i01, 1.5*i02, 1.5*i03, xc+1, 2*y0]
        
       % if(a<3) 
        %    start = [3, fwhm, 10, i02, 10, xc, y0]
        %    lower = [2.5, fwhm-1, -5, 0.5*i02, -5, xc-1, 0]
        %    upper = [5, fwhm+1, 10, 1.5*i02, 10, xc+1, 2*y0]
        %end

        
        fo = fitoptions('Method','NonlinearLeastSquares','StartPoint',start,'Lower',lower,'Upper',upper);
        
        ft = fittype('y0 + (i01/(fwhm * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(x-(xc-a))^2/(fwhm^2/log(4))) + (i02/(fwhm * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(x-xc)^2/(fwhm^2/log(4))) + (i03/(fwhm * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(x-(xc+a))^2/(fwhm^2/log(4)))',...
                        'options',fo);

        % load desired depth file and update main grasp GUI
        % the +1 is necessary as file 1 represents a sum
        status_flags.selector.fd = img_num+1;
        grasp_update;

        % plot I vs xi and create the current figure, store handle
        radial_average_callbacks('averaging_control','azimuthal');
        mf_fitter.handles.plot_handle = gcf;
        grasp_update;
            
        angles = plot_info.export_data(:,1);
        intensity = plot_info.export_data(:,2);
        intensity_err = plot_info.export_data(:,3);
                    
        [fitobject, gof] = fit(angles,intensity,ft);
        ci = confint(fitobject);
        
        mf_fitter.fit_data.background(img_num,1) = fitobject.y0;
        mf_fitter.fit_data.intensity1(img_num,1) = fitobject.i01;
        mf_fitter.fit_data.intensity2(img_num,1) = fitobject.i02;
        mf_fitter.fit_data.intensity3(img_num,1) = fitobject.i03;
        mf_fitter.fit_data.center1(img_num,1) = fitobject.xc - fitobject.a;
        mf_fitter.fit_data.center2(img_num,1) = fitobject.xc;
        mf_fitter.fit_data.center3(img_num,1) = fitobject.xc + fitobject.a; % + fitobject.b;
        mf_fitter.fit_data.fwhm(img_num,1) = fitobject.fwhm;
  
        
           case 'set_free'
        % sets curve fit window
        % fixes neither the fwhm nor the centers
        % groups background, fwhm
        % runs on the reference peaks (one or two peak fit)

        mf_fitter_callbacks('set_params',num_peaks,img_num);
        
        status_flags.fitter.function_info_1d.group = repmat([1,0,0,1],1,num_peaks);
        status_flags.fitter.function_info_1d.err_values = repmat([0,0,0,0],1,num_peaks);
        
        if( num_peaks==1 )
            % fitter capable of autoguessing
            status_flags.fitter.function_info_1d.fix = [0,0,0,0];
            status_flags.fitter.function_info_1d.values = [0,0,0,0];
        elseif (num_peaks==2)
            % grab guess values from 1st reference file
            % since this was a one peak fit, the data is stored under peak2
            % only need to set background and fwhm once since grouped
            % order: [yo, io, xc, fwhm]
            status_flags.fitter.function_info_1d.fix = [0,0,0,0,1,0,0,1];
            status_flags.fitter.function_info_1d.values = repmat([0],1,8);
            status_flags.fitter.function_info_1d.values(1) = mf_fitter.fit_data.background(mf_fitter.good_peaks.inner_num,1);
            status_flags.fitter.function_info_1d.values(2) = mf_fitter.fit_data.intensity2(mf_fitter.good_peaks.inner_num,1);
            status_flags.fitter.function_info_1d.values(3) = mf_fitter.fit_data.center2(mf_fitter.good_peaks.inner_num,1) - 5;
            status_flags.fitter.function_info_1d.values(4) = mf_fitter.fit_data.fwhm(mf_fitter.good_peaks.inner_num,1);
            status_flags.fitter.function_info_1d.values(6) = mf_fitter.fit_data.intensity2(mf_fitter.good_peaks.inner_num,1);
            status_flags.fitter.function_info_1d.values(7) = mf_fitter.fit_data.center2(mf_fitter.good_peaks.inner_num,1)+ 5;
        else
            status_flags.fitter.function_info_1d.fix = [0,0,0,0,1,0,0,1,1,0,0,1];
            status_flags.fitter.function_info_1d.values = repmat([0],1,12);
            status_flags.fitter.function_info_1d.values(1) = mf_fitter.fit_data.background(mf_fitter.good_peaks.inner_num,1);
            status_flags.fitter.function_info_1d.values(2) = mf_fitter.fit_data.intensity2(mf_fitter.good_peaks.inner_num,1);
            status_flags.fitter.function_info_1d.values(3) = mf_fitter.fit_data.center2(mf_fitter.good_peaks.inner_num,1) - 5;
            status_flags.fitter.function_info_1d.values(4) = mf_fitter.fit_data.fwhm(mf_fitter.good_peaks.inner_num,1);
            status_flags.fitter.function_info_1d.values(6) = mf_fitter.fit_data.intensity2(mf_fitter.good_peaks.inner_num,1);
            status_flags.fitter.function_info_1d.values(7) = mf_fitter.fit_data.center2(mf_fitter.good_peaks.inner_num,1);
            status_flags.fitter.function_info_1d.values(10) = mf_fitter.fit_data.intensity2(mf_fitter.good_peaks.inner_num,1);
            status_flags.fitter.function_info_1d.values(11) = mf_fitter.fit_data.center2(mf_fitter.good_peaks.inner_num,1)+ 5;
                
        end

        
       case 'set_free2'
        % sets curve fit window
        % fixes neither the fwhm nor the centers
        % groups background, fwhm
        % runs on the reference peaks (one or two peak fit)

        mf_fitter_callbacks('set_params',num_peaks,img_num);
        
        status_flags.fitter.function_info_1d.group = repmat([1,0,0,1],1,num_peaks);
        status_flags.fitter.function_info_1d.err_values = repmat([0,0,0,0],1,num_peaks);
        
        if( num_peaks==1 )
            % fitter capable of autoguessing
            status_flags.fitter.function_info_1d.fix = [0,0,0,0];
            status_flags.fitter.function_info_1d.values = [0,0,0,0];
        elseif( num_peaks==2)
            % grab guess values from 1st reference file
            % since this was a one peak fit, the data is stored under peak2
            % only need to set background and fwhm once since grouped
            % order: [yo, io, xc, fwhm]
            status_flags.fitter.function_info_1d.fix = [0,0,0,0,1,0,0,1];
            status_flags.fitter.function_info_1d.values = repmat([0],1,8);
            status_flags.fitter.function_info_1d.values(1) = mf_fitter.fit_data.background(mf_fitter.good_peaks.inner_num,1);
            status_flags.fitter.function_info_1d.values(2) = mf_fitter.fit_data.intensity2(mf_fitter.good_peaks.inner_num,1);
            status_flags.fitter.function_info_1d.values(3) = mf_fitter.fit_data.center2(mf_fitter.good_peaks.inner_num,1) - 5;
            status_flags.fitter.function_info_1d.values(4) = mf_fitter.fit_data.fwhm(mf_fitter.good_peaks.inner_num,1);
            status_flags.fitter.function_info_1d.values(6) = mf_fitter.fit_data.intensity2(mf_fitter.good_peaks.inner_num,1);
            status_flags.fitter.function_info_1d.values(7) = mf_fitter.fit_data.center2(mf_fitter.good_peaks.inner_num,1)+ 5;
        else
             status_flags.fitter.function_info_1d.fix = [0,0,0,0,1,0,0,1,1,0,0,1];
            status_flags.fitter.function_info_1d.values = repmat([0],1,12);
            status_flags.fitter.function_info_1d.values(1) = mf_fitter.fit_data.background(mf_fitter.good_peaks.inner_num,1);
            status_flags.fitter.function_info_1d.values(2) = mf_fitter.fit_data.intensity1(mf_fitter.good_peaks.outers_num,1);
            status_flags.fitter.function_info_1d.values(3) = mf_fitter.fit_data.center1(mf_fitter.good_peaks.outers_num,1);
            status_flags.fitter.function_info_1d.values(4) = mf_fitter.fit_data.fwhm(mf_fitter.good_peaks.inner_num,1);
            status_flags.fitter.function_info_1d.values(6) = mf_fitter.fit_data.intensity2(mf_fitter.good_peaks.inner_num,1);
            status_flags.fitter.function_info_1d.values(7) = mf_fitter.fit_data.center2(mf_fitter.good_peaks.inner_num,1);
            status_flags.fitter.function_info_1d.values(10) = mf_fitter.fit_data.intensity3(mf_fitter.good_peaks.outers_num,1);
            status_flags.fitter.function_info_1d.values(11) = mf_fitter.fit_data.center3(mf_fitter.good_peaks.outers_num,1);
        end
        
        
    case 'set_centers'
        % sets curve fit window
        % fixes the centers, but not the fwhm
        % groups background, fwhm
        % prepares for 3-peak fit
        
        mf_fitter_callbacks('set_params',num_peaks,img_num);
        
        status_flags.fitter.function_info_1d.group = repmat([1,0,0,1],1,num_peaks);
        status_flags.fitter.function_info_1d.fix = [0,0,1,0,1,0,1,1,1,0,1,1];
        status_flags.fitter.function_info_1d.err_values = repmat([0,0,0,0],1,num_peaks);
        status_flags.fitter.function_info_1d.values = repmat([0],1,12);
        
        % centers are being fixed to values determined from reference files
        status_flags.fitter.function_info_1d.values(3) = mf_fitter.fit_data.center1(mf_fitter.good_peaks.outers_num,1);
        status_flags.fitter.function_info_1d.values(7) = mf_fitter.fit_data.center2(mf_fitter.good_peaks.inner_num,1);
        status_flags.fitter.function_info_1d.values(11) = mf_fitter.fit_data.center3(mf_fitter.good_peaks.outers_num,1);
        
        % fwhm and background guess values come from averages
        status_flags.fitter.function_info_1d.values(1) = mf_fitter.averages.background;
        status_flags.fitter.function_info_1d.values(4) = mf_fitter.averages.fwhm;
        
        % grab guess values from previous files if fit performed
        % otherwise use reference files (also for img_num 1)
        % only need to set background and fwhm once since grouped
        % order: [yo, io, xc, fwhm]
        if( img_num==1 || mf_fitter.did_fit(img_num-1)==0 )
            status_flags.fitter.function_info_1d.values(2) = mf_fitter.fit_data.intensity2(mf_fitter.good_peaks.inner_num,1);
            status_flags.fitter.function_info_1d.values(6) = mf_fitter.fit_data.intensity2(mf_fitter.good_peaks.inner_num,1);
            status_flags.fitter.function_info_1d.values(10) = mf_fitter.fit_data.intensity2(mf_fitter.good_peaks.inner_num,1);
            
        else
            status_flags.fitter.function_info_1d.values(2) = mf_fitter.fit_data.intensity1(img_num-1,1);
            status_flags.fitter.function_info_1d.values(6) = mf_fitter.fit_data.intensity2(img_num-1,1);
            status_flags.fitter.function_info_1d.values(10) = mf_fitter.fit_data.intensity3(img_num-1,1);
            
        end
        
            
    case 'set_fwhm'
        % sets curve fit window
        % fixes the fwhm, but not the centers
        % groups background, fwhm
        % prepares for 3-peak fit
        
        mf_fitter_callbacks('set_params',num_peaks,img_num);
        
        status_flags.fitter.function_info_1d.group = repmat([1,0,0,1],1,num_peaks);
        status_flags.fitter.function_info_1d.err_values = repmat([0,0,0,0],1,num_peaks);
        status_flags.fitter.function_info_1d.values = repmat([0],1,4*num_peaks);
        
        if (num_peaks==1)
            status_flags.fitter.function_info_1d.fix = [0,0,0,1];
            % centers
            status_flags.fitter.function_info_1d.values(3) = mf_fitter.fit_data.center2(mf_fitter.good_peaks.outers_num,1);
            
        elseif (num_peaks==2)
            status_flags.fitter.function_info_1d.fix = [0,0,0,1,1,0,0,1];
            % centers
            status_flags.fitter.function_info_1d.values(3) = mf_fitter.fit_data.center1(mf_fitter.good_peaks.outers_num,1);
            status_flags.fitter.function_info_1d.values(7) = mf_fitter.fit_data.center3(mf_fitter.good_peaks.inner_num,1);
        else
            status_flags.fitter.function_info_1d.fix = [0,0,0,1,1,0,0,1,1,0,0,1];
            % centers
            status_flags.fitter.function_info_1d.values(3) = mf_fitter.fit_data.center1(mf_fitter.good_peaks.outers_num,1);
            status_flags.fitter.function_info_1d.values(7) = mf_fitter.fit_data.center2(mf_fitter.good_peaks.inner_num,1);
            status_flags.fitter.function_info_1d.values(11) = mf_fitter.fit_data.center3(mf_fitter.good_peaks.outers_num,1);
        end
              
        % fwhm fix values and background guess values come from averages
        status_flags.fitter.function_info_1d.values(1) = mf_fitter.averages.background;
        status_flags.fitter.function_info_1d.values(4) = mf_fitter.averages.fwhm;
        
        % grab guess values from previous files if fit performed
        % otherwise use reference files (also for img_num 1)
        % only need to set background and fwhm once since grouped
        % order: [yo, io, xc, fwhm]
        if( img_num==1 || mf_fitter.did_fit(img_num-1)==0 )
            status_flags.fitter.function_info_1d.values(2) = mf_fitter.fit_data.intensity2(mf_fitter.good_peaks.inner_num,1);
            if( num_peaks~=1 )
                status_flags.fitter.function_info_1d.values(6) = mf_fitter.fit_data.intensity2(mf_fitter.good_peaks.inner_num,1);
                if( num_peaks~=2 )
                    status_flags.fitter.function_info_1d.values(10) = mf_fitter.fit_data.intensity2(mf_fitter.good_peaks.inner_num,1);
                end
            end
        else
            if(num_peaks==1)
                status_flags.fitter.function_info_1d.values(2) = mf_fitter.fit_data.intensity2(img_num-1,1);
            elseif(num_peaks==2)
                status_flags.fitter.function_info_1d.values(2) = mf_fitter.fit_data.intensity1(img_num-1,1);
                status_flags.fitter.function_info_1d.values(6) = mf_fitter.fit_data.intensity3(img_num-1,1);
            else
                status_flags.fitter.function_info_1d.values(2) = mf_fitter.fit_data.intensity1(img_num-1,1);
                status_flags.fitter.function_info_1d.values(6) = mf_fitter.fit_data.intensity2(img_num-1,1);
                status_flags.fitter.function_info_1d.values(10) = mf_fitter.fit_data.intensity3(img_num-1,1);
            end
        end
        

    case 'set_both'
        % sets curve fit window
        % fixes both the fwhm and the centers
        % groups background, fwhm
        % prepares for 3-peak fit
        
        mf_fitter_callbacks('set_params',num_peaks,img_num);
        
        status_flags.fitter.function_info_1d.group = repmat([1,0,0,1],1,num_peaks);
        status_flags.fitter.function_info_1d.fix = [0,0,1,1,1,0,1,1,1,0,1,1]
        status_flags.fitter.function_info_1d.err_values = repmat([0,0,0,0],1,num_peaks);
        status_flags.fitter.function_info_1d.values = repmat([0],1,12);
        
        % centers are being fixed as the average value from previous fit
        status_flags.fitter.function_info_1d.values(3) = mf_fitter.averages.center1;
        status_flags.fitter.function_info_1d.values(7) = mf_fitter.averages.center2;
        status_flags.fitter.function_info_1d.values(11) = mf_fitter.averages.center3;
        
        % fwhm fix values and background guess values come from averages
        status_flags.fitter.function_info_1d.values(1) = mf_fitter.averages.background;
        status_flags.fitter.function_info_1d.values(4) = mf_fitter.averages.fwhm;
        
        % grab guess values from previous files if fit performed
        % otherwise use reference files (also for img_num 1)
        % only need to set background and fwhm once since grouped
        % order: [yo, io, xc, fwhm]
        if( img_num==1 || mf_fitter.did_fit(img_num-1)==0 )
            status_flags.fitter.function_info_1d.values(2) = mf_fitter.fit_data.intensity2(mf_fitter.good_peaks.inner_num,1);
            status_flags.fitter.function_info_1d.values(6) = mf_fitter.fit_data.intensity2(mf_fitter.good_peaks.inner_num,1);
            status_flags.fitter.function_info_1d.values(10) = mf_fitter.fit_data.intensity2(mf_fitter.good_peaks.inner_num,1);
            
        else
            status_flags.fitter.function_info_1d.values(2) = mf_fitter.fit_data.intensity1(img_num-1,1);
            status_flags.fitter.function_info_1d.values(6) = mf_fitter.fit_data.intensity2(img_num-1,1);
            status_flags.fitter.function_info_1d.values(10) = mf_fitter.fit_data.intensity3(img_num-1,1);
            
        end
           

    case 'set_fwhm_mvc'
        % sets curve fit window
        % fixes the fwhm, frees GS centers, fixes MS to average
        % groups background, fwhm
        % prepares for 3-peak fit
        
        mf_fitter_callbacks('set_params',num_peaks,img_num);
        
        status_flags.fitter.function_info_1d.group = repmat([1,0,0,1],1,num_peaks);
        status_flags.fitter.function_info_1d.err_values = repmat([0,0,0,0],1,num_peaks);
        status_flags.fitter.function_info_1d.values = repmat([0],1,4*num_peaks);
        
        if (num_peaks==1)
            status_flags.fitter.function_info_1d.fix = [0,0,1,1];
            % centers
            status_flags.fitter.function_info_1d.values(3) = mf_fitter.averages.center2;
            
        elseif (num_peaks==2)
            status_flags.fitter.function_info_1d.fix = [0,0,0,1,1,0,0,1];
            % centers
            status_flags.fitter.function_info_1d.values(3) = mf_fitter.fit_data.center1(img_num-1,1);
            status_flags.fitter.function_info_1d.values(7) = mf_fitter.fit_data.center3(img_num-1,1);
        else
            status_flags.fitter.function_info_1d.fix = [0,0,0,1,1,0,1,1,1,0,0,1];
            % centers
            status_flags.fitter.function_info_1d.values(3) = mf_fitter.fit_data.center1(img_num-1,1);
            status_flags.fitter.function_info_1d.values(7) =  mf_fitter.averages.center2;
            status_flags.fitter.function_info_1d.values(11) = mf_fitter.fit_data.center3(img_num-1,1);
        end
              
        % fwhm fix values and background guess values come from averages
        status_flags.fitter.function_info_1d.values(1) = mf_fitter.averages.background;
        status_flags.fitter.function_info_1d.values(4) = mf_fitter.averages.fwhm; 
        
        % grab guess values from previous files if fit performed
        % otherwise use reference files (also for img_num 1)
        % only need to set background and fwhm once since grouped
        % order: [yo, io, xc, fwhm]
        if( img_num==1 || mf_fitter.did_fit(img_num-1)==0 )
            status_flags.fitter.function_info_1d.values(2) = mf_fitter.fit_data.intensity2(mf_fitter.good_peaks.inner_num,1);
            if( num_peaks~=1 )
                status_flags.fitter.function_info_1d.values(6) = mf_fitter.fit_data.intensity2(mf_fitter.good_peaks.inner_num,1);
                if( num_peaks~=2 )
                    status_flags.fitter.function_info_1d.values(10) = mf_fitter.fit_data.intensity2(mf_fitter.good_peaks.inner_num,1);
                end
            end
        else
            if(num_peaks==1)
                % work with previous file number since this callback is to
                % fix a bad fit
                status_flags.fitter.function_info_1d.values(2) = mf_fitter.fit_data.intensity2(img_num-1,1);  
            elseif(num_peaks==2)
                status_flags.fitter.function_info_1d.values(2) = mf_fitter.fit_data.intensity1(img_num-1,1);
                status_flags.fitter.function_info_1d.values(6) = mf_fitter.fit_data.intensity3(img_num-1,1);
            else
                status_flags.fitter.function_info_1d.values(2) = mf_fitter.fit_data.intensity1(img_num-1,1);
                status_flags.fitter.function_info_1d.values(6) = mf_fitter.fit_data.intensity2(img_num-1,1);
                status_flags.fitter.function_info_1d.values(10) = mf_fitter.fit_data.intensity3(img_num-1,1);
            end
        end
        
        
    case 'set_centers_mvc'
        % sets curve fit window
        % fixes the MS center to the average, fixes the GS centers, and frees the fwhm
        % groups background, fwhm
        % prepares for 3-peak fit
        
        mf_fitter_callbacks('set_params',num_peaks,img_num);
        
        status_flags.fitter.function_info_1d.group = repmat([1,0,0,1],1,num_peaks);
        status_flags.fitter.function_info_1d.err_values = repmat([0,0,0,0],1,num_peaks);
        status_flags.fitter.function_info_1d.values = repmat([0],1,4*num_peaks);
        
        if (num_peaks==1)
            status_flags.fitter.function_info_1d.fix = [0,0,1,0];
            % centers
            status_flags.fitter.function_info_1d.values(3) = mf_fitter.averages.center2;
            
        elseif (num_peaks==2)
            status_flags.fitter.function_info_1d.fix = [0,0,1,0,1,0,1,1];
            % centers
            status_flags.fitter.function_info_1d.values(3) = mf_fitter.fit_data.center1(img_num,1);
            status_flags.fitter.function_info_1d.values(7) = mf_fitter.fit_data.center3(img_num,1);
        else
            status_flags.fitter.function_info_1d.fix = [0,0,1,0,1,0,1,1,1,0,1,1];
            % centers
            status_flags.fitter.function_info_1d.values(3) = mf_fitter.fit_data.center1(img_num,1);
            status_flags.fitter.function_info_1d.values(7) = mf_fitter.averages.center2;
            status_flags.fitter.function_info_1d.values(11) = mf_fitter.fit_data.center3(img_num,1);
        end
        
        % grab guess values from previous files if fit performed
        % otherwise use reference files (also for img_num 1)
        % only need to set background and fwhm once since grouped
        % order: [yo, io, xc, fwhm]
        if( img_num==1 || mf_fitter.did_fit(img_num-1)==0 )
            status_flags.fitter.function_info_1d.values(2) = mf_fitter.fit_data.intensity2(img_num,1);
            if( num_peaks~=1 )
                status_flags.fitter.function_info_1d.values(6) = mf_fitter.fit_data.intensity2(img_num,1);
                if( num_peaks~=2 )
                    status_flags.fitter.function_info_1d.values(10) = mf_fitter.fit_data.intensity2(img_num,1);
                end
            end
        else
            if(num_peaks==1)
                status_flags.fitter.function_info_1d.values(2) = mf_fitter.fit_data.intensity2(img_num,1);
            elseif(num_peaks==2)
                status_flags.fitter.function_info_1d.values(2) = mf_fitter.fit_data.intensity1(img_num,1);
                status_flags.fitter.function_info_1d.values(6) = mf_fitter.fit_data.intensity3(img_num,1);
            else
                status_flags.fitter.function_info_1d.values(2) = mf_fitter.fit_data.intensity1(img_num,1);
                status_flags.fitter.function_info_1d.values(6) = mf_fitter.fit_data.intensity2(img_num,1);
                status_flags.fitter.function_info_1d.values(10) = mf_fitter.fit_data.intensity3(img_num,1);
            end
        end
        
        % MS center fix values and fwhm, background guess values come from averages
        status_flags.fitter.function_info_1d.values(1) = mf_fitter.averages.background;
        status_flags.fitter.function_info_1d.values(4) = mf_fitter.averages.fwhm;
        

    case 'set_fwhm&centers_mvc'
        % sets curve fit window
        % fixes the centers and the fwhm
        % groups background, fwhm
        % prepares for 3-peak fit
        % determines the final intensities!
      
        mf_fitter_callbacks('set_params',num_peaks,img_num);
        
        status_flags.fitter.function_info_1d.group = repmat([1,0,0,1],1,num_peaks);
        status_flags.fitter.function_info_1d.err_values = repmat([0,0,0,0],1,num_peaks);
        status_flags.fitter.function_info_1d.values = repmat([0],1,4*num_peaks);
        
        if (num_peaks==1)
            status_flags.fitter.function_info_1d.fix = [0,0,1,1];
            status_flags.fitter.function_info_1d.values(2) = mf_fitter.fit_data.intensity2(img_num,1);
            status_flags.fitter.function_info_1d.values(3) = mf_fitter.fit_data.center3(img_num,1);
            
        elseif (num_peaks==2)
            status_flags.fitter.function_info_1d.fix = [0,0,1,1,1,0,1,1];
            status_flags.fitter.function_info_1d.values(2) = mf_fitter.fit_data.intensity1(img_num,1);
            status_flags.fitter.function_info_1d.values(6) = mf_fitter.fit_data.intensity3(img_num,1);
            status_flags.fitter.function_info_1d.values(3) = mf_fitter.fit_data.center1(img_num,1);
            status_flags.fitter.function_info_1d.values(7) = mf_fitter.fit_data.center3(img_num,1);
        else
            status_flags.fitter.function_info_1d.fix = [0,0,1,1,1,0,1,1,1,0,1,1];
            status_flags.fitter.function_info_1d.values(2) = mf_fitter.fit_data.intensity1(img_num,1);
            status_flags.fitter.function_info_1d.values(6) = mf_fitter.fit_data.intensity2(img_num,1);
            status_flags.fitter.function_info_1d.values(10) = mf_fitter.fit_data.intensity3(img_num,1);
            status_flags.fitter.function_info_1d.values(3) = mf_fitter.fit_data.center1(img_num,1);
            status_flags.fitter.function_info_1d.values(7) = mf_fitter.fit_data.center2(img_num,1);
            status_flags.fitter.function_info_1d.values(11) = mf_fitter.fit_data.center3(img_num,1);
        end

        %fwhm and background guess values come from averages
        status_flags.fitter.function_info_1d.values(1) = mf_fitter.averages.background;
        status_flags.fitter.function_info_1d.values(4) = mf_fitter.averages.fwhm; 
        
end 


end

