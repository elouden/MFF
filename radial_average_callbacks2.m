function radial_average_callbacks2(to_do,options,options2)

%Does something different I think?

% v 7.4
% 7/22/14 MFF Allan

global displayimage
global status_flags
global grasp_handles
global inst_params
global grasp_data
global grasp_env
global tof_iq_store


if nargin <3; options2 = []; end
if nargin <2; options = ''; end
if nargin <1; to_do = ''; end


%Radial Average Window & Parameter Callbacks
switch to_do
    
    case 'tof_iq'
        status_flags.analysis_modules.radial_average.tof_iq = get(gcbo,'value');
        
%     case 'depth_frame_start'
%         temp = str2num(get(gcbo,'string'));
%         if not(isempty(temp));
%             status_flags.analysis_modules.radial_average.depth_frame_start = temp;
%         end
%         
%     case 'depth_frame_end'
%         temp = str2num(get(gcbo,'string'));
%         if not(isempty(temp));
%             status_flags.analysis_modules.radial_average.depth_frame_end = temp;
%         end
%         
        
    
    case 'd33_tof_combine_check'
        status_flags.analysis_modules.radial_average.d33_tof_combine = get(gcbo,'value');
        if status_flags.analysis_modules.radial_average.d33_tof_combine ==1;
            %open rebin window
            rebin_window;
        end
        
    case 'close'
        grasp_handles.window_modules.radial_average.window = [];
        return
    case 'q_bin_pixels'
        status_flags.analysis_modules.radial_average.q_bin_units = 'pixels';
        
    case 'q_bin_absolute'
        status_flags.analysis_modules.radial_average.q_bin_units = 'absolute';
        status_flags.analysis_modules.radial_average.q_bin_absolute_scale = options;
        
    case 'q_bin_resolution'
        status_flags.analysis_modules.radial_average.q_bin_units = 'resolution';
        
    case 'theta_bin_pixels'
        status_flags.analysis_modules.radial_average.theta_bin_units = 'pixels';
        
    case 'theta_bin_absolute'
        status_flags.analysis_modules.radial_average.theta_bin_units = 'absolute';
        status_flags.analysis_modules.radial_average.theta_bin_absolute_scale = options;
        
    case 'theta_bin_resolution'
        status_flags.analysis_modules.radial_average.theta_bin_units = 'resolution';
        
    case 'azimuth_bin_absolute'
        status_flags.analysis_modules.radial_average.azimuth_bin_units = 'absolute';
        
    case 'azimuth_bin_smart'
        status_flags.analysis_modules.radial_average.azimuth_bin_units = 'smart';
        
    case 'q_bin'
        temp = str2num(get(grasp_handles.window_modules.radial_average.q_bin,'string'));
        if not(isempty(temp));
            if strcmp(status_flags.analysis_modules.radial_average.q_bin_units,'pixels')
                status_flags.analysis_modules.radial_average.q_bin_pixels = temp;
            elseif strcmp(status_flags.analysis_modules.radial_average.q_bin_units,'absolute')
                status_flags.analysis_modules.radial_average.q_bin_absolute = temp;
            elseif strcmp(status_flags.analysis_modules.radial_average.q_bin_units,'resolution')
                status_flags.analysis_modules.radial_average.q_bin_resolution = temp;
            end
        end
        
    case 'theta_bin'
        temp = str2num(get(grasp_handles.window_modules.radial_average.theta_bin,'string'));
        if not(isempty(temp));
            if strcmp(status_flags.analysis_modules.radial_average.theta_bin_units,'pixels')
                status_flags.analysis_modules.radial_average.theta_bin_pixels = temp;
            elseif strcmp(status_flags.analysis_modules.radial_average.theta_bin_units,'absolute')
                status_flags.analysis_modules.radial_average.theta_bin_absolute = temp;
            elseif strcmp(status_flags.analysis_modules.radial_average.theta_bin_units,'resolution')
                status_flags.analysis_modules.radial_average.theta_bin_resolution = temp;
            end
        end
        
    case 'azimuth_bin'
        temp = str2num(get(grasp_handles.window_modules.radial_average.azimuth_bin,'string'));
        if not(isempty(temp));
            if strcmp(status_flags.analysis_modules.radial_average.azimuth_bin_units,'smart')
                
            elseif strcmp(status_flags.analysis_modules.radial_average.azimuth_bin_units,'absolute')
                status_flags.analysis_modules.radial_average.azimuth_bin_absolute = temp;
            end
        end
        
    case 'sector_mask_chk'
        status_flags.analysis_modules.radial_average.sector_mask_chk = not(status_flags.analysis_modules.radial_average.sector_mask_chk);
        
    case 'ellipse_mask_chk'
        status_flags.analysis_modules.radial_average.ellipse_mask_chk = not(status_flags.analysis_modules.radial_average.ellipse_mask_chk);
        
    case 'strip_mask_chk'
        status_flags.analysis_modules.radial_average.strip_mask_chk = not(status_flags.analysis_modules.radial_average.strip_mask_chk);
        
    case 'single_radio'
        status_flags.analysis_modules.radial_average.single_depth_radio = 0;
        set(grasp_handles.window_modules.radial_average.depth_combine_text,'visible','off');
        set(grasp_handles.window_modules.radial_average.depth_combine_check,'visible','off');
        
    case 'depth_radio'
        status_flags.analysis_modules.radial_average.single_depth_radio = 1;
        set(grasp_handles.window_modules.radial_average.depth_combine_text,'visible','on');
        set(grasp_handles.window_modules.radial_average.depth_combine_check,'visible','on');
        status_flags.selector.depth_range_chk = 1; %Switch on depth range in main display
        hide_stuff;
        
        
        
    case 'direct_to_file_check'
        status_flags.analysis_modules.radial_average.direct_to_file = get(gcbo,'value');
        
    case 'averaging_control'
        
        %temporary stuff
        if not(isfield(status_flags.analysis_modules.radial_average,'tof_iq'));
           status_flags.analysis_modules.radial_average.tof_iq =0;
        end
        if status_flags.analysis_modules.radial_average.tof_iq == 1;
            tof_iq
            return
        end
        
        

        
        
        
        if status_flags.analysis_modules.radial_average.single_depth_radio == 1; %i.e. do depth
            
            %Turn off command window parameter update for the boxing -
            %Doesnt make much different in speed ~ 3 % for text display & 10 % for graphic display
            %remember_display_params_status = status_flags.command_window.display_params;
            %remember_display_refresh_status = status_flags.display.refresh;
            %status_flags.command_window.display_params= 0;
            %status_flags.display.refresh = 0; %Turn off the 2D display for speed
            
            %Churn through the depth doing the averages
            index = data_index(status_flags.selector.fw);
            foreground_depth = status_flags.selector.fdpth_max - grasp_data(index).sum_allow;
            depth_start = status_flags.selector.fd; %remember the initial foreground depth
            
            disp(['Averaging Workhseets through Depth']);
            tof_iq_store = []; %Empty array
            
            if status_flags.selector.depth_range_min<= foreground_depth;
                d_start = status_flags.selector.depth_range_min;
            else
                d_start = foreground_depth; status_flags.selector.depth_range_min = foreground_depth;
            end
            if status_flags.selector.depth_range_max<= foreground_depth;
                d_end = status_flags.selector.depth_range_max;
            else
                d_end = foreground_depth; status_flags.selector.depth_range_max = foreground_depth;
            end
            
            if strcmp(status_flags.analysis_modules.radial_average.display_update,'on');
                status_flags.command_window.display_params=1;
                status_flags.display.refresh = 1;
            else
                status_flags.command_window.display_params=0;
                status_flags.display.refresh = 0;
            end
            
            for n = d_start:d_end
                message_handle = grasp_message(['Averaging Worksheets through Depth: ' num2str(n) ' of ' num2str(foreground_depth)]);
                status_flags.selector.fd = n+grasp_data(index).sum_allow;
                main_callbacks('depth_scroll'); %Scroll all linked depths and update
                radial_average_callbacks(options,'on',(n-d_start+1)); % On is the hold status
            end
            if ishandle(message_handle);
                delete(message_handle);
            end
            status_flags.display.refresh = 1;
            status_flags.command_window.display_params = 1;
            %Check for Depth (TOF) rebin
            if status_flags.analysis_modules.radial_average.single_depth_radio == 1 && status_flags.analysis_modules.radial_average.d33_tof_combine == 1;
                disp(' ');
                disp('Secondary Binning (e.g. TOF) of Data');
                radial_average_callbacks('d33_tof_rebin');
            end
            
            %Turn on command window parameter update for the boxing
            %status_flags.command_window.display_params = remember_display_params_status;
            %status_flags.display.refresh = remember_display_refresh_status;
            
            %Set depth selector back to as before
            status_flags.selector.fd = depth_start;
            grasp_update
            
        else %just a single average of the current display
            radial_average_callbacks(options,'off',1);
            %Check for Complext D33 TOF rebin
            if status_flags.analysis_modules.radial_average.single_depth_radio == 1 && status_flags.analysis_modules.radial_average.d33_tof_combine == 1;
                disp(' ');
                disp('Secondary Binning (e.g. TOF) of Data');
                radial_average_callbacks('d33_tof_rebin');
            end
          
        end
        
    case 'd33_tof_rebin'
        %This code is the same (almost) as found in grasp_plot_callbacks
        
        n_curves = length(tof_iq_store)
        
        %Compile all data points into a single list ready for re-binning
        iq_store = []; res_kernel_store_weight = []; res_kernel_store_x = [];
        for n = 1:n_curves
            %note:  exdat is just the gaussian equivalent or classic FWHM gaussian resolution
            iq_store = [iq_store; [tof_iq_store{n}.xdat, tof_iq_store{n}.ydat, tof_iq_store{n}.edat, tof_iq_store{n}.exdat]];
            res_kernel_store_weight = [res_kernel_store_weight, tof_iq_store{n}.resolution_kernels.weight];
            res_kernel_store_x = [res_kernel_store_x, tof_iq_store{n}.resolution_kernels.x];
        end
        
        %Calculate delta_q / q for the iq list
        iq_store(:,5) = iq_store(:,4)./iq_store(:,1);
        data_count_max = length(iq_store);
        
        q_min = min(iq_store(:,1)); q_max = max(iq_store(:,1));
        n_bins =  status_flags.analysis_modules.rebin.n_bins;
        
        %Bin Spacing
        if strcmp(status_flags.analysis_modules.rebin.bin_spacing,'log');
            %Log Bin Spacing
            logq_span =  ceil(log10(q_max)) - floor(log10(q_min));
            bin_step = logq_span/n_bins;
            log_edges = floor(log10(q_min)):bin_step:ceil(log10(q_max));
            bin_edges = 10.^log_edges;
        elseif strcmp(status_flags.analysis_modules.rebin.bin_spacing,'linear');
            %Linear Bin Spacing
            q_span = q_max - q_min;
            bin_step = q_span /n_bins;
            bin_edges = q_min:bin_step:q_max;
        end
        
        %delta_q/q regrouping bands
        dqq =  status_flags.analysis_modules.rebin.regroup_bands;
        
        %Display parameters in command window
        disp(' ')
        disp('Re-grouping Data')
        disp(['Bin spacing = '  status_flags.analysis_modules.rebin.bin_spacing]);
        disp(['# bins = ' num2str(status_flags.analysis_modules.rebin.n_bins)]);
        disp(['dq/q regouping bands = ' num2str(status_flags.analysis_modules.rebin.regroup_bands)]);
        disp(['Weighting average to di/i^' num2str(status_flags.analysis_modules.rebin.dii_power)]);
        disp(['Weighting average to dq/q^' num2str(status_flags.analysis_modules.rebin.dqq_power)]);
        disp(' ')
        
        data_count = 0;
        
        for n = 1:(length(dqq)-1);
            
            %find IQ points that fall within the resolution band
            temp = find(iq_store(:,5)>= dqq(n) & iq_store(:,5)<dqq(n+1));
            if not(isempty(temp));
                
                data_count = data_count + numel(temp);
                
                iq_store2 = iq_store(temp,:); %List of IQ points within the resolution band
                
                res_kernels.x = res_kernel_store_x(temp);
                res_kernels.weight = res_kernel_store_weight(temp);
                
                %Rebin
                temp = tof_rebin(iq_store2,res_kernels,bin_edges);
                
                %Plot  D33-TOF-Rebinned Data
                column_labels = ['Mod_Q   ' char(9) 'I       ' char(9) 'Err_I   ' char(9) 'FWHM_Q'];
                plot_info = struct(....
                    'plot_type','plot',....
                    'hold_graph','on',....
                    'plot_title',['Radial Re-grouping: |q|'],....
                    'x_label',['|q| (' char(197) '^{-1})'],....
                    'y_label',displayimage.units,....
                    'legend_str',['#' num2str(displayimage.params1(128))],....
                    'params',displayimage.params1,....
                    'parsub',displayimage.subtitle,....
                    'info',displayimage.info,....
                    'column_labels',column_labels);
                %plot_info.history = local_history;
                
                plot_info.plot_data.resolution_kernels.weight = temp.res_kernels.weight;
                plot_info.plot_data.resolution_kernels.x = temp.res_kernels.x;
                plot_info.plot_data.resolution_kernels.history = [];
                
                plot_info.plot_data.xdat = temp.array(:,1);
                plot_info.plot_data.ydat = temp.array(:,2);
                plot_info.plot_data.edat = temp.array(:,3);
                plot_info.plot_data.exdat = temp.res_kernels.fwhm; %Gaussian Equivalent of the weighted averaged kernels
                
                
                %Use old history to pass through
                plot_info.history = [];
                
                if strcmp(status_flags.subfigure.show_resolution,'on')
                    column_format = 'xyhe'; %Show resolution
                    %Note, need to swap colums around for the ploterr fn.  Need [x,y,err_x,err_y]
                    plotdata = [plot_info.plot_data.xdat,plot_info.plot_data.ydat,plot_info.plot_data.exdat,plot_info.plot_data.edat];
                else
                    column_format = 'xye'; %Do not show resolution.
                    %Note, uses Matlab errorbar fn.  Need [x,y,err_y]
                    plotdata = [plot_info.plot_data.xdat,plot_info.plot_data.ydat,plot_info.plot_data.edat];
                end
                
                export_data = [plot_info.plot_data.xdat,plot_info.plot_data.ydat,plot_info.plot_data.edat,plot_info.plot_data.exdat];

                plot_info.export_data = export_data; % replace the export data with the new math op data
                
                plot_info.legend_str = [num2str(dqq(n)) '-' num2str(dqq(n+1)) ' dq/q Resolution'];
                
                
                %Plot Radial Averaged Curves or Direct export to file
                if status_flags.analysis_modules.radial_average.direct_to_file == 0; %plot curves
                    grasp_plot(plotdata,column_format,plot_info);
                    
                else % Direct to file
                    disp('Exporting Radial Average Direct to File')
                    
                    %The code below is copied and modified from the export data
                    %routine in grasp_plot_menu_callbacks.
                    
                    %In the future a better single routine should be called for
                    %exporting data
                    
                    %Use different line terminators for PC or unix
                    if ispc; newline = 'pc'; terminator_str = [char(13) char(10)]; %CR/LF
                    else newline = 'unix'; terminator_str = [char(10)]; %LF
                    end
                    
                    %ONLY use Auto file numbering for 'direct to file'
                    %***** Build Output file name *****
                    numor_str = num2str(plot_info.params(128));
                    a = length(numor_str);
                    if a ==1; addzeros = '00000';
                    elseif a==2; addzeros = '0000';
                    elseif a==3; addzeros = '000';
                    elseif a==4; addzeros = '00';
                    elseif a==5; addzeros = '0';
                    elseif a==6; addzeros = ''; end
                    
                    fname = [addzeros numor_str '_tof.dat']; %options2 is the depth number
                    
                    %Open file for writing
                    disp(['Exporting data: '  grasp_env.path.project_dir fname]);
                    fid=fopen([grasp_env.path.project_dir fname],'wt');
                    
                    %Check if to include history header
                    if strcmp(status_flags.subfigure.export.data_history,'on');
                        history = plot_info.history;
                        
                        for m = 1:length(history)
                            textstring = history{m};
                            fprintf(fid,'%s \n',textstring);
                        end
                        fprintf(fid,'%s \n','');
                        fprintf(fid,'%s \n','');
                    end
                    
                    %Check if to include column labels
                    if strcmp(status_flags.subfigure.export.column_labels,'on')
                        if isfield(plot_info,'column_labels');
                            %Convert column labels to hwhm or fwhm if necessary
                            if strcmp(status_flags.subfigure.export.resolution_format,'hwhm') %Convert to hwhm
                                plot_info.column_labels = strrep(plot_info.column_labels,'FWHM_Q','HWHM_Q');
                            elseif strcmp(status_flags.subfigure.export.resolution_format,'sigma') %Convert to sigma
                                plot_info.column_labels = strrep(plot_info.column_labels,'FWHM_Q','Sigma_Q');
                            end
                            fprintf(fid,'%s \n',[plot_info.column_labels terminator_str]);
                            fprintf(fid,'%s \n','');
                        end
                    end
                    
                    %Strip out any Nans
                    temp = find(not(isnan(export_data(:,1))));
                    export_data = export_data(temp,:);
                    
                    %Check if to include q-reslution (4th column)
                    if strcmp(status_flags.subfigure.export.include_resolution,'on')
                        %Check what format of q-resolution, sigma, hwhm, fwhm
                        %Default coming though from Grasp is sigma
                        if strcmp(status_flags.subfigure.export.resolution_format,'hwhm') %Convert to hwhm
                            export_data(:,4) = export_data(:,4)/2; %hwhm
                        elseif strcmp(status_flags.subfigure.export.resolution_format,'sigma') %Convert to sigma
                            export_data(:,4) = export_data(:,4)/ (2 * sqrt(2 * log(2)));%fwhm
                            
                        end
                    else
                        disp('help here:  radial_average_callbacks 409 & grasp_plot_menu_callbacks line 189')
                    end
                    dlmwrite([grasp_env.path.project_dir fname],export_data,'delimiter','\t','newline',newline,'-append','precision',6);
                    fclose(fid);
                end
            end
        end
        redundancy = 1-data_count/data_count_max;
        disp(['Re-bin data redundancy: ' num2str(redundancy*100) '[%]'])
        disp([num2str(data_count) ' counts used'])
        disp([num2str(data_count_max-data_count) ' counts thrown in the dustbin'])


        
        
    case 'radial_q'
        
        local_history = displayimage.history;  %Take a copy of the history to modify though this process
        local_history = [local_history, {['***** Analysis *****']}];
        
        %Check for Sector Masks
        smask = [];
        if status_flags.analysis_modules.radial_average.sector_mask_chk ==1;
            %Check sector window is still open
            if ishandle(grasp_handles.window_modules.sector.window);
                smask = sector_callbacks('build_sector_mask');
            else
                status_flags.analysis_modules.radial_average.sector_mask_chk =0;
            end
        end
        
        %Check for Strip Masks
        strip_mask = [];
        if status_flags.analysis_modules.radial_average.strip_mask_chk ==1;
            %Check if strip window is still open
            if ishandle(grasp_handles.window_modules.strips.window);
                strip_mask = strips_callbacks('build_strip_mask');
            else
                status_flags.analysis_modules.radial_average.strip_mask_chk =0;
            end
        end
        
        %Radial and Azimuthal average takes the data directly from displayimage
        iq_data = []; %Final iq for all detectors appended together
        iq_big_list = [];
        
        for det = 1:inst_params.detectors
            
            if status_flags.display.(['axis' num2str(det) '_onoff']) ==1; %i.e. Detector is Active
                
                %Check current displayimage is not empty
                if sum(sum(displayimage.(['data' num2str(det)])))==0 && sum(displayimage.(['params' num2str(det)]))==0
                    disp(['Detector ' num2str(det) ' data and parameters are empty']);
                else
                    
                    %Prepare any aditional masks, e.g. sector & strip masks.
                    mask = displayimage.(['mask' num2str(det)]);  %This is the combined user & instrument mask
                    %Add the Sector Mask
                    if not(isempty(smask));
                        mask = mask.*smask.(['det' num2str(det)]);
                    end
                    %Add the Strip Mask
                    if not(isempty(strip_mask));
                        mask = mask.*strip_mask.(['det' num2str(det)]);
                    end
                    
                    %***** Turn 2D detector data into list(s) for re-binning *****
                    %Turn 2D data into a list for re-binning
                    temp = displayimage.(['qmatrix' num2str(det)])(:,:,5); %mod_q
                    temp2 = displayimage.(['qmatrix' num2str(det)])(:,:,13); %delta_q (FWHM) - Classic Resolution
                    temp3 = displayimage.(['qmatrix' num2str(det)])(:,:,11); %delta_q_lambda (FWHM)
                    temp4 = displayimage.(['qmatrix' num2str(det)])(:,:,12); %delta_q_theta (FWHM)
                    temp5 = displayimage.(['qmatrix' num2str(det)])(:,:,18); %delta_q_pixel (FWHM)
                    temp6 = displayimage.(['qmatrix' num2str(det)])(:,:,19); %delta_q_sample_aperture (FWHM)
                    
                    iq_list = [];
                    iq_list(:,1) = temp(logical(mask)); %mod q
                    iq_list(:,2) = displayimage.(['data' num2str(det)])(logical(mask)); %Intensity
                    iq_list(:,3) = displayimage.(['error' num2str(det)])(logical(mask)); %err_Intensity
                    iq_list(:,4) = temp2(logical(mask)); %delta_q FWHM
                    iq_list(:,5) = temp3(logical(mask)); %delta_q_lambda FWHM
                    iq_list(:,6) = temp4(logical(mask)); %delta_q_theta FWHM
                    iq_list(:,7) = temp5(logical(mask)); %delta_q_pixel FWHM
                    iq_list(:,8) = temp6(logical(mask)); %delta_q_sample_aperture FWHM
                    
                    iq_big_list = [iq_big_list; iq_list];
                    
                    
                    
                    
                end
            end
        end
        
        iq_list = iq_big_list;
        
        
        %Generate Bin_Edges
        x_min = min(iq_list(:,1)); x_max = max(iq_list(:,1));
        
        det =1; %So below takes parameters from main detector
        
        if strcmp(status_flags.analysis_modules.radial_average.q_bin_units,'pixels')
            bin_step = status_flags.analysis_modules.radial_average.q_bin_pixels;
            local_history = [local_history, {['Averaging I vs. Q.  Bin size:  ' num2str(status_flags.analysis_modules.radial_average.q_bin_pixels) ' [Pixel(s)]']}];
            
            %Calculate bin edges based on pixel steps across the detector *****
            %Calculate delta_q across 1 pixel at q=0
            delta_q = (4*pi / displayimage.(['params' num2str(det)])(inst_params.vectors.wav)) * ((inst_params.(['detector' num2str(det)]).pixel_size(1) *1e-3 * bin_step)/displayimage.(['params' num2str(det)])(inst_params.vectors.det))/2;
            bin_edges = x_min; bin_edge = x_min;
            while bin_edge < x_max
                bin_edge = bin_edge + delta_q;
                bin_edges = [bin_edges, bin_edge];
            end
            
            
        elseif strcmp(status_flags.analysis_modules.radial_average.q_bin_units,'absolute')
            bin_step = status_flags.analysis_modules.radial_average.q_bin_absolute;
            local_history = [local_history, {['Averaging I vs. Q.  Bin size:  ' num2str(status_flags.analysis_modules.radial_average.q_bin_absolute) ' [�-1]  '  status_flags.analysis_modules.radial_average.q_bin_absolute_scale]}];
            
            %Check if using linear or log bins
            if strcmp(status_flags.analysis_modules.radial_average.q_bin_absolute_scale,'linear')
                %Constant bin size across data_range
                bin_edges = x_min:bin_step:x_max;
            elseif strcmp(status_flags.analysis_modules.radial_average.q_bin_absolute_scale,'log10')
                log_edges = floor(log10(x_min)):bin_step:ceil(log10(x_max));
                bin_edges = 10.^log_edges;
            end
            
        elseif strcmp(status_flags.analysis_modules.radial_average.q_bin_units,'resolution');
            qmin = min(iq_list(:,1)); qmax = max(iq_list(:,1));
            if qmin==0; qmin = eps; end %this avoids an error when the beam centre has not been set and is left on 64,64
            bin_edges = [qmin];
            while bin_edges(length(bin_edges)) < qmax;
                %Find the closest data q-point to this q
                temp = iq_list(:,1) - bin_edges(length(bin_edges));
                temp = abs(temp);
                [temp, itemp] = min(temp);
                delta_q_fraction = iq_list(itemp,4) / status_flags.analysis_modules.radial_average.q_bin_resolution;
                bin_edges = [bin_edges, bin_edges(length(bin_edges))+delta_q_fraction];
                local_history = [local_history, {['Averaging I vs. Q.  Bin size:  ' num2str(status_flags.analysis_modules.radial_average.q_bin_resolution) ' [Fractional Resolution]  '  status_flags.analysis_modules.radial_average.q_bin_absolute_scale]}];
            end
        end
        
        
        
        if length(bin_edges) <2;
            disp('Error generating Bin_Edges - not enough Bins')
            disp('Please check re-binning paramters');
        end
        %***** Now re-bin *****
        if length(bin_edges) >=2;
            %temp = rebin([iq_list(:,1),iq_list(:,2),iq_list(:,3),iq_list(:,4)],bin_edges); %[q,I,errI,delta_q,pixel_count]
            temp = rebin(iq_list,bin_edges);
            iq_data = [iq_data; temp]; %append the iq data from the different detectors together
            %Note: Note, output from rebin is:
            %iq_data(:,1) = mod_q
            %iq_data(:,2) = Intensity
            %iq_data(:,3) = Err_Intensity
            %iq_data(:,4) = delta_q Classic q-resolution
            %iq_data(:,5) = delta_q Lambda
            %iq_data(:,6) = delta_q Theta
            %iq_data(:,7) = delta_q Detector Pixels
            %iq_data(:,8) = delta_q Sample Aperture
            
            %iq_data(:,9) = delta_q Binning (FWHM Square) - always next to last
            %iq_data(:,10) = # elements - always last
            
            %If enabled (ticked) add the binning resolution (FWHM) to the classic q-resolution (FWHM)
            if status_flags.resolution_control.binning_check == 1;
                %convert both back to sigma before adding in quadrature
                sigma1 = iq_data(:,4)/2.3548; %Came as a Gaussian FWHM
                sigma2 = iq_data(:,9)/3.4; %Came as a Square FWHM
                sigma = sqrt( sigma1.^2 + sigma2.^2 ); %Gaussian Equivalent
                fwhm = sigma * 2.3548;
                iq_data(:,4) = fwhm;
            end
            
        end
        
        %Check all the detector data wasn't empty
        if isempty(iq_data);
            disp(['All detector data was empty:  Nothing to rebin']);
            return
        end
        
        %***** Build Resolution Kernels for every q point *****
        kernel_data.fwhmwidth = status_flags.resolution_control.fwhmwidth;
        kernel_data.finesse = status_flags.resolution_control.finesse * status_flags.resolution_control.fwhmwidth;
        if not(isodd(kernel_data.finesse)); kernel_data.finesse = kernel_data.finesse+1; end %Finesse should be ODD number
        kernel_data.classic_res.fwhm = iq_data(:,4); kernel_data.classic_res.shape = 'Gaussian';
        
        kernel_data.lambda.fwhm = iq_data(:,5); kernel_data.lambda.shape = '';
        kernel_data.theta.fwhm = iq_data(:,6); kernel_data.theta.shape = 'tophat';
        kernel_data.pixel.fwhm = iq_data(:,7); kernel_data.pixel.shape = 'tophat';
        kernel_data.binning.fwhm = iq_data(:,9); kernel_data.binning.shape = 'tophat';
        kernel_data.aperture.fwhm = iq_data(:,8); kernel_data.aperture.shape = 'tophat';
        
        kernel_data.cm = displayimage.cm.(['det' num2str(det)]); %Send real beam profile in to use for resolution smearing
        %Build the kernels
        
        resolution_kernels = build_resolution_kernels(iq_data(:,1), kernel_data);
        
        
        %Prepare the Data to Plot, Export, or keep for D33_TOF_Rebin
        plot_data.xdat = iq_data(:,1);
        plot_data.ydat = iq_data(:,2);
        plot_data.edat = iq_data(:,3);
        plot_data.exdat = resolution_kernels.fwhm;
        plot_data.resolution_kernels = resolution_kernels;
        plot_data.no_elements = iq_data(:,9);
        
        %Check for Complext D33 TOF rebin
        if status_flags.analysis_modules.radial_average.single_depth_radio == 1 && status_flags.analysis_modules.radial_average.d33_tof_combine == 1;
            disp(['Depth Rebin - Building Pre-averages : Frame #' num2str(options2)])
            tof_iq_store{options2} = plot_data;
        else
            
            %***** Plot I vs Q ****
            column_labels = ['Mod_Q   ' char(9) 'I       ' char(9) 'Err_I   ' char(9) 'FWHM_Q'];
            
            export_data = iq_data(:,1:4); %[q, I, err_I, dq resolutuion(fwhm)]
            %classic resolution by default.
            %or
            %replace by gaussian equivalent resolution
            if status_flags.resolution_control.convolution_type == 1 || status_flags.resolution_control.convolution_type == 2; %Real shape kernel & gaussian equivalent
                export_data(:,4) = resolution_kernels.fwhm(:,1);
            end
            
            if strcmp(status_flags.subfigure.show_resolution,'on')
                column_format = 'xyhe'; %Show resolution
                %Note, need to swap colums around for the ploterr fn.  Need [x,y,err_x,err_y]
                plotdata = [iq_data(:,1:2),iq_data(:,4),iq_data(:,3)];% IvsQ, Horz Delta_q FWHM and Vert Err_I error bars
            else
                column_format = 'xye'; %Do not show resolution.
                %Note, uses Matlab errorbar fn.  Need [x,y,err_y]
                plotdata = [iq_data(:,1:3)];
            end
            
            plot_info = struct(....
                'plot_type','plot',....
                'hold_graph',options,....
                'plot_title',['Radial Re-grouping: |q|'],....
                'x_label',['|q| (' char(197) '^{-1})'],....
                'y_label',displayimage.units,....
                'legend_str',['#' num2str(displayimage.params1(128))],....
                'params',displayimage.params1,....
                'parsub',displayimage.subtitle,....
                'export_data',export_data,....
                'export_column_format','xyhe',....
                'plot_data',plot_data,....
                'info',displayimage.info,....
                'column_labels',column_labels);
            plot_info.history = local_history;
            %Plot Radial Averaged Curves or Direct export to file
            if status_flags.analysis_modules.radial_average.direct_to_file == 0; %plot curves
                
                grasp_plot(plotdata,column_format,plot_info);
                
            else % Direct to file
                disp('Exporting Radial Average Direct to File')
                
                %The code below is copied and modified from the export data
                %routine in grasp_plot_menu_callbacks.
                
                %In the future a better single routine should be called for
                %exporting data
                
                %Use different line terminators for PC or unix
                if ispc; newline = 'pc'; terminator_str = [char(13) char(10)]; %CR/LF
                else newline = 'unix'; terminator_str = [char(10)]; %LF
                end
                
                %ONLY use Auto file numbering for 'direct to file'
                %***** Build Output file name *****
                numor_str = num2str(plot_info.params(128));
                a = length(numor_str);
                if a ==1; addzeros = '00000';
                elseif a==2; addzeros = '0000';
                elseif a==3; addzeros = '000';
                elseif a==4; addzeros = '00';
                elseif a==5; addzeros = '0';
                elseif a==6; addzeros = ''; end
                
                fname = [addzeros numor_str '_' num2str(options2) '.dat']; %options2 is the depth number
                
                %Open file for writing
                disp(['Exporting data: '  grasp_env.path.project_dir fname]);
                fid=fopen([grasp_env.path.project_dir fname],'wt');
                
                %Check if to include history header
                if strcmp(status_flags.subfigure.export.data_history,'on');
                    history = plot_info.history;
                    
                    for m = 1:length(history)
                        textstring = history{m};
                        fprintf(fid,'%s \n',textstring);
                    end
                    fprintf(fid,'%s \n','');
                    fprintf(fid,'%s \n','');
                end
                
                export_data = plot_info.export_data;
                %Check if to include column labels
                if strcmp(status_flags.subfigure.export.column_labels,'on')
                    if isfield(plot_info,'column_labels');
                        %Convert column labels to hwhm or fwhm if necessary
                        if strcmp(status_flags.subfigure.export.resolution_format,'hwhm') %Convert to hwhm
                            plot_info.column_labels = strrep(plot_info.column_labels,'FWHM_Q','HWHM_Q');
                        elseif strcmp(status_flags.subfigure.export.resolution_format,'sigma') %Convert to sigma
                            plot_info.column_labels = strrep(plot_info.column_labels,'FWHM_Q','Sigma_Q');
                        end
                        fprintf(fid,'%s \n',[plot_info.column_labels terminator_str]);
                        fprintf(fid,'%s \n','');
                    end
                end
                
                %Strip out any Nans
                temp = find(not(isnan(export_data(:,1))));
                export_data = export_data(temp,:);
                
                %Check if to include q-reslution (4th column)
                if strcmp(status_flags.subfigure.export.include_resolution,'on')
                    %Check what format of q-resolution, sigma, hwhm, fwhm
                    %Default coming though from Grasp is sigma
                    if strcmp(status_flags.subfigure.export.resolution_format,'hwhm') %Convert to hwhm
                        export_data(:,4) = export_data(:,4)/2; %hwhm
                    elseif strcmp(status_flags.subfigure.export.resolution_format,'sigma') %Convert to sigma
                        export_data(:,4) = export_data(:,4)/ (2 * sqrt(2 * log(2)));%fwhm
                        
                    end
                else
                    disp('help here:  radial_average_callbacks 409 & grasp_plot_menu_callbacks line 189')
                end
                dlmwrite([grasp_env.path.project_dir fname],export_data,'delimiter','\t','newline',newline,'-append','precision',6);
                fclose(fid);
            end
        end
        
        
    case 'radial_theta'
        
        local_history = displayimage.history;  %Take a copy of the history to modify though this process
        local_history = [local_history, {['***** Analysis *****']}];
        
        %Check for Sector Masks
        smask = [];
        if status_flags.analysis_modules.radial_average.sector_mask_chk ==1;
            %Check sector window is still open
            if ishandle(grasp_handles.window_modules.sector.window);
                smask = sector_callbacks('build_sector_mask');
            else
                status_flags.analysis_modules.radial_average.sector_mask_chk =0;
            end
        end
        
        %Check for Strip Masks
        strip_mask = [];
        if status_flags.analysis_modules.radial_average.strip_mask_chk ==1;
            %Check if strip window is still open
            if ishandle(grasp_handles.window_modules.strips.window);
                strip_mask = strips_callbacks('build_strip_mask');
            else
                status_flags.analysis_modules.radial_average.strip_mask_chk =0;
            end
        end
        
        %Radial and Azimuthal average takes the data directly from displayimage
        iq_data = []; %Final iq for all detectors appended together
        for det = 1:inst_params.detectors
            
                        if status_flags.display.(['axis' num2str(det) '_onoff']) ==1; %i.e. Detector is Active

            
            %Check current displayimage is not empty
            if sum(sum(displayimage.(['data' num2str(det)])))==0 && sum(displayimage.(['params' num2str(det)]))==0
                disp(['Detector ' num2str(det) ' data and parameters are empty']);
            else
                
                %Prepare any aditional masks, e.g. sector & strip masks.
                mask = displayimage.(['mask' num2str(det)]);  %This is the combined user & instrument mask
                %Add the Sector Mask
                if not(isempty(smask));
                    mask = mask.*smask.(['det' num2str(det)]);
                end
                %Add the Strip Mask
                if not(isempty(strip_mask));
                    mask = mask.*strip_mask.(['det' num2str(det)]);
                end
                
                %***** Turn 2D detector data into list(s) for re-binning *****
                %Turn 2D data into a list for re-binning
                temp = displayimage.(['qmatrix' num2str(det)])(:,:,9); %mod_2theta
                temp2 = displayimage.(['qmatrix' num2str(det)])(:,:,17)*2; %delta_2theta (FWHM)
                iq_list = [];
                iq_list(:,1) = temp(logical(mask)); %mod_2theta
                iq_list(:,2) = displayimage.(['data' num2str(det)])(logical(mask)); %Intensity
                iq_list(:,3) = displayimage.(['error' num2str(det)])(logical(mask)); %err_Intensity
                iq_list(:,4) = temp2(logical(mask)); %delta_theta
                
                %Generate Bin_Edges
                x_min = min(iq_list(:,1)); x_max = max(iq_list(:,1));
                
                if strcmp(status_flags.analysis_modules.radial_average.theta_bin_units,'pixels')
                    bin_step = status_flags.analysis_modules.radial_average.theta_bin_pixels;
                    local_history = [local_history, {['Averaging I vs. 2Theta.  Bin size:  ' num2str(status_flags.analysis_modules.radial_average.theta_bin_pixels) ' [Pixel(s)]']}];
                    
                    %Calculate bin edges based on pixel steps across the detector *****
                    %Calculate delta_2theta across 1 pixel at q=0
                    delta_2theta = (180/pi)*((inst_params.(['detector' num2str(det)]).pixel_size(1) *1e-3 * bin_step)/displayimage.(['params' num2str(det)])(inst_params.vectors.det));
                    bin_edges = x_min; bin_edge = x_min;
                    while bin_edge < x_max
                        bin_edge = bin_edge + delta_2theta;
                        bin_edges = [bin_edges, bin_edge];
                    end
                    
                elseif strcmp(status_flags.analysis_modules.radial_average.theta_bin_units,'absolute')
                    bin_step = status_flags.analysis_modules.radial_average.theta_bin_absolute;
                    local_history = [local_history, {['Averaging I vs. 2Theta.  Bin size:  ' num2str(status_flags.analysis_modules.radial_average.theta_bin_absolute) ' [�-1]  '  status_flags.analysis_modules.radial_average.theta_bin_absolute_scale]}];
                    %Check if using linear or log bins
                    if strcmp(status_flags.analysis_modules.radial_average.theta_bin_absolute_scale,'linear')
                        %Constant bin size across data_range
                        bin_edges = x_min:bin_step:x_max;
                    elseif strcmp(status_flags.analysis_modules.radial_average.theta_bin_absolute_scale,'log10')
                        log_edges = floor(log10(x_min)):bin_step:ceil(log10(x_max));
                        bin_edges = 10.^log_edges;
                    end
                    
                elseif strcmp(status_flags.analysis_modules.radial_average.theta_bin_units,'resolution');
                    qmin = min(iq_list(:,1)); qmax = max(iq_list(:,1));
                    if qmin==0; qmin = eps; end %this avoids an error when the beam centre has not been set and is left on 64,64
                    bin_edges = [qmin];
                    while bin_edges(length(bin_edges)) < qmax;
                        %Find the closest data q-point to this q
                        temp = iq_list(:,1) - bin_edges(length(bin_edges));
                        temp = abs(temp);
                        [temp, itemp] = min(temp);
                        delta_q_fraction = iq_list(itemp,4) / status_flags.analysis_modules.radial_average.q_bin_resolution;
                        bin_edges = [bin_edges, bin_edges(length(bin_edges))+delta_q_fraction];
                        local_history = [local_history, {['Averaging I vs. Q.  Bin size:  ' num2str(status_flags.analysis_modules.radial_average.q_bin_resolution) ' [Fractional Resolution]  '  status_flags.analysis_modules.radial_average.q_bin_absolute_scale]}];
                    end
                end
                
                if length(bin_edges) <2;
                    disp('Error generating Bin_Edges - not enough Bins')
                    disp('Please check re-binning paramters');
                end
                
                
                %***** Now re-bin *****
                if length(bin_edges) >=2;
                    temp = rebin([iq_list(:,1),iq_list(:,2),iq_list(:,3),iq_list(:,4)],bin_edges); %[two_theta,I,errI,delta_two_theta,pixel_count]
                    iq_data = [iq_data; temp]; %append the iq data from the different detectors together
                end
            end
                        end
                        
        end
        
        
        %Check all the detector data wasn't empty
        if isempty(iq_data);
            disp(['All detector data was empty:  Nothing to rebin']);
            return
        end
        
        %***** Plot I vs Q ****
        column_labels = ['Two_Theta   ' char(9) 'I       ' char(9) 'Err_I   ' char(9) 'FWHM 2Theta'];
        export_data = iq_data(:,1:4); %[2theta, I, err_I, d2theta, pixelcount]
        if strcmp(status_flags.subfigure.show_resolution,'on')
            column_format = 'xyhe'; %Show resolution
            plotdata = [iq_data(:,1:2),iq_data(:,4),iq_data(:,3)];% IvsQ, Horz Delta_q and Vert Err_I error bars
        else
            column_format = 'xye'; %Do not show resolution
            plotdata = [iq_data(:,1:3)];
        end
        
        plot_info = struct(....
            'plot_type','plot',....
            'hold_graph',options,....
            'plot_title',['Radial Re-grouping: 2Theta'],....
            'x_label',['2Theta [degrees]'],....
            'y_label',displayimage.units,....
            'legend_str',['#' num2str(displayimage.params1(128))],....
            'params',displayimage.params1,....
            'parsub',displayimage.subtitle,....
            'export_data',export_data,....
            'info',displayimage.info,....
            'column_labels',column_labels);
        plot_info.history = local_history;
        grasp_plot(plotdata,column_format,plot_info);
        
        
        
        
        case 'azimuthal'
            
            local_history = displayimage.history;  %Take a copy of the history to modify though this process
            local_history = [local_history, {['***** Analysis *****']}];
            
            %Check for Sector Masks
            smask = [];
            if status_flags.analysis_modules.radial_average.sector_mask_chk ==1;
                %Check sector window is still open
               % if ishandle(grasp_handles.window_modules.sector.window);
                    smask = sector_callbacks('build_sector_mask');
               % else
                    status_flags.analysis_modules.radial_average.sector_mask_chk =0;
                %end
            end
            
            %Check for Strip Masks
            strip_mask = [];
            if status_flags.analysis_modules.radial_average.strip_mask_chk ==1;
                %Check if strip window is still open
                if ishandle(grasp_handles.window_modules.strips.window);
                    strip_mask = strips_callbacks('build_strip_mask');
                else
                    status_flags.analysis_modules.radial_average.strip_mask_chk =0;
                end
            end
            
            %         %Check if Nan Mask exists (generate by the 'divide' data sets algorythm)
            %         if isfield(displayimage,'nan_mask');
            %             mask = mask .* displayimage.mask;
            %         end
            
            
            %Radial and Azimuthal average takes the data directly from displayimage
            iq_data = []; %Final iq for all detectors appended together
            for det = 1:inst_params.detectors
                
                %Check current displayimage is not empty
                if sum(sum(displayimage.(['data' num2str(det)])))==0 && sum(displayimage.(['params' num2str(det)]))==0
                    disp(['Detector ' num2str(det) ' data and parameters are empty']);
                else
                    
                    %Prepare any aditional masks, e.g. sector & strip masks.
                    mask = displayimage.(['mask' num2str(det)]);  %This is the combined user & instrument mask
                    %Add the Sector Mask
                    if not(isempty(smask));
                        mask = mask.*smask.(['det' num2str(det)]);
                    end
                    %Add the Strip Mask
                    if not(isempty(strip_mask));
                        mask = mask.*strip_mask.(['det' num2str(det)]);
                    end
                    
                    %***** Turn 2D detector data into list(s) for re-binning *****
                    %Turn 2D data into a list for re-binning
                    temp = displayimage.(['qmatrix' num2str(det)])(:,:,6); %q_angle
                    %temp2 = displayimage.(['qmatrix' num2str(det)])(:,:,13); %delta_q (sigma)
                    iq_list = [];
                    iq_list(:,1) = temp(logical(mask)); %mod q
                    iq_list(:,2) = displayimage.(['data' num2str(det)])(logical(mask)); %Intensity
                    iq_list(:,3) = displayimage.(['error' num2str(det)])(logical(mask)); %err_Intensity
                    %iq_list(:,4) = temp2(logical(mask)); %delta_q_angle
                    
                    
                    %Generate azimuthal bins
                    azimuthal_min = min(iq_list(:,3)); azimuthal_max = max(iq_list(:,3));
                    if strcmp(status_flags.analysis_modules.radial_average.azimuth_bin_units,'absolute')
                        bin_step = status_flags.analysis_modules.radial_average.azimuth_bin_absolute;
                        bin_edges = 0:bin_step:360;
                    end
                    local_history = [local_history, {['Averaging I vs. Azimuth.  Bin size:  ' num2str(status_flags.analysis_modules.radial_average.azimuth_bin_absolute) ' [Degree(s)]']}];
                    
                    if length(bin_edges) <2;
                        disp('Error generating Bin_Edges - not enough Bins')
                        disp('Please check re-binning paramters');
                    end
                    
                    %***** Now re-bin *****
                    if length(bin_edges) >=2;
                        temp = rebin([iq_list(:,1),iq_list(:,2),iq_list(:,3)],bin_edges); %[q,I,errI,delta_q,pixel_count]
                        iq_data = [iq_data; temp]; %append the iq data from the different detectors together
                    end
                end
            end
            
            
            
            %Check all the detector data wasn't empty
            if isempty(iq_data);
                disp(['All detector data was empty:  Nothing to rebin']);
                return
            end
            
            %***** Plot I vs Angle ****
            column_labels = ['Az_Angle   ' char(9) 'I       ' char(9) 'Err_I   '];
            export_data = iq_data(:,1:4); %[Az_Angle, I, err_I, pixelcount]
            if strcmp(status_flags.subfigure.show_resolution,'on')
                column_format = 'xye'; %Show resolution
                plotdata = [iq_data(:,1:2),iq_data(:,3)];% IvsQ, Horz Delta_q and Vert Err_I error bars
            else
                column_format = 'xye'; %Do not show resolution
                plotdata = [iq_data(:,1:3)];
            end
            
            plot_info = struct(....
                'plot_type','plot',....
                'hold_graph',options,....
                'plot_title',['Azimuthal Re-grouping: |degrees|'],....
                'x_label',['Azimuthal Angle [degrees]'],....
                'y_label',displayimage.units,....
                'legend_str',['#' num2str(displayimage.params1(128))],....
                'params',displayimage.params1,....
                'parsub',displayimage.subtitle,....
                'export_data',export_data,....
                'info',displayimage.info,....
                'column_labels',column_labels);
            plot_info.history = local_history;
            grasp_plot(plotdata,column_format,plot_info);
end




%General update stuff if the Radial average window exists
if isfield(grasp_handles.window_modules.radial_average,'window')
    if ishandle(grasp_handles.window_modules.radial_average.window)
        %Update displayed Radial Average objects
        if strcmp(status_flags.analysis_modules.radial_average.q_bin_units,'pixels')
            bin_string = 'Radial Bin (pixels):';
            bin_size = status_flags.analysis_modules.radial_average.q_bin_pixels;
        elseif strcmp(status_flags.analysis_modules.radial_average.q_bin_units,'absolute')
            bin_string = 'Radial Bin (q):';
            bin_size = status_flags.analysis_modules.radial_average.q_bin_absolute;
        elseif strcmp(status_flags.analysis_modules.radial_average.q_bin_units,'resolution')
            bin_string = 'Radial Bin (D_q):';
            bin_size = status_flags.analysis_modules.radial_average.q_bin_resolution;
        end
        if ishandle(grasp_handles.window_modules.radial_average.q_bin_text)
            set(grasp_handles.window_modules.radial_average.q_bin_text,'string',bin_string);
        end
        if ishandle(grasp_handles.window_modules.radial_average.q_bin)
            set(grasp_handles.window_modules.radial_average.q_bin,'string',num2str(bin_size));
        end
        
        
        if strcmp(status_flags.analysis_modules.radial_average.theta_bin_units,'pixels')
            bin_string = 'Radial Bin (pixels):';
            bin_size = status_flags.analysis_modules.radial_average.theta_bin_pixels;
        elseif strcmp(status_flags.analysis_modules.radial_average.theta_bin_units,'absolute')
            bin_string = 'Radial Bin (2theta):';
            bin_size = status_flags.analysis_modules.radial_average.theta_bin_absolute;
        elseif strcmp(status_flags.analysis_modules.radial_average.theta_bin_units,'resolution')
            bin_string = 'Radial Bin (D_2T):';
            bin_size = status_flags.analysis_modules.radial_average.theta_bin_resolution;
        end
        if ishandle(grasp_handles.window_modules.radial_average.theta_bin_text)
            set(grasp_handles.window_modules.radial_average.theta_bin_text,'string',bin_string);
        end
        if ishandle(grasp_handles.window_modules.radial_average.theta_bin)
            set(grasp_handles.window_modules.radial_average.theta_bin,'string',num2str(bin_size));
        end
        
        
        if strcmp(status_flags.analysis_modules.radial_average.azimuth_bin_units,'smart')
            bin_string = 'Angle Bin (smart):';
            bin_size = 0;
        elseif strcmp(status_flags.analysis_modules.radial_average.azimuth_bin_units,'absolute')
            bin_string = 'Angle Bin (degs):';
            bin_size = status_flags.analysis_modules.radial_average.azimuth_bin_absolute;
        end
        if ishandle(grasp_handles.window_modules.radial_average.azimuth_bin_text)
            set(grasp_handles.window_modules.radial_average.azimuth_bin_text,'string',bin_string);
        end
        if ishandle(grasp_handles.window_modules.radial_average.azimuth_bin)
            set(grasp_handles.window_modules.radial_average.azimuth_bin,'string',num2str(bin_size));
        end
        
        if ishandle(grasp_handles.window_modules.radial_average.sector_mask_chk)
            set(grasp_handles.window_modules.radial_average.sector_mask_chk,'value',status_flags.analysis_modules.radial_average.sector_mask_chk);
        end
        if ishandle(grasp_handles.window_modules.radial_average.strip_mask_chk)
            set(grasp_handles.window_modules.radial_average.strip_mask_chk,'value',status_flags.analysis_modules.radial_average.strip_mask_chk);
        end
        %set(grasp_handles.window_modules.radial_average.ellipse_mask_chk,'value',status_flags.analysis_modules.radial_average.ellipse_mask_chk);
        
        %Single or Depth & Combine
        if status_flags.analysis_modules.radial_average.single_depth_radio == 1; %Depth
            single_radio_value = 0; depth_radio_value = 1;
            %combine_visible = 'on';
        else
            single_radio_value = 1; depth_radio_value = 0;
            %combine_visible = 'off';
        end
        if ishandle(grasp_handles.window_modules.radial_average.radio_single)
            set(grasp_handles.window_modules.radial_average.radio_single,'value',single_radio_value);
        end
        if ishandle(grasp_handles.window_modules.radial_average.radio_depth)
            set(grasp_handles.window_modules.radial_average.radio_depth,'value',depth_radio_value);
        end
        %         if ishandle(grasp_handles.window_modules.radial_average.combine_check)
        %             set(grasp_handles.window_modules.radial_average.combine_check,'visible',combine_visible,'value',status_flags.analysis_modules.radial_average.depth_combine);
        %         end
        
    end
end

