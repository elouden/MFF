function [phi, Int, Int_err] = getData(img_num)
        global status_flags
        global mf_fitter
        global plot_info
        global grasp_handles
    
        % get binning
        curAB = get(grasp_handles.window_modules.radial_average.azimuth_bin,'String');
       
        % load desired depth file and update main grasp GUI
        % the +1 is necessary as file 1 represents a sum
        status_flags.selector.fd = img_num+1;
        
        % set binning to 0.1 to get raw data, data will be rebinned later
        if(mf_fitter.smoothing.fwhm ~= 1)
            set(grasp_handles.window_modules.radial_average.azimuth_bin,'String',0.1);
        end
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

