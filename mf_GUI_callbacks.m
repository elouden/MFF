function [ ] = mf_GUI_callbacks( to_do )

%Callbacks for mf_GUI_window

% initialize
% fit_dropdown
% plot_dropdown
% good_peaks: 1-peak, 2-peak, 3-peak file
% change_int_cutoff
% change_dphi_cutoff
% go
% apply_sector
% change_ext
% change_cyc

% save & view params have their own functions

% v. 8.0
% 11/12/2016 MFF Liz

global mf_fitter
global status_flags
global displayimage
global grasp_handles
global grasp_env

switch to_do

    case 'initialize'
        
        notification = msgbox('Loading files...');
        % Before loading, first initialize in case previous names still loaded
        mf_fitter.fit_data.names = [];
        mf_fitter.fit_data.cycles = [];
        mf_fitter.int_cutoff = 0.12;
        mf_fitter.dphi_cutoff = 2;
        mf_fitter.depth = (status_flags.selector.fdpth_max - 1);
        
        %initialize numor array for future convenience
        i = 1;
        while(i<=mf_fitter.depth)
            
            status_flags.selector.fd = i+1;  %the +1 is necessary because position 1 corresponds to having all of the data files summed
            grasp_update;

            mf_fitter.fit_data.names(i,1) = displayimage.params1(128);
            
            i= i + 1;
        end
        
        %Set initial good peaks as first and last files
        mf_fitter.good_peaks.inner = mf_fitter.fit_data.names(1,1);
        mf_fitter.good_peaks.outers = mf_fitter.fit_data.names(mf_fitter.depth, 1);
        
        % Set default smoothing values
        mf_fitter.smoothing.fwhm = 1.25;
        mf_fitter.smoothing.step = 1.5;
        
        % Cycle Views
        mf_fitter.cycleview.view1 = 0;
        mf_fitter.cycleview.view2 = 0;
        mf_fitter.cycleview.view3 = 0;
        mf_fitter.cycleview.view4 = 0;
        mf_fitter.cycleview.view5 = 1;
        
        % Save files
        mf_fitter.save.fig = 0;
        mf_fitter.save.jpg = 0;
        mf_fitter.save.eps = 0;
        mf_fitter.save.pdf = 0;
        
        % old code to play "everything is awesome"
        %[mf_fitter.awesome.y, mf_fitter.awesome.Fs] = audioread('sound_fx_package.mp3');
        close(notification);
        
        mf_fitter.extension = [];
        mf_fitter.folder = [];
    
        
    case 'fit_dropdown'
        disp('setting fit paradigm')
        
        % Fitting Paradigms
        if ( get(mf_fitter.handles.fit_dropdown,'Value') == 1 )
            disp('Please select a fitting option')
            
        % Peaks Fixed
        elseif( get(mf_fitter.handles.fit_dropdown,'Value') == 2 )
            disp('Peak centers are fixed');
        
        % 2K Peaks Grow
        elseif( get(mf_fitter.handles.fit_dropdown,'Value') == 3 )
            disp('GS peaks nucleate in final positions and get larger - fits 3 peaks')
            
        % 14K Peaks Rotate
        elseif( get(mf_fitter.handles.fit_dropdown,'Value') == 4 )
            disp('MS peaks rotate in to final position - fits 2 peaks')
            
            % Fitting cycles to view window
            fig_position = [1000, 580, 350, 150];
            grasp_handles.window_modules.fit_cycles_4.window = figure(....
                'units','pixels',....
                'Position',fig_position,....
                'Name','View Cycles' ,....
                'NumberTitle', 'off',....
                'Tag','mf_cyc_window',...
                'color',grasp_env.background_color,....
                'menubar','none',....
                'resize','off');
            
            handle = grasp_handles.window_modules.fit_cycles_4.window;
            
            %Question           
            uicontrol(handle,'units','normalized','Position',[0.05 0.55 0.75 0.4],'FontName',grasp_env.font,'FontSize',2*grasp_env.fontsize,'HorizontalAlignment','left','Style','text','String','Which fitting cycles would you like to view?','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
            
            %Checkboxes
            val =  mf_fitter.cycleview.view1;
            uicontrol(handle,'units','normalized','Position',[0.2 0.3 0.3 0.25],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','left','Style','text','String','Ref Files?','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
            mf_fitter.handles.view.view1 = uicontrol(handle,'units','normalized','Position',[0.2 0.25 0.05 0.1],'tooltip','Check if you would like to view plots from this fitting cycle','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','checkbox','HorizontalAlignment','right','Tag','angular_bin','Visible','on','Value', val,'Callback','mf_GUI_callbacks(''view_cycles'')');

            val =  mf_fitter.cycleview.view2;
            uicontrol(handle,'units','normalized','Position',[0.35 0.3 0.3 0.25],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','left','Style','text','String','Cyc 1?','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
            mf_fitter.handles.view.view2 = uicontrol(handle,'units','normalized','Position',[0.35 0.25 0.05 0.1],'tooltip','Check if you would like to view plots from this fitting cycle','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','checkbox','HorizontalAlignment','right','Tag','angular_bin','Visible','on','Value', val,'Callback','mf_GUI_callbacks(''view_cycles'')');
        
            val = mf_fitter.cycleview.view3;
            uicontrol(handle,'units','normalized','Position',[0.5 0.3 0.3 0.25],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','left','Style','text','String','Cyc 2?','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
            mf_fitter.handles.view.view3 = uicontrol(handle,'units','normalized','Position',[0.5 0.25 0.05 0.1],'tooltip','Check if you would like to view plots from this fitting cycle','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','checkbox','HorizontalAlignment','right','Tag','angular_bin','Visible','on','Value', val,'Callback','mf_GUI_callbacks(''view_cycles'')');
            
            val =  mf_fitter.cycleview.view4;
            uicontrol(handle,'units','normalized','Position',[0.65 0.3 0.3 0.25],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','left','Style','text','String','Cyc 3?','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
            mf_fitter.handles.view.view4 = uicontrol(handle,'units','normalized','Position',[0.65 0.25 0.05 0.1],'tooltip','Check if you would like to view plots from this fitting cycle','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','checkbox','HorizontalAlignment','right','Tag','angular_bin','Visible','on','Value', val,'Callback','mf_GUI_callbacks(''view_cycles'')');
           
            val =  mf_fitter.cycleview.view5;
            uicontrol(handle,'units','normalized','Position',[0.8 0.3 0.3 0.25],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','left','Style','text','String','Cyc 4?','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
            mf_fitter.handles.view.view5 = uicontrol(handle,'units','normalized','Position',[0.8 0.25 0.05 0.1],'tooltip','Check if you would like to view plots from this fitting cycle','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','checkbox','HorizontalAlignment','right','Tag','angular_bin','Visible','on','Value', val,'Callback','mf_GUI_callbacks(''view_cycles'')');

            
        % Raster Scan
        elseif( get(mf_fitter.handles.fit_dropdown,'Value') == 5 )
            disp('Perform a 2D analysis')

        % Angle Binning Scan
        elseif( get(mf_fitter.handles.fit_dropdown,'Value') == 6 )
                choice = questdlg('Would you like to perform MFF over several angle binning values?  Note: you should enter the experiment applied cycles before running this.', ...
                    'Angle Binning Iteration?','Yes','No', 'No');
                switch choice 
                    case 'Yes'
                         angle_binning_iteration('open');
                    case 'No'
                        disp('ok');
                end
                
        % I vs Q
        elseif( get(mf_fitter.handles.fit_dropdown,'Value') == 7 )
                disp('Fitting I vs Q')
                
        % Error
        else
            disp('error - please select an option')           
        end
        
        
    case 'plot_dropdown'
        
        % Plot Options
        if ( get(mf_fitter.handles.plot_dropdown,'Value') == 1 )
            disp('Please select an option')
            
        % FMS
        elseif( get(mf_fitter.handles.plot_dropdown,'Value') == 2 )
            mf_fitter_callbacks('fraction');
            
        % Peak Centers
        elseif( get(mf_fitter.handles.plot_dropdown,'Value') == 3 )
            mf_fitter_callbacks('center_plot');
           
        % Peak Separation
        elseif( get(mf_fitter.handles.plot_dropdown,'Value') == 4 )
            mf_fitter_callbacks('center_separation');
            x = mf_fitter.fit_data.cycles;
            y  = abs(mf_fitter.fit_data.center3(:,1) - mf_fitter.fit_data.center1(:,1));
            y_err = sqrt((mf_fitter.fit_data.center3(:,2)).^2+(mf_fitter.fit_data.center1(:,2)).^2);
            plotData(x, y, y_err, 'Applied || AC Cycles', '\Delta \phi (degrees)', ['Peak Separation - ' mf_fitter.folder])
            set(gca,'xScale','log')
                  
        % Intensity Colormaps
        elseif( get(mf_fitter.handles.plot_dropdown,'Value') == 5 )
            %mf_fitter_callbacks('colormap'); 
            
if(mf_fitter.handles.smoothing.switch  & isempty(mf_fitter.SmoothedData) )
    disp(['Perform MRE smoothing with FWHM: ' num2str(mf_fitter.smoothing.fwhm)]);
    for i = 1:length(mf_fitter.fit_data.names)
        [phi, Int(i,:), Int_err(i,:)] = getData(i);
        
        
        phisym = status_flags.analysis_modules.sectors.theta;
        dTheta = status_flags.analysis_modules.sectors.delta_theta;
        ss = mf_fitter.smoothing.step;
        
        %phirng = [phisym - dTheta/2; phisym + dTheta/2; ss]   % because we are centering, do not need the phisym value
        phirng = [-dTheta/2; dTheta/2; ss]
        [phiS, IntS(i,:), Int_errS(i,:)] = smooth(phi, Int(i,:), Int_err(i,:), i, phisym, 0, phirng, mf_fitter.smoothing.fwhm)
        
        fig_h = plotData(phiS, IntS(i,:), Int_errS(i,:), '\phi  - \phi_0', 'Intensity (arb. units)', ['Smoothed Data Numor - ' num2str(i)]);
        close(fig_h)
    end
    
    mf_fitter.SmoothedData.phi = phiS;
    mf_fitter.SmoothedData.Int = IntS;
    mf_fitter.SmoothedData.Int_err = Int_errS;
end
            phi = mf_fitter.SmoothedData.phi;
            Int = mf_fitter.SmoothedData.Int;
            cm_x = mf_fitter.fit_data.cycles;
            CM(phi, Int, cm_x);
        
        % Peak Decay
        elseif( get(mf_fitter.handles.plot_dropdown,'Value') == 6 )
           % mf_fitter_callbacks('peak_decay');
           mf_fitter_peakDecay(-20,20) 
           
        else
            disp('error')           
        end
        
        
        
        
    case 'good_peaks'
        %Changes good peaks
        temp1 = str2double(get(mf_fitter.handles.good_inner, 'String'));
        if not(isempty(temp1))
            mf_fitter.good_peaks.inner = temp1;
        end
        temp2 = str2double(get(mf_fitter.handles.good_outers, 'String'));
        if not(isempty(temp2))
            mf_fitter.good_peaks.outers = temp2;
        end
        
        
    case 'change_int_cutoff'      
        %Changes intensity cutoff        
        temp = str2double(get(mf_fitter.handles.int_cutoff, 'String'));
        if not(isempty(temp)) 
            mf_fitter.int_cutoff = temp;
        end

        
     case 'change_dphi_cutoff'        
        %Changes multiplicative cutoff factor for peak separation        
        temp = str2num(get(mf_fitter.handles.dphi_cutoff, 'String'));
        if not(isempty(temp))            
            mf_fitter.dphi_cutoff = temp;           
        end
        
        
    case 'go'
        % Double check the user (cough cough Liz) has correctly changed all the settings
        choice = questdlg('Have you adjusted all the settings?  experiment name, beam center, fit paragdigm, cycles, extension name, sector, angle binning, colormap limits','Double Check','Yes','No', 'No');
        switch choice         
            case 'Yes'
            disp('setting fit paradigm')

            % Fitting Paradigms
            if ( get(mf_fitter.handles.fit_dropdown,'Value') == 1 )
                disp('Please select a fitting option')

            % Peaks Fixed
            elseif( get(mf_fitter.handles.fit_dropdown,'Value') == 2 )
                disp('Peak centers are fixed');
                %mf_fitter_callbacks('fraction');

            % 2K Peaks Grow
            elseif( get(mf_fitter.handles.fit_dropdown,'Value') == 3 )
                disp('GS peaks nucleate in final positions and get larger - fits 3 peaks')
                %mf_fitter_callbacks('center_plot');
                %mf_fitter_mpf2();
                mf_fitter_mpf2K();

            % 14K Peaks Rotate
            elseif( get(mf_fitter.handles.fit_dropdown,'Value') == 4 )
                disp('TEST MS peaks rotate in to final position - fits 2 peaks')
                mf_fitter_mpf3();

            % Raster Scan
            elseif( get(mf_fitter.handles.fit_dropdown,'Value') == 5 )
                disp('Perform a 2D analysis')
                %mf_fitter_callbacks('peak_decay');

            % Angle Binning Scan
            elseif( get(mf_fitter.handles.fit_dropdown,'Value') == 6 )
                    choice = questdlg('Would you like to perform MFF over several angle binning values?  Note: you should enter the experiment applied cycles before running this.', ...
                        'Angle Binning Iteration?','Yes','No', 'No');
                    switch choice 
                        case 'Yes'
                             angle_binning_iteration('open');
                        case 'No'
                            disp('ok');
                    end
            
            elseif( get(mf_fitter.handles.fit_dropdown,'Value') == 7 )
                mf_fitter_mpf_q();
                
        else
            disp('error')           
        end
        
            case 'No'
                disp('Fix your settings!');
        end
                
        
        
    case 'apply_sector'
        %Applies sector and opens plot
        status_flags.analysis_modules.radial_average.sector_mask_chk = 1;
        radial_average_callbacks2('averaging_control', 'azimuthal' );
        
        
     case 'change_ext'        
        %Changes extension name for file saving    
        temp = get(mf_fitter.handles.extension, 'String');
        if not(isempty(temp))            
            mf_fitter.extension = temp;           
        end
        
        
    case 'change_folder'
        %Changes folder name for file saving    
        temp = get(mf_fitter.handles.folder, 'String');
        if not(isempty(temp))            
            mf_fitter.folder = temp;           
        end
        
        
    case 'change_cycs'
        %Changes fit data cycles for plots & data tabulation      
        temp = get(mf_fitter.handles.cycles, 'String')
        temp = str2num(temp)
        if not(isempty(temp))    
            mf_fitter.fit_data.cycles = temp           
        end
        
        
    case 'smoothing_fwhm'
        %Changes fwhm for smoothing  
        temp = get(mf_fitter.handles.smoothing.fwhm, 'String')
        temp = str2num(temp)
        if not(isempty(temp))    
            mf_fitter.smoothing.fwhm = temp           
        end
        
        % Clear old Smoothed Data
        mf_fitter.SmoothedData = [];
            
    case 'smoothing_stepsize'
        %Changes fwhm for smoothing  
        temp = get(mf_fitter.handles.smoothing.step, 'String')
        temp = str2num(temp)
        if not(isempty(temp))    
            mf_fitter.smoothing.step = temp           
        end
        
        % Clear old SmoothedData
        mf_fitter.SmoothedData = [];
        
    case 'save_files'
        %checks save file types from check boxes
        mf_fitter.save.fig = get(mf_fitter.handles.figFile,'Value');
        mf_fitter.save.jpg = get(mf_fitter.handles.jpgFile,'Value');
        mf_fitter.save.eps = get(mf_fitter.handles.epsFile,'Value');
        mf_fitter.save.pdf = get(mf_fitter.handles.pdfFile,'Value');
        
    case 'view_cycles'
        %checks which cycles user wants to view
        mf_fitter.cycleview.view1 = get(mf_fitter.handles.view.view1,'Value');
        mf_fitter.cycleview.view2 = get(mf_fitter.handles.view.view2,'Value');
        mf_fitter.cycleview.view3 = get(mf_fitter.handles.view.view3,'Value');
        mf_fitter.cycleview.view4 = get(mf_fitter.handles.view.view4,'Value');
        mf_fitter.cycleview.view5 = get(mf_fitter.handles.view.view5,'Value');
            
        
end

