function [ ] = mf_GUI_callbacks( to_do, opt1 )
% v 9.2 2/28/2018 E R Louden

% This function contains the callbacks for the main MFF GUI window.
% It can have up to 2 string inputs, to_do designates which callback to
% perform and opt1 handles othe miscellaneous options for the callbacks.
% The callbacks include:
%       apply_sector
%       change_user_input
%       fit_dropdown
%       initialize
%       go
%       plot_dropdown
%       view_cycles

% view params have their own functions - should this be moved?


%% Initialize Global Variables
global mf_fitter
global status_flags
global displayimage
global grasp_handles
global grasp_env

%% Wrong Number of Inputs
% not all callbacks require an option

if(nargin < 2)
    opt1 = 'none';
end

%% Callbacks

switch to_do
    
            
    case 'apply_sector'
        % Applies sector and opens plot
        status_flags.analysis_modules.radial_average.sector_mask_chk = 1;
        radial_average_callbacks2('averaging_control', 'azimuthal' );        

        
    
    case 'change_user_input'
        if(strcmp(opt1,'none'))
           disp('please specifiy which user_input needs to be changed using the opt1 variable') 
        else
            % check if passing more than just one field
            if(strfind(opt1,'.')) 
                opt1a = opt1(1:(strfind(opt1,'.')-1));
                opt1b = opt1((strfind(opt1,'.')+1):length(opt1));
                temp = str2double(get(mf_fitter.handles.(char(opt1a)).(char(opt1b)),'string'));
                if not(isempty(temp)) 
                    mf_fitter.user_inputs.(char(opt1a)).(char(opt1b)) = temp;
                end
            % only one field    
            else
                temp = str2double(get(mf_fitter.handles.(char(opt1)),'string'));
                if not(isempty(temp)) 
                    mf_fitter.user_inputs.(char(opt1)) = temp;
                end
            end
        end

     
        
    case 'fit_dropdown'
        disp('setting fittting algorithm')
        
        % Fitting Paradigms
        if ( get(mf_fitter.handles.fit_dropdown,'Value') == 1 )
            disp('Please select a fitting option')
            
        % Peaks Fixed - initial fitting paradigm, not currenlty used
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
            mf_fitter.handles.window_modules.fit_cycles_4.window = figure(....
                'units','pixels',....
                'Position',fig_position,....
                'Name','View Cycles' ,....
                'NumberTitle', 'off',....
                'Tag','mf_cyc_window',...
                'color',grasp_env.background_color,....
                'menubar','none',....
                'resize','off');
            
            handle = mf_fitter.handles.window_modules.fit_cycles_4.window;
            
            %Question           
            uicontrol(handle,'units','normalized','Position',[0.05 0.55 0.75 0.4],'FontName',grasp_env.font,'FontSize',2*grasp_env.fontsize,'HorizontalAlignment','left','Style','text','String','Which fitting cycles would you like to view?','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
            
            %Checkboxes
            val =  mf_fitter.algorithm_options.cycle_view.view1;
            uicontrol(handle,'units','normalized','Position',[0.2 0.3 0.3 0.25],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','left','Style','text','String','Ref Files?','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
            mf_fitter.handles.window_modules.fit_cycles_4.view1 = uicontrol(handle,'units','normalized','Position',[0.2 0.25 0.05 0.1],'tooltip','Check if you would like to view plots from this fitting cycle','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','checkbox','HorizontalAlignment','right','Tag','angular_bin','Visible','on','Value', val,'Callback','mf_GUI_callbacks(''view_cycles'')');

            val =  mf_fitter.algorithm_options.cycle_view.view2;
            uicontrol(handle,'units','normalized','Position',[0.35 0.3 0.3 0.25],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','left','Style','text','String','Cyc 1?','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
            mf_fitter.mf_fitter.handles.window_modules.fit_cycles_4.view2 = uicontrol(handle,'units','normalized','Position',[0.35 0.25 0.05 0.1],'tooltip','Check if you would like to view plots from this fitting cycle','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','checkbox','HorizontalAlignment','right','Tag','angular_bin','Visible','on','Value', val,'Callback','mf_GUI_callbacks(''view_cycles'')');
        
            val = mf_fitter.algorithm_options.cycle_view.view3;
            uicontrol(handle,'units','normalized','Position',[0.5 0.3 0.3 0.25],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','left','Style','text','String','Cyc 2?','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
            mf_fitter.handles.window_modules.fit_cycles_4.view3 = uicontrol(handle,'units','normalized','Position',[0.5 0.25 0.05 0.1],'tooltip','Check if you would like to view plots from this fitting cycle','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','checkbox','HorizontalAlignment','right','Tag','angular_bin','Visible','on','Value', val,'Callback','mf_GUI_callbacks(''view_cycles'')');
            
            val =  mf_fitter.algorithm_options.cycle_view.view4;
            uicontrol(handle,'units','normalized','Position',[0.65 0.3 0.3 0.25],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','left','Style','text','String','Cyc 3?','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
            mf_fitter.handles.window_modules.fit_cycles_4.view4 = uicontrol(handle,'units','normalized','Position',[0.65 0.25 0.05 0.1],'tooltip','Check if you would like to view plots from this fitting cycle','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','checkbox','HorizontalAlignment','right','Tag','angular_bin','Visible','on','Value', val,'Callback','mf_GUI_callbacks(''view_cycles'')');
           
            val =  mf_fitter.algorithm_options.cycle_view.view5;
            uicontrol(handle,'units','normalized','Position',[0.8 0.3 0.3 0.25],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','left','Style','text','String','Cyc 4?','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
            mf_fitter.handles.window_modules.fit_cycles_4.view5 = uicontrol(handle,'units','normalized','Position',[0.8 0.25 0.05 0.1],'tooltip','Check if you would like to view plots from this fitting cycle','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','checkbox','HorizontalAlignment','right','Tag','angular_bin','Visible','on','Value', val,'Callback','mf_GUI_callbacks(''view_cycles'')');
            
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
        

        
    case 'go'
        % Double check the user (cough cough Liz) has correctly changed all the settings
        choice = questdlg('Have you adjusted all the settings?  experiment name, beam center, fitting algorithm, control parameter, extension name, sector, smoothing','Double Check','Yes','No', 'No');
        
        switch choice         
            case 'Yes'
                disp('setting fit paradigm')

                % Fitting Paradigms
                if ( get(mf_fitter.handles.fit_dropdown,'Value') == 1 )
                    disp('Please select a fitting option')

                % Peaks Fixed
                elseif( get(mf_fitter.handles.fit_dropdown,'Value') == 2 )
                    disp('Peak centers are fixed');
                    % algorithm script

                % 2K Peaks Grow
                elseif( get(mf_fitter.handles.fit_dropdown,'Value') == 3 )
                    disp('ES peaks nucleate in final positions and get larger - fits 3 peaks')
                    mf_fitter_mpf2K();

                % 14K Peaks Rotate
                elseif( get(mf_fitter.handles.fit_dropdown,'Value') == 4 )
                    disp('MS peaks rotate in to final position - fits 2 peaks')
                    mf_fitter_mpf14K();

                % Raster Scan
                elseif( get(mf_fitter.handles.fit_dropdown,'Value') == 5 )
                    disp('Perform a 2D analysis')
                    % algorithm script

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

                % I vs Q fitting
                elseif( get(mf_fitter.handles.fit_dropdown,'Value') == 7 )
                    mf_fitter_mpf_q();

            else
                disp('error')           
            end
        
            case 'No'
                disp('Fix your settings!');
        end
            
        
        
    case 'initialize'
        
        notification = msgbox('Loading files...');
        
        % Before loading, first initialize in case previous run still loaded
        
        % MFF data structure initialization
        mf_fitter.algorithm_options = [];
        mf_fitter.fit_data = [];
        mf_fitter.handles = [];
        mf_fitter.numors = [];
        mf_fitter.save_options = [];
        mf_fitter.user_inputs = [];       
        
        % Default values
        mf_fitter.depth = (status_flags.selector.fdpth_max - 1);
        mf_fitter.user_inputs.int_cutoff = 0.1;
        mf_fitter.user_inputs.dphi_cutoff = 2; 

        
        % Numor array
        for(i = 1:mf_fitter.depth)
        
            status_flags.selector.fd = i+1;  %the +1 is necessary because position 1 corresponds to having all of the data files summed
            grasp_update;
            mf_fitter.numors(i,1) = displayimage.params1(128);

        end
        
        %Set initial reference files as the first and last files
        mf_fitter.user_inputs.reference_files.one = mf_fitter.numors(1,1);
        mf_fitter.user_inputs.reference_files.two = mf_fitter.numors(mf_fitter.depth, 1);
        
        % Set default smoothing values
        mf_fitter.user_inputs.smoothing.fwhm = 1.25;
        mf_fitter.user_inputs.smoothing.step = 1.5;
        
        % Cycle Views
        mf_fitter.algorithm_options.cycle_view.view1 = 0;
        mf_fitter.algorithm_options.cycle_view.view2 = 0;
        mf_fitter.algorithm_options.cycle_view.view3 = 0;
        mf_fitter.algorithm_options.cycle_view.view4 = 0;
        mf_fitter.algorithm_options.cycle_view.view5 = 1;
        
        % Save files
        mf_fitter.save_options.file_type.fig = 0;
        mf_fitter.save_options.file_type.jpg = 0;
        mf_fitter.save_options.file_type.eps = 0;
        mf_fitter.save_options.file_type.pdf = 0;
        
        % old code to play "everything is awesome" - I'm too sentimental to delete this
        %[mf_fitter.awesome.y, mf_fitter.awesome.Fs] = audioread('sound_fx_package.mp3');
        close(notification);
        
        mf_fitter.save_options.extension = [];
        mf_fitter.save_options.folder = [];
    
        
    
    case 'plot_dropdown'
        
        % Plot Options
        if ( get(mf_fitter.handles.plot_dropdown,'Value') == 1 )
            disp('Please select an option')
            
        % FMS
        elseif( get(mf_fitter.handles.plot_dropdown,'Value') == 2 )
            mf_fitter_callbacks('fraction');
            %mf_fitter_plotMaker('fraction')
            
        % Peak Centers
        elseif( get(mf_fitter.handles.plot_dropdown,'Value') == 3 )
            %mf_fitter_callbacks('center_plot');
            mf_fitter_plotMaker('center_plot');
           
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
            % mf_fitter_callbacks('colormap');
            % mf_fitter_plotMaker('colormap');
            
            % code to make the colormap, this does NOT belong here - move to plotMaker
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
           mf_fitter_peakDecay(-20,20); % again, move to plotMaker 
           
        else
            disp('error')           
        end
       

      
  
        
      %%%%  
           
    case 'view_cycles'
        % checks & stores which cycles user wants to view
        mf_fitter.algorithm_options.cycle_view.view1 = get(mf_fitter.handles.window_modules.fit_cycles_4.view1,'Value');
        mf_fitter.algorithm_options.cycle_view.view2 = get(mf_fitter.handles.window_modules.fit_cycles_4.view2,'Value');
        mf_fitter.algorithm_options.cycle_view.view3 = get(mf_fitter.handles.window_modules.fit_cycles_4.view3,'Value');
        mf_fitter.algorithm_options.cycle_view.view4 = get(mf_fitter.handles.window_modules.fit_cycles_4.view4,'Value');
        mf_fitter.algorithm_options.cycle_view.view5 = get(mf_fitter.handles.window_modules.fit_cycles_4.view5,'Value');
            
        
end

