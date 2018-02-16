function [ ] = mf_fitter_NEWcallbacks( to_do, num_peaks, img_num, variable)
% List of callbacks in this script:
    % FitCheck              -   callback to run after every fit cycle
    % MakeFitWindow         -   generates window to pause fitting
    % Refit - Grasp *       -   pulls up 
    % Update *
    
    % * consider moving these callbacks to "mf_fitter_table_callbacks.m"

% v 9.1
% 5/12/2017 MFF Liz

%% Initialize 

    global grasp_handles;
    global grasp_env;
    global mf_fitter;
    global status_flags;

    if(nargin<2)
    elseif(nargin<2)
        num_peaks = 0;
    elseif(nargin<3)
        img_num = 0;
    elseif(nargin<4)
        variable = 0;
    end
    
%% Callbacks
switch to_do

    case 'FitCheck'
        %Pop up window asking if ready to move on
        mf_fitter_NEWcallbacks('MakeFitWindow')
        
        % generate plots to see if things look as expected
        [cm_han, ps_han, int_han, int_tot_han, ad_han] = mf_fitter_plotMaker();
        
        % pull up fit table
        % here user can check that values are ok and pull up questionable fits
        mf_fitter_table;
        
        % nothing else will happen until the continue screen is closed
        waitfor(grasp_handles.window_modules.fit_check.window)

        % close all windows
        handles = [cm_han, ps_han, int_han, int_tot_han, ad_han, mf_fitter.handles.table]
        for h = 1:length(handles)
            try
                close(handles(h))
            catch
                disp('That figure is already closed or plot was never generated.')
            end
        end 
        
        
      
    case 'MakeFitWindow'
        fig_position = [1000, 580, 350, 150];
            grasp_handles.window_modules.fit_check.window = figure(....
                    'units','pixels',....
                    'Position',fig_position,....
                    'Name','Check Fit' ,....
                    'NumberTitle', 'off',....
                    'Tag','mf_cyc_window',...
                    'color',grasp_env.background_color,....
                    'menubar','none',....
                    'resize','off');

        handle = grasp_handles.window_modules.fit_check.window;

        %Question           
        uicontrol(handle,'units','normalized','Position',[0.05 0.55 0.75 0.4],'FontName',grasp_env.font,'FontSize',2*grasp_env.fontsize,'HorizontalAlignment','left','Style','text','String','Please check that fits look reasonable so far.','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
        uicontrol(handle,'units','normalized','Position',[0.05 0.25 0.75 0.4],'FontName',grasp_env.font,'FontSize',2*grasp_env.fontsize,'HorizontalAlignment','left','Style','text','String','Fitting will not continue until this window is closed.','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);

        
        
    case 'refit'    
        % in terms of offering an interactive fitting option, 
        % GRASP is better
        
        status_flags.fitter.number1d = num_peaks;
        status_flags.fitter.function_info_1d.no_parameters = 4;
        status_flags.fitter.function_info_1d.variable_names = repmat({'y0', 'i0', 'xc', 'fwhm'},1,num_peaks);
        status_flags.fitter.function_info_1d.long_names = repmat({'Background', 'Integrated Intensity', 'Center', 'FWHM'},1,num_peaks);
        status_flags.fitter.function_info_1d.group = repmat([1,0,0,1],1,num_peaks);
        status_flags.fitter.function_info_1d.values = repmat([0,0,0,0],1,num_peaks);
        status_flags.fitter.function_info_1d.err_values = repmat([0,0,0,0],1,num_peaks);
        status_flags.fitter.function_info_1d.fix = repmat([1],1,4*num_peaks)
      
        %First give guess values to be overwritten if should be fixed              
        % Check whether matlab or GRASP used for fittings to adjust fit values accordingly     
        if(mf_fitter.fitter.type{mf_fitter.fitter.cycles} == 'M')
            % Data was smoothed and offset using the sector box value
            offset = status_flags.analysis_modules.sectors.theta;
            
            % Matlab only fits one background, GRASP fits one per peak
            divisor = 3;
        else
            % No offset/divisor necessary
            offset = 0;
            divisor = 1;
        end
        
        if(num_peaks~=1)
           status_flags.fitter.function_info_1d.values(1) = mf_fitter.fit_data.background(img_num,1) / divisor;
           status_flags.fitter.function_info_1d.values(4) = mf_fitter.fit_data.fwhm(img_num,1);
           status_flags.fitter.function_info_1d.values(3) = mf_fitter.fit_data.center1(img_num,1)+offset;
           status_flags.fitter.function_info_1d.values(7) = mf_fitter.fit_data.center3(img_num,1)+offset;
           status_flags.fitter.function_info_1d.values(2) = mf_fitter.fit_data.intensity1(img_num,1);
           status_flags.fitter.function_info_1d.values(6) = mf_fitter.fit_data.intensity3(img_num,1);
            if(num_peaks==3)
                status_flags.fitter.function_info_1d.values(7) = mf_fitter.fit_data.center2(img_num,1);
                status_flags.fitter.function_info_1d.values(11) = mf_fitter.fit_data.center3(img_num,1)+offset;
                status_flags.fitter.function_info_1d.values(6) = mf_fitter.fit_data.intensity2(img_num,1);
                status_flags.fitter.function_info_1d.values(10) = mf_fitter.fit_data.intensity3(img_num,1);
            end
        end
           
         grasp_update();
         grasp_plot_fit_callbacks('update_curve_fit_window');
         grasp_plot_fit_window()
         mf_fitter_callbacks('fit',3,img_num)


    case 'update'    
        % save data/errors into data struct
        
        if(mf_fitter.fitter.type{mf_fitter.fitter.cycles} == 'M')
            % Data was smoothed and offset using the sector box value
            offset = status_flags.analysis_modules.sectors.theta;
            
            % Matlab only fits one background, GRASP fits one per peak
            multiplier = 3;
        else
            % No offset/multipier necessary
            offset = 0;
            multiplier = 1;
        end
         
        varNames = {'background', 'intensity1', 'center1', 'fwhm', 'intensity2', 'center2', 'intensity3', 'center3'};
        varNums = [1, 2, 3, 4, 6, 7, 10, 11];
  
        for vN = 1:length(varNames)
            if(strcmp(varNames{vN},'background'))
                mult = multiplier;
                off = 0;
            elseif(strcmp(varNames{vN},'center1') || strcmp(varNames{vN},'center3'))
                mult = 1;
                off = offset;
            else
                mult = 1;
                off = 0;
            end
           mf_fitter.fit_data.(varNames{vN})(img_num,1) = status_flags.fitter.function_info_1d.values(varNums(vN))*mult-off;
           mf_fitter.fit_data.(varNames{vN})(img_num,2) = status_flags.fitter.function_info_1d.err_values(varNums(vN)) 
        end
        
        % regenerate the table
        %mf_fitter_table();
   
end

end