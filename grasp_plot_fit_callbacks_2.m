function output = grasp_plot_fit_callbacks_2(to_do,option)

%Same as grasp_plot-fit_callback, but uses temp var to store chi^2 of fit

% 7/22/14 MFF Allan
% 2/15/2018 MFF v. 9.1 E R Louden

global status_flags
global grasp_handles
global grasp_env
global mf_fitter

if nargin <2; option = []; end

output = [];

%Some useful handles to the current Grasp_Plot
curve_handles = [];
grasp_plot_figure = findobj('tag','grasp_plot');
if not(isempty(grasp_plot_figure));
    grasp_plot_figure = grasp_plot_figure(1);
    grasp_plot_handles = get(grasp_plot_figure,'userdata');
    curve_handles = get(grasp_plot_handles.axis,'userdata');
end


%some editable parameters
guess_color = 'red';
fit_color = 'green';

switch(to_do)
    
    case 'delete_onoff_check'
        status_flags.fitter.delete_curves_check = get(gcbo,'value');
        
    
    case 'average_fit_params'
        
        %Place 1D Curve Fit History on Figure
        n_curves = length(curve_handles);
        
        %Find which curves have fits
        curve_fits = [];
        for n =1:n_curves
            plot_info = get(curve_handles{n}(1),'userdata');
            if isfield(plot_info,'fit1d_parameters'); %then curve has a fit
                curve_fits = [curve_fits, n];
            end
        end
        
        if not(isempty(curve_fits)); %See if any curve fits exist
            
            %Get plot_info for first curve (that has a fit)
            plot_info = get(curve_handles{curve_fits(1)}(1),'userdata');
            n_params = plot_info.fit1d_parameters.no_parameters;
            n_func = plot_info.fit1d_parameters.number_functions;
            pnames = plot_info.fit1d_parameters.long_names;
            fn_name = plot_info.fit1d_parameters.name;
            values_store = []; errors_store = [];
            
            disp(' ')
            disp(['***** Average over Fit Parameters [' num2str(length(curve_fits)) ' curves] *****'])
            disp(['Curve Fit Function is ' fn_name]);
            for n = 1:length(curve_fits)
                plot_info = get(curve_handles{curve_fits(n)}(1),'userdata');
                values_store(n,:) = plot_info.fit1d_parameters.values;
                errors_store(n,:) = plot_info.fit1d_parameters.err_values;
            end
            
            %Make the average of the parameters
            for m = 1:n_func;
                for n = 1:n_params;
                    average = average_error([values_store(:,n*m),errors_store(:,n*m)]);
                    disp(['Average ' pnames{n} ' = ' num2str(average(1)) ' err ' num2str(average(2))]);
                end
            end
            disp('**********************************************************************')
            disp(' ')
            
        else
            disp('No Fit Parameters Availaible to Average')
        end
        
    case 'include_resolution_check'
        status_flags.fitter.include_res_check = get(gcbo,'value');
        grasp_plot_fit_callbacks_2('draw_fn',guess_color); %update the function
        
    case 'curve_number'
        grasp_plot_fit_callbacks_2('draw_fn',guess_color); %update the function
        
    case 'simultaneous_fit_check'
        status_flags.fitter.simultaneous_check = get(gcbo,'value');
        grasp_plot_fit_callbacks_2('draw_fn',guess_color); %update the function
        
        
    case 'auto_guess'
        
        if strcmp(option,'point_click');
            code = status_flags.fitter.function_info_1d.point_click_code;
        else
            code = status_flags.fitter.function_info_1d.auto_guess_code;
            if status_flags.fitter.number1d >1;
                disp(' ');
                disp('Sorry: Cannot Auto-Guess for multiple curves');
                disp('Please enter start parameters manually or use ''Point & Click''');
                disp(' ');
                return
            end
        end
        
        if not(isempty(code));
            
           try
                
                figure(grasp_plot_figure); %Make sure the current grasp_plot is the active figure so mouse input works on that figure
                fitdata = grasp_plot_fit_callbacks_2('get_fit_data');
                x = fitdata.xdat_all; y = fitdata.ydat_all;
                
                final_guess = [];
                for n_fn = 1:status_flags.fitter.number1d
                    for line = 1:length(code)
                        %code{line}
                        eval([code{line} ';'])
                    end
                    if exist('guess_values');
                        final_guess = [final_guess,guess_values];
                    end
                end
                
                %Place guess parameters back into fitter
                %Only place parameters that are NOT fixed
                for n = 1:length(final_guess)
                    if status_flags.fitter.function_info_1d.fix(n) ~=1; %i.e. not fixed
                        status_flags.fitter.function_info_1d.values(n) = final_guess(n);
                        status_flags.fitter.function_info_1d.err_values(n) = 0;
                    end
                end
                
                grasp_plot_fit_callbacks_2('update_curve_fit_window');
                grasp_plot_fit_callbacks_2('draw_fn',guess_color);
                
            catch
                disp(' ')
                disp(['There was an error in evaluating the AutoGuessCode in your Fit Function']);
                disp('Reform your Matlab code and Try again.');
                disp('Make sure that the variable ''guess_values = [    ]'' is defined as a final result of the auto guess procedure');
                disp(' ')
                return
            end
            
        else
            disp(' ');
            disp('Sorry:  No Auto-Guess or Point & Click code is availaible for this function');
            disp('Please enter start parameters manually');
            disp(' ');
        end
        
    case 'copy_to_clipboard'
        
        str = [];
        for n = 1:length(status_flags.fitter.function_info_1d.values)
            str = [str status_flags.fitter.function_info_1d.long_names{n} char(9) num2str(status_flags.fitter.function_info_1d.values(n),'%0.5g') char(9) num2str(status_flags.fitter.function_info_1d.err_values(n),'%0.5g') char(13) char(10)];
        end
        clipboard('copy',str);
        
        
    case 'build_curve_number'
        
        n_str = [];
        for n =1:length(curve_handles);
            if not(isempty(curve_handles{n}));
                plot_info = get(curve_handles{n}(1),'userdata');
                curve_number = plot_info.curve_number;
                n_str = [n_str, {num2str(curve_number)}];
            end
        end
        
        if isfield(grasp_handles.window_modules.curve_fit1d,'curve_number')
            if ishandle(grasp_handles.window_modules.curve_fit1d.curve_number);
                set(grasp_handles.window_modules.curve_fit1d.curve_number,'string',n_str);
                value = get(grasp_handles.window_modules.curve_fit1d.curve_number,'value');
                if value > length(n_str);
                    set(grasp_handles.window_modules.curve_fit1d.curve_number,'value',length(n_str));
                end
                
            end
        end
        
        
    case 'delete_curves'
        for n = 1:length(curve_handles);
            if not(isempty(curve_handles{n}))
                plot_info = get(curve_handles{n}(1),'userdata');
                if isfield(plot_info,'fit_curve_handle');
                    if ishandle(plot_info.fit_curve_handle);
                        delete(plot_info.fit_curve_handle);
                        if isfield(plot_info,'fit1d_history')
                            plot_info = rmfield(plot_info,'fit1d_history');
                            plot_info = rmfield(plot_info,'fit1d_parameters');
                        end
                        set(curve_handles{n}(1),'userdata',plot_info);
                    end
                end
            end
        end
        
        
    case 'build_fn_list'
        
        fn_name_list = [];
        fid=fopen('functions1d.fn');
        if not(fid==-1)
            disp('Loading 1D Fitting List');
            while feof(fid) ==0;
                line = fgetl(fid);
                if strcmpi(line,'<FnName>')
                    fn_name = fgetl(fid);
                    fn_name_list = [fn_name_list, {fn_name}];
                end
                if strcmpi(line,'<Divider>')
                    fn_name = fgetl(fid);
                    fn_name_list = [fn_name_list, {fn_name}];
                end
            end
            fclose(fid);
        else
            disp('No 1D Fitting Functions File Found');
        end
        status_flags.fitter.fn_list1d = fn_name_list; %Keep a store of the fitting functions in the status_flags
        set(grasp_handles.window_modules.curve_fit1d.fn_selector,'string',fn_name_list,'value',status_flags.fitter.fn1d);
        
        
    case 'retrieve_fn'
        
        if isempty(option);
            if status_flags.fitter.fn1d > length(status_flags.fitter.fn_list1d); status_flags.fitter.fn1d = 1;  close(grasp_handles.window_modules.curve_fit1d.window); grasp_plot_fit_window; end
            fn.name = status_flags.fitter.fn_list1d{status_flags.fitter.fn1d};
        else
            fn.name = option;
        end
        fn.variable_names = [];
        fn.long_names = [];
        fn.math_code = [];
        fn.auto_guess_code = [];
        fn.point_click_code = [];
        values = [];
  
        %Open the function file and find the fit parameters and function
        loop_exit = 0;
        fid=fopen('functions1d.fn');
        if not(fid==-1)
            disp(['Loading 1D Fit Function: ' fn.name]);
            if not(isempty(strfind(fn.name,'----'))); %Check for Divider line
                status_flags.fitter.function_info_1d = fn; %Put the blank functions parameters in to status flags
                grasp_plot_fit_callbacks_2('update_curve_fit_window'); %update the window to be blank
                return;
            end
            while feof(fid) ==0 || loop_exit ==0;
                line = fgetl(fid);
                %Search for FnName declarations
                if strcmpi(line,'<Function>')
                    while not(strcmpi(line,'</Function>')) || loop_exit ==0
                        line = fgetl(fid);
                        if strcmpi(line,'<FnName>')
                            line = fgetl(fid);
                            
                            %Search for Specific Fn
                            if strcmpi(line,fn.name)
                                
                                %Search for Parameters
                                while not(strcmpi(line,'</Function>'))
                                    line = fgetl(fid);
                                    
                                    %Read Param variable names
                                    if strcmpi(line,'<Params>')
                                        while not(strcmpi(line,'</Params>'))
                                            line = fgetl(fid);
                                            if not(strcmpi(line,'</Params>')) && not(strcmpi(line(1),'%'))
                                                fn.variable_names = [fn.variable_names, {line}];
                                            end
                                        end
                                    end
                                    
                                    %Read Param Start Values
                                    if strcmpi(line,'<StartValues>')
                                        while not(strcmpi(line,'</StartValues>'))
                                            line = fgetl(fid);
                                            if not(strcmpi(line,'</StartValues>')) && not(strcmpi(line(1),'%'))
                                                values = [values, {line}];
                                            end
                                            %Convert cell of values to numbers
                                            for n =1:length(values)
                                                fn.values(1,n) = str2num(values{n});
                                            end
                                        end
                                        
                                        
                                    end
                                    
                                    %Read Long Names
                                    if strcmpi(line,'<ParamNames>')
                                        while not(strcmpi(line,'</ParamNames>'))
                                            line = fgetl(fid);
                                            if not(strcmpi(line,'</ParamNames>')) && not(strcmpi(line(1),'%'))
                                                fn.long_names = [fn.long_names, {line}];
                                            end
                                        end
                                    end
                                    
                                    %Read the math function string
                                    if strcmpi(line,'<FnCode>')
                                        while not(strcmpi(line,'</FnCode>'))
                                            line = fgetl(fid);
                                            if not(strcmpi(line,'</FnCode>')) && not(strcmpi(line(1),'%'))
                                                fn.math_code = [fn.math_code, {line}];
                                            end
                                        end
                                    end
                                    
                                    %Read the AutoGuess Code
                                    if strcmpi(line,'<AutoGuessCode>')
                                        while not(strcmpi(line,'</AutoGuessCode>'))
                                            line = fgetl(fid);
                                            if not(strcmpi(line,'</AutoGuessCode>')) && not(strcmpi(line(1),'%'))
                                                fn.auto_guess_code = [fn.auto_guess_code, {line}];
                                            end
                                        end
                                    end
                                    
                                    %Read the Point Click parameter guide
                                    if strcmpi(line,'<PointClickCode>')
                                        while not(strcmpi(line,'</PointClickCode>'))
                                            line = fgetl(fid);
                                            if not(strcmpi(line,'</PointClickCode>')) && not(strcmpi(line(1),'%'))
                                                fn.point_click_code = [fn.point_click_code, {line}];
                                            end
                                        end
                                    end
                                end
                                loop_exit =1;
                            end
                        end
                    end
                end
            end
            fclose(fid);
        else
            disp('No 1D Fitting Functions File Found');
        end
        fn.fix = zeros(1,length(fn.variable_names)); %Fix parameter
        fn.group = zeros(1,length(fn.variable_names)); %Group parameter
        fn.err_values = zeros(1,length(fn.values)); %Err_Values
        fn.no_parameters = length(fn.values); %Intrinsic number of parameters BEFORE multiplexing
        %Put the fn into the status flags
        status_flags.fitter.function_info_1d = fn; %some of these get added to in the multiplex process below
        
        %***** Function Multiplex *****
        %Multiply up the function i.e. n of the same functions, with
        %different or grouped parameters added together
        for n = 2:status_flags.fitter.number1d;
            status_flags.fitter.function_info_1d.fix = [status_flags.fitter.function_info_1d.fix, fn.fix];
            status_flags.fitter.function_info_1d.group = [status_flags.fitter.function_info_1d.group, fn.group];
            status_flags.fitter.function_info_1d.values = [status_flags.fitter.function_info_1d.values, fn.values];
            status_flags.fitter.function_info_1d.err_values = [status_flags.fitter.function_info_1d.err_values, fn.err_values];
            status_flags.fitter.function_info_1d.variable_names = [status_flags.fitter.function_info_1d.variable_names, fn.variable_names];
            status_flags.fitter.function_info_1d.long_names = [status_flags.fitter.function_info_1d.long_names, fn.long_names];
            %status_flags.fitter.function_info_1d.point_click_code = [status_flags.fitter.function_info_1d.point_click_code, fn.point_click_code];
        end
        
        %Also allow output of the fit function (e.g. for Ancos2)
        output = status_flags.fitter.function_info_1d;
        
        
    case 'update_curve_fit_window'
        
        %Delete old displayed parameters
        if isfield(grasp_handles,'fitter');
            if isfield(grasp_handles.fitter,'fit1d_handles');
                for n = 1:length(grasp_handles.fitter.fit1d_handles);
                    if ishandle(grasp_handles.fitter.fit1d_handles(n)); delete(grasp_handles.fitter.fit1d_handles(n)); end
                end
            end
        end
        %Draw parameter entry boxes and check-boxes - dependent on the number and parameter names returned by the function
        grasp_handles.fitter.fit1d_handles = [];
        handles_store = [];
        for n = 1:length(status_flags.fitter.function_info_1d.long_names)
            %Parameter Names
            handle = uicontrol(grasp_handles.window_modules.curve_fit1d.window, 'units','normalized','Position',[0.02,(0.8-(n*0.03)),0.3, 0.028],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','right','Style','text','String',[status_flags.fitter.function_info_1d.long_names{n} ': ' status_flags.fitter.function_info_1d.variable_names{n} ' : ' ],'BackgroundColor', grasp_env.sub_figure_background_color, 'ForegroundColor', [1 1 1]);
            handles_store = [handles_store, handle];
            %Fix Check
            handle = uicontrol(grasp_handles.window_modules.curve_fit1d.window, 'units','normalized','Position',[0.37,(0.8-(n*0.03)),0.038,0.028],'ToolTip',['Fix ' status_flags.fitter.function_info_1d.long_names{n}],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','checkbox','HorizontalAlignment','left','backgroundcolor',grasp_env.sub_figure_background_color,'value',status_flags.fitter.function_info_1d.fix(n),'userdata',n,'callback','grasp_plot_fit_callbacks_2(''fn_value_fix_check'');');
            handles_store = [handles_store, handle];
            
            if status_flags.fitter.function_info_1d.group(n) == 1 && n > status_flags.fitter.function_info_1d.no_parameters;
                visible = 'off';
            else
                visible = 'on';
            end
            %Parameter Value
            handle = uicontrol(grasp_handles.window_modules.curve_fit1d.window, 'visible',visible, 'units','normalized','Position',[0.45,(0.8-(n*0.03)),0.15,0.028],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','edit','HorizontalAlignment','left','string',num2str(status_flags.fitter.function_info_1d.values(n)),'userdata',n,'callback','grasp_plot_fit_callbacks_2(''fn_value_enter'');');
            handles_store = [handles_store, handle];
            %Parameter Value Error
            handle = uicontrol(grasp_handles.window_modules.curve_fit1d.window, 'visible',visible, 'units','normalized','Position',[0.62,(0.8-(n*0.03)),0.15,0.028],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','text','HorizontalAlignment','left','string',num2str(status_flags.fitter.function_info_1d.err_values(n)));
            handles_store = [handles_store, handle];
            if n <= status_flags.fitter.function_info_1d.no_parameters;
                %Group Check
                handle = uicontrol(grasp_handles.window_modules.curve_fit1d.window, 'units','normalized','Position',[0.81,(0.8-(n*0.03)),0.038,0.028],'ToolTip',['Group ' status_flags.fitter.function_info_1d.long_names{n} '`s'],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','checkbox','HorizontalAlignment','left','backgroundcolor',grasp_env.sub_figure_background_color,'value',status_flags.fitter.function_info_1d.group(n),'userdata',n,'callback','grasp_plot_fit_callbacks_2(''fn_group_check'');');
                handles_store = [handles_store, handle];
            end
        end
        handle = uicontrol(grasp_handles.window_modules.curve_fit1d.window, 'units','normalized','Position',[0.02,0.1,0.96,0.05],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','text','HorizontalAlignment','left','backgroundcolor',grasp_env.sub_figure_background_color,'foregroundcolor','white','string',status_flags.fitter.function_info_1d.math_code);
        handles_store = [handles_store, handle];
        
        grasp_handles.fitter.fit1d_handles = handles_store;
        
        
    case 'draw_fn'
        %Last check to see that there is a plot open with a valid curve
        if isempty(curve_handles); return; end

        %Get the xdata over the zoomed range
        fitdata = grasp_plot_fit_callbacks_2('get_fit_data'); %fitdata.curve_number contains the curve number
        
        %Loop though the number of simultaneous curves
        for n = 1:fitdata.n_curves
            
            %Evaluate the function in higher resolution (1000 points) if needed
            curve_res = 1000; %High-resolution re-drawing of function
            if length(fitdata.xdat.(['curve' num2str(n)]))<curve_res;
                
                xmin = min(fitdata.xdat.(['curve' num2str(n)])); xmax = max(fitdata.xdat.(['curve' num2str(n)]));
                xspan = xmax-xmin;
                %Make the x-data in higher resolution
                x = xmin:xspan/curve_res:xmax;
                
                %Build the resolution data in higher resolution
                %exdat - simple x-resolution
                if not(isempty(fitdata.exdat.(['curve' num2str(n)])))
                    [temp, i] = unique(fitdata.xdat.(['curve' num2str(n)]));
                    extra.exdat_all = interp1(fitdata.xdat.(['curve' num2str(n)])(i),fitdata.exdat.(['curve' num2str(n)])(i),x);
                else extra.exdat_all = [];
                end
                
                %Build new resolution_kernels
                
                %Need to first interpolate all the sigma's for the various resolution components, lambda, theta, pixel, binning
                if isfield(fitdata,'resolution_kernels')
                    kernel_data = fitdata.resolution_kernels.(['curve' num2str(n)]); %Defaults before interpolation below
                    [temp, i] = unique(fitdata.xdat.(['curve' num2str(n)]));
                    
                    %re-find kernel resolution (assuming it doesn't exist in the kernel_data)
                    k_res = numel(fitdata.resolution_kernels.(['curve' num2str(n)]).x{1});
                    k_mid_index = ceil(k_res/2);

                    if isfield(kernel_data,'x');
                        %Loop though all the high-res x's and interpolate the already built x & weight kernel from the resolution kernel
                        for m = 1:length(fitdata.resolution_kernels.(['curve' num2str(n)]).x);
                            kernel_x_centres(m) = fitdata.resolution_kernels.(['curve' num2str(n)]).x{m}(k_mid_index);
                            kernel_x_matrix(:,m) = fitdata.resolution_kernels.(['curve' num2str(n)]).x{m};
                            kernel_weight_matrix(:,m) = fitdata.resolution_kernels.(['curve' num2str(n)]).weight{m};
                        end
                        [temp,j] = unique(kernel_x_centres);
                        for m = 1:k_res;
                            kernel_x_matrix_new(m,:) = interp1(kernel_x_centres(j),kernel_x_matrix(m,j),x);
                            kernel_weight_matrix_new(m,:) = interp1(kernel_x_centres(j),kernel_weight_matrix(m,j),x);
                        end
                        
                        %Rebuild into the correct form for kernel_data - now interpolated kernels in high res for drawing
                        for m = 1:length(x);
                            kernel_data.x{m} = rot90(kernel_x_matrix_new(:,m),3);
                            kernel_data.weight{m} = rot90(kernel_weight_matrix_new(:,m),3);
                        end
                    end


                    
                    if isfield(kernel_data,'lambda');
                        kernel_data.lambda.fwhm = interp1(fitdata.xdat.(['curve' num2str(n)])(i),fitdata.resolution_kernels.(['curve' num2str(n)]).lambda.fwhm(i),x);
                    end
                    if isfield(kernel_data,'theta');
                        kernel_data.theta.fwhm = interp1(fitdata.xdat.(['curve' num2str(n)])(i),fitdata.resolution_kernels.(['curve' num2str(n)]).theta.fwhm(i),x);
                    end
                    if isfield(kernel_data,'pixel');
                        kernel_data.pixel.fwhm = interp1(fitdata.xdat.(['curve' num2str(n)])(i),fitdata.resolution_kernels.(['curve' num2str(n)]).pixel.fwhm(i),x);
                    end
                    if isfield(kernel_data,'fwhm');
                        kernel_data.binning.fwhm = interp1(fitdata.xdat.(['curve' num2str(n)])(i),fitdata.resolution_kernels.(['curve' num2str(n)]).binning.fwhm(i),x);
                    end
                    if isfield(kernel_data,'classic_res');
                        kernel_data.classic_res.fwhm = interp1(fitdata.xdat.(['curve' num2str(n)])(i),fitdata.resolution_kernels.(['curve' num2str(n)]).classic_res.fwhm(i),x);
                    end
                    if isfield(kernel_data,'aperture');
                        kernel_data.aperture.fwhm = interp1(fitdata.xdat.(['curve' num2str(n)])(i),fitdata.resolution_kernels.(['curve' num2str(n)]).aperture.fwhm(i),x);
                    end
                    
                    %resolution_kernels = build_resolution_kernels(x, kernel_data);
                    %extra.resolution_kernels_all = resolution_kernels;
                    extra.resolution_kernels_all = kernel_data;
                end
            else
                temp = [fitdata.xdat.(['curve' num2str(n)]), fitdata.exdat.(['curve' num2str(n)])];
                temp = sort(temp,1);
                x = temp(:,1); extra.exdat_all = temp(:,2);
                extra.resolution_kernels = fitdata.resolution_kernels.(['curve' num2str(n)]);
            end
 
            %Generate function
            y = pseudo_fn(x,status_flags.fitter.function_info_1d.values,extra);
                
            
            %Delete old guess curve
            plot_info = get(curve_handles{fitdata.curve_number(n)}(1),'userdata');
            if status_flags.fitter.delete_curves_check == 0;
                if isfield(plot_info,'fit_curve_handle')
                    if ishandle(plot_info.fit_curve_handle); delete(plot_info.fit_curve_handle); end
                end
            end
            
            set(grasp_plot_handles.axis,'nextplot','add') %Same as hold on
            handle = plot(grasp_plot_handles.axis, rot90(x,3),y,'color',option);
                        
            %Store the curve fit handle & Fit Curve Data in the curve info.
            if status_flags.fitter.delete_curves_check == 0;
                plot_info.fit_curve_handle = handle;
            else
                if isfield(plot_info,'fit_curve_handle');
                    plot_info.fit_curve_handle = [handle, plot_info.fit_curve_handle];
                else
                    plot_info.fit_curve_handle = handle;
                end
            end
            plot_info.fit1d_curve_data = [rot90(x,3),y];
            set(curve_handles{fitdata.curve_number(n)}(1),'userdata',plot_info);
        end
        
        
    case 'toggle_fn'
        status_flags.fitter.fn1d = get(gcbo,'value');
        grasp_plot_fit_callbacks_2('retrieve_fn');
        grasp_plot_fit_callbacks_2('update_curve_fit_window');
        
        
        
    case 'fn_value_enter'
        param_number = get(gcbo,'userdata');
        value_str = get(gcbo,'string');
        value = str2num(value_str);
        if not(isempty(value))
            status_flags.fitter.function_info_1d.values(1,param_number) = value;
        end
        grasp_plot_fit_callbacks_2('draw_fn',guess_color); %update the function
        
    case 'fn_value_fix_check'
        param_number = get(gcbo,'userdata');
        status_flags.fitter.function_info_1d.fix(param_number) = get(gcbo,'value');
        
    case 'fn_group_check'
        param_number = get(gcbo,'userdata');
        status_flags.fitter.function_info_1d.group(param_number) = get(gcbo,'value');
        while param_number <= length(status_flags.fitter.function_info_1d.group);
            status_flags.fitter.function_info_1d.group(param_number) = get(gcbo,'value');
            %also set the fix check for the grouped copies
            if param_number > status_flags.fitter.function_info_1d.no_parameters;
                status_flags.fitter.function_info_1d.fix(param_number) = 1;
            end
            param_number = param_number + status_flags.fitter.function_info_1d.no_parameters;
        end
        grasp_plot_fit_callbacks_2('update_curve_fit_window');
        
    case 'number_of_fn'
        status_flags.fitter.number1d = get(gcbo,'value');
        grasp_plot_fit_callbacks_2('retrieve_fn');
        grasp_plot_fit_callbacks_2('update_curve_fit_window');
        
    case 'get_fit_data'
        %Check if Plot (data) exist at all
        if isempty(curve_handles); return; end
        
        if status_flags.fitter.simultaneous_check == 1; %fit all curves
            n_curves = length(curve_handles);
            curve = 1:n_curves;
            %Make empty arrays
            for n = 1:n_curves
                xdat.(['curve' num2str(n)]) = [];
                ydat.(['curve' num2str(n)]) = [];
                edat.(['curve' num2str(n)]) = [];
                exdat.(['curve' num2str(n)]) = [];
                resolution_kernels.(['curve' num2str(n)]) = [];
            end
            
            %Loop though all the curves and make individual lists of all the data
            for n = 1:n_curves;
                plot_info = get(curve_handles{n}(1),'userdata');
                xdat.(['curve' num2str(n)]) = [xdat.(['curve' num2str(n)]); plot_info.plot_data.xdat];
                ydat.(['curve' num2str(n)]) = [ydat.(['curve' num2str(n)]); plot_info.plot_data.ydat];
                if isfield(plot_info.plot_data,'edat') && not(isempty(plot_info.plot_data.edat));
                    edat.(['curve' num2str(n)]) = [edat.(['curve' num2str(n)]); plot_info.plot_data.edat];
                else
                    edat.(['curve' num2str(n)]) = [edat.(['curve' num2str(n)]); zeros(size(ydat.(['curve' num2str(n)])))];
                end
                if isfield(plot_info.plot_data,'exdat');
                    exdat.(['curve' num2str(n)]) = [exdat.(['curve' num2str(n)]); plot_info.plot_data.exdat];
                end
                if isfield(plot_info.plot_data,'resolution_kernels');
                    resolution_kernels.(['curve' num2str(n)]) = plot_info.plot_data.resolution_kernels;
                end
            end
        else
            n_curves = 1;
            
            %Make empty arrays
            xdat.(['curve1']) = [];
            ydat.(['curve1']) = [];
            edat.(['curve1']) = [];
            exdat.(['curve1']) = [];
            resolution_kernels.(['curve1']) = [];
            
            %Find which curve from the curve number selector
            value = get(grasp_handles.window_modules.curve_fit1d.curve_number,'value');
            str = get(grasp_handles.window_modules.curve_fit1d.curve_number,'string');
            curve = str2num(str{value});
            
            plot_info = get(curve_handles{curve}(1),'userdata');
            xdat.(['curve1']) = plot_info.plot_data.xdat;
            ydat.(['curve1']) = plot_info.plot_data.ydat;
            if isfield(plot_info.plot_data,'edat') && not(isempty(plot_info.plot_data.edat));
                edat.(['curve1']) = plot_info.plot_data.edat;
            else
                edat.(['curve1']) = zeros(size(ydat));
            end
            if isfield(plot_info.plot_data,'exdat');
                exdat.(['curve1']) = plot_info.plot_data.exdat;
            end
            if isfield(plot_info.plot_data,'resolution_kernels');
                resolution_kernels.(['curve1']) = plot_info.plot_data.resolution_kernels;
            end
        end
        
        %Chop the fit data so as only to use data over the zoomed area
        xlim = get(grasp_plot_handles.axis,'xlim');
        for n = 1:n_curves;
            data_mask = xdat.(['curve' num2str(n)]) > xlim(1) & xdat.(['curve' num2str(n)]) < xlim(2);
            temp = find(data_mask);
            xdat.(['curve' num2str(n)]) = xdat.(['curve' num2str(n)])(temp);
            ydat.(['curve' num2str(n)]) = ydat.(['curve' num2str(n)])(temp);
            edat.(['curve' num2str(n)]) = edat.(['curve' num2str(n)])(temp);
            if not(isempty(exdat.(['curve' num2str(n)])));
                exdat.(['curve' num2str(n)]) = exdat.(['curve' num2str(n)])(temp);
            end
            
            if not(isempty(resolution_kernels.(['curve' num2str(n)])));
                if isfield(resolution_kernels.(['curve' num2str(n)]),'x');
                    resolution_kernels.(['curve' num2str(n)]).x = resolution_kernels.(['curve' num2str(n)]).x(temp);
                end
                if isfield(resolution_kernels.(['curve' num2str(n)]),'weight');
                    resolution_kernels.(['curve' num2str(n)]).weight = resolution_kernels.(['curve' num2str(n)]).weight(temp);
                end
                if isfield(resolution_kernels.(['curve' num2str(n)]),'fwmh');
                    resolution_kernels.(['curve' num2str(n)]).fwhm = resolution_kernels.(['curve' num2str(n)]).fwhm(temp);
                end
                if isfield(resolution_kernels.(['curve' num2str(n)]),'lambda');
                    resolution_kernels.(['curve' num2str(n)]).lambda.fwhm = resolution_kernels.(['curve' num2str(n)]).lambda.fwhm(temp);
                end
                if isfield(resolution_kernels.(['curve' num2str(n)]),'theta');
                    resolution_kernels.(['curve' num2str(n)]).theta.fwhm = resolution_kernels.(['curve' num2str(n)]).theta.fwhm(temp);
                end
                if isfield(resolution_kernels.(['curve' num2str(n)]),'pixel');
                    resolution_kernels.(['curve' num2str(n)]).pixel.fwhm = resolution_kernels.(['curve' num2str(n)]).pixel.fwhm(temp);
                end
                if isfield(resolution_kernels.(['curve' num2str(n)]),'binning');
                    resolution_kernels.(['curve' num2str(n)]).binning.fwhm = resolution_kernels.(['curve' num2str(n)]).binning.fwhm(temp);
                end
                if isfield(resolution_kernels.(['curve' num2str(n)]),'aperture');
                    resolution_kernels.(['curve' num2str(n)]).aperture.fwhm = resolution_kernels.(['curve' num2str(n)]).aperture.fwhm(temp);
                end
                if isfield(resolution_kernels.(['curve' num2str(n)]),'classic_res');
                    resolution_kernels.(['curve' num2str(n)]).classic_res.fwhm = resolution_kernels.(['curve' num2str(n)]).classic_res.fwhm(temp);
                end
                if isfield(resolution_kernels.(['curve' num2str(n)]),'cm');
                    resolution_kernels.(['curve' num2str(n)]).cm = resolution_kernels.(['curve' num2str(n)]).cm;
                end
            end
        end

        %***** Check for any zero counts or errors coming though.  These should be allocated an error of 1 in order to fit correctly. *****
        for n = 1:n_curves
            %temp = find(ydat.(['curve' num2str(n)])==0);
            %if not(isempty(temp)); ydat.(['curve' num2str(n)])(temp) = 1; end
            
            %find max of all errors in data
            %max_error = max(edat.(['curve' num2str(n)]))
            min_error = min(edat.(['curve' num2str(n)])(edat.(['curve' num2str(n)])>0));
            temp = find(edat.(['curve' num2str(n)])==0);
            if not(isempty(temp)); edat.(['curve' num2str(n)])(temp) = min_error; end
        end
        
        %Prevent NaN's coming though (from divide by zeros in data), Delete them from the fit data
        for n = 1:n_curves
            temp = find(not(isnan(ydat.(['curve' num2str(n)]))));
            if not(isempty(temp));
                xdat.(['curve' num2str(n)]) = xdat.(['curve' num2str(n)])(temp);
                ydat.(['curve' num2str(n)]) = ydat.(['curve' num2str(n)])(temp);
                edat.(['curve' num2str(n)]) = edat.(['curve' num2str(n)])(temp);
                if not(isempty(exdat.(['curve' num2str(n)])))
                    exdat.(['curve' num2str(n)]) = exdat.(['curve' num2str(n)])(temp);
                end
                if not(isempty(resolution_kernels.(['curve' num2str(n)])));
                    if isfield(resolution_kernels.(['curve' num2str(n)]),'x');
                    resolution_kernels.(['curve' num2str(n)]).x = resolution_kernels.(['curve' num2str(n)]).x(temp);
                    end
                    if isfield(resolution_kernels.(['curve' num2str(n)]),'weight');
                        resolution_kernels.(['curve' num2str(n)]).weight = resolution_kernels.(['curve' num2str(n)]).weight(temp);
                    end
                    if isfield(resolution_kernels.(['curve' num2str(n)]),'fwmh');
                        resolution_kernels.(['curve' num2str(n)]).fwhm = resolution_kernels.(['curve' num2str(n)]).fwhm(temp);
                    end
                    if isfield(resolution_kernels.(['curve' num2str(n)]),'lambda');
                        resolution_kernels.(['curve' num2str(n)]).lambda.fwhm = resolution_kernels.(['curve' num2str(n)]).lambda.fwhm(temp);
                    end
                    if isfield(resolution_kernels.(['curve' num2str(n)]),'theta');
                        resolution_kernels.(['curve' num2str(n)]).theta.fwhm = resolution_kernels.(['curve' num2str(n)]).theta.fwhm(temp);
                    end
                    if isfield(resolution_kernels.(['curve' num2str(n)]),'pixel');
                        resolution_kernels.(['curve' num2str(n)]).pixel.fwhm = resolution_kernels.(['curve' num2str(n)]).pixel.fwhm(temp);
                    end
                    if isfield(resolution_kernels.(['curve' num2str(n)]),'binning');
                        resolution_kernels.(['curve' num2str(n)]).binning.fwhm = resolution_kernels.(['curve' num2str(n)]).binning.fwhm(temp);
                    end
                    if isfield(resolution_kernels.(['curve' num2str(n)]),'aperture');
                        resolution_kernels.(['curve' num2str(n)]).aperture.fwhm = resolution_kernels.(['curve' num2str(n)]).aperture.fwhm(temp);
                    end
                    if isfield(resolution_kernels.(['curve' num2str(n)]),'classic_res');
                        resolution_kernels.(['curve' num2str(n)]).classic_res.fwhm = resolution_kernels.(['curve' num2str(n)]).classic_res.fwhm(temp);
                    end
                end
            end
        end
        
        % Check number of data points is not less than 10
        for n = 1:n_curves
            no_points = length(xdat.(['curve' num2str(n)]));
            if no_points <10;
                disp(' ');
                disp('Warning:  Number of data points to fit is less than 10');
                disp(['Warning:  Fit data array (curve #' num2str(n) ') will be padded using 0 counts with 10e99 errors']);
                disp(' ');
                need_points = 20-no_points;
                xadd = zeros(need_points,1);
                yadd = zeros(need_points,1);
                exadd = ones(need_points,1)*1e99;
                erradd = ones(need_points,1)*1e99;
                xdat.(['curve' num2str(n)]) = [xdat.(['curve' num2str(n)]);xadd];
                ydat.(['curve' num2str(n)]) = [ydat.(['curve' num2str(n)]);yadd];
                edat.(['curve' num2str(n)]) = [edat.(['curve' num2str(n)]);erradd];
                if not(isempty(exdat.(['curve' num2str(n)])))
                    exdat.(['curve' num2str(n)]) = [exdat.(['curve' num2str(n)]); exadd];
                end
                disp('Help:  What to do in grasp_plot_fit_callbacks_2')
                disp('Need to pad-out the resolution_kernel also??')
            end
        end

        %Prepare output data & Sort into ascending x order (q)
        for n = 1:n_curves
            %Individual curve data for drawing seperately etc.
            [temp, i] = sort(xdat.(['curve' num2str(n)])); %Generate the sort intex, i
            %[i] = find(xdat.(['curve' num2str(n)])); %Leave unsorted
            output.xdat.(['curve' num2str(n)]) = xdat.(['curve' num2str(n)])(i);
            output.ydat.(['curve' num2str(n)]) = ydat.(['curve' num2str(n)])(i);
            output.edat.(['curve' num2str(n)]) = edat.(['curve' num2str(n)])(i);
            if not(isempty(exdat.(['curve' num2str(n)])));
                output.exdat.(['curve' num2str(n)]) = exdat.(['curve' num2str(n)])(i);
            else
                output.exdat.(['curve' num2str(n)]) = [];
            end
            if not(isempty(resolution_kernels.(['curve' num2str(n)])));
                if isfield(resolution_kernels.(['curve' num2str(n)]),'x');
                    output.resolution_kernels.(['curve' num2str(n)]).x = resolution_kernels.(['curve' num2str(n)]).x(i);
                end
                if isfield(resolution_kernels.(['curve' num2str(n)]),'weight');
                    output.resolution_kernels.(['curve' num2str(n)]).weight = resolution_kernels.(['curve' num2str(n)]).weight(i);
                end
                if isfield(resolution_kernels.(['curve' num2str(n)]),'fwhm');
                    output.resolution_kernels.(['curve' num2str(n)]).fwhm = resolution_kernels.(['curve' num2str(n)]).fwhm(i);
                end
                if isfield(resolution_kernels.(['curve' num2str(n)]),'lambda');
                    output.resolution_kernels.(['curve' num2str(n)]).lambda.fwhm = resolution_kernels.(['curve' num2str(n)]).lambda.fwhm(i);
                    output.resolution_kernels.(['curve' num2str(n)]).lambda.shape = resolution_kernels.(['curve' num2str(n)]).lambda.shape;
                end
                if isfield(resolution_kernels.(['curve' num2str(n)]),'theta');
                    output.resolution_kernels.(['curve' num2str(n)]).theta.fwhm = resolution_kernels.(['curve' num2str(n)]).theta.fwhm(i);
                    output.resolution_kernels.(['curve' num2str(n)]).theta.shape = resolution_kernels.(['curve' num2str(n)]).theta.shape;
                end
                if isfield(resolution_kernels.(['curve' num2str(n)]),'pixel');
                    output.resolution_kernels.(['curve' num2str(n)]).pixel.fwhm = resolution_kernels.(['curve' num2str(n)]).pixel.fwhm(i);
                    output.resolution_kernels.(['curve' num2str(n)]).pixel.shape = resolution_kernels.(['curve' num2str(n)]).pixel.shape;
                end
                if isfield(resolution_kernels.(['curve' num2str(n)]),'binning');
                    output.resolution_kernels.(['curve' num2str(n)]).binning.fwhm = resolution_kernels.(['curve' num2str(n)]).binning.fwhm(i);
                    output.resolution_kernels.(['curve' num2str(n)]).binning.shape = resolution_kernels.(['curve' num2str(n)]).binning.shape;
                end
                if isfield(resolution_kernels.(['curve' num2str(n)]),'aperture');
                    output.resolution_kernels.(['curve' num2str(n)]).aperture.fwhm = resolution_kernels.(['curve' num2str(n)]).aperture.fwhm(i);
                    output.resolution_kernels.(['curve' num2str(n)]).aperture.shape = resolution_kernels.(['curve' num2str(n)]).aperture.shape;
                end
                if isfield(resolution_kernels.(['curve' num2str(n)]),'theta');
                    output.resolution_kernels.(['curve' num2str(n)]).classic_res.fwhm = resolution_kernels.(['curve' num2str(n)]).classic_res.fwhm(i);
                    output.resolution_kernels.(['curve' num2str(n)]).classic_res.shape = resolution_kernels.(['curve' num2str(n)]).classic_res.shape;
                end
                if isfield(resolution_kernels.(['curve' num2str(n)]),'cm');
                    output.resolution_kernels.(['curve' num2str(n)]).cm = resolution_kernels.(['curve' num2str(n)]).cm;
                end
                if isfield(resolution_kernels.(['curve' num2str(n)]),'fwhmwidth');
                    output.resolution_kernels.(['curve' num2str(n)]).fwhmwidth = resolution_kernels.(['curve' num2str(n)]).fwhmwidth;
                end
                if isfield(resolution_kernels.(['curve' num2str(n)]),'finesse');
                    output.resolution_kernels.(['curve' num2str(n)]).finesse = resolution_kernels.(['curve' num2str(n)]).finesse;
                end
                if isfield(resolution_kernels.(['curve' num2str(n)]),'history');
                    output.resolution_kernels.(['curve' num2str(n)]).history = resolution_kernels.(['curve' num2str(n)]).history;
                end
            end
        end
        
        %Make an appended list of all fit data for simultaneous fitting of all curves
        output.xdat_all = []; output.ydat_all = []; output.edat_all = []; output.exdat_all = []; 
        %if not(isempty(exdat.(['curve' num2str(n)]))); output.exdat_all = []; end
        for n = 1:n_curves
            %temp = find(output.edat.(['curve' num2str(n)])~=1e99); %do not include the padded points in the list
            %IF USING THIS CODE THEN ALSO SWAP BELOW FOR EXDAT_ALL
            %AND Swap the large block between the '(temp)' or not.
            %output.xdat_all = [output.xdat_all; output.xdat.(['curve' num2str(n)])(temp)];
            %output.ydat_all = [output.ydat_all; output.ydat.(['curve' num2str(n)])(temp)];
            %output.edat_all = [output.edat_all; output.edat.(['curve' num2str(n)])(temp)];
            
            %instead of above commented, use all points, included padded points
            output.xdat_all = [output.xdat_all; output.xdat.(['curve' num2str(n)])];
            output.ydat_all = [output.ydat_all; output.ydat.(['curve' num2str(n)])];
            output.edat_all = [output.edat_all; output.edat.(['curve' num2str(n)])];
            

            
            if not(isempty(exdat.(['curve' num2str(n)])));
                %output.exdat_all = [output.exdat_all; output.exdat.(['curve' num2str(n)])(temp)];
                output.exdat_all = [output.exdat_all; output.exdat.(['curve' num2str(n)])];
            else
                output.exdat_all = [];
            end
            
            if not(isempty(resolution_kernels.(['curve' num2str(n)])));
                if n ==1;
                    output.resolution_kernels_all = output.resolution_kernels.(['curve' num2str(n)]);
                else
                    output.resolution_kernels_all.x = [output.resolution_kernels_all.x output.resolution_kernels.(['curve' num2str(n)]).x];
                    output.resolution_kernels_all.weight = [output.resolution_kernels_all.weight output.resolution_kernels.(['curve' num2str(n)]).weight];
                    output.resolution_kernels_all.fwhm = [output.resolution_kernels_all.fwhm; output.resolution_kernels.(['curve' num2str(n)]).fwhm];
                    output.resolution_kernels_all.lambda.fwhm = [output.resolution_kernels_all.lambda.fwhm; output.resolution_kernels.(['curve' num2str(n)]).lambda.fwhm];
                    output.resolution_kernels_all.lambda.shape = output.resolution_kernels.(['curve' num2str(n)]).lambda.shape;
                    output.resolution_kernels_all.theta.fwhm = [output.resolution_kernels_all.theta.fwhm; output.resolution_kernels.(['curve' num2str(n)]).theta.fwhm];
                    output.resolution_kernels_all.theta.shape = output.resolution_kernels.(['curve' num2str(n)]).theta.shape;
                    output.resolution_kernels_all.pixel.fwhm = [output.resolution_kernels_all.pixel.fwhm; output.resolution_kernels.(['curve' num2str(n)]).pixel.fwhm];
                    output.resolution_kernels_all.pixel.shape = output.resolution_kernels.(['curve' num2str(n)]).pixel.shape;
                    output.resolution_kernels_all.binning.fwhm = [output.resolution_kernels_all.binning.fwhm; output.resolution_kernels.(['curve' num2str(n)]).binning.fwhm];
                    output.resolution_kernels_all.binning.shape = output.resolution_kernels.(['curve' num2str(n)]).binning.shape;
                    output.resolution_kernels_all.aperture.fwhm = [output.resolution_kernels_all.aperture.fwhm; output.resolution_kernels.(['curve' num2str(n)]).binning.fwhm];
                    output.resolution_kernels_all.aperture.shape = output.resolution_kernels.(['curve' num2str(n)]).aperture.shape;
                    output.resolution_kernels_all.classic_res.fwhm = [output.resolution_kernels_all.classic_res.fwhm; output.resolution_kernels.(['curve' num2str(n)]).classic_res.fwhm];
                    output.resolution_kernels_all.classic_res.shape = output.resolution_kernels.(['curve' num2str(n)]).classic_res.shape;
                    output.resolution_kernels_all.cm = output.resolution_kernels.(['curve' num2str(n)]).cm;
                    output.resolution_kernels_all.fwhmwidth = output.resolution_kernels.(['curve' num2str(n)]).fwhmwidth;
                    output.resolution_kernels_all.finesse = output.resolution_kernels.(['curve' num2str(n)]).finesse;
                    
%                     output.resolution_kernels_all.x = [output.resolution_kernels_all.x output.resolution_kernels.(['curve' num2str(n)]).x(temp)];
%                     output.resolution_kernels_all.weight = [output.resolution_kernels_all.weight output.resolution_kernels.(['curve' num2str(n)]).weight(temp)];
%                     output.resolution_kernels_all.fwhm = [output.resolution_kernels_all.fwhm; output.resolution_kernels.(['curve' num2str(n)]).fwhm(temp)];
%                     output.resolution_kernels_all.lambda.fwhm = [output.resolution_kernels_all.lambda.fwhm; output.resolution_kernels.(['curve' num2str(n)]).lambda.fwhm(temp)];
%                     output.resolution_kernels_all.lambda.shape = output.resolution_kernels.(['curve' num2str(n)]).lambda.shape;
%                     output.resolution_kernels_all.theta.fwhm = [output.resolution_kernels_all.theta.fwhm; output.resolution_kernels.(['curve' num2str(n)]).theta.fwhm(temp)];
%                     output.resolution_kernels_all.theta.shape = output.resolution_kernels.(['curve' num2str(n)]).theta.shape;
%                     output.resolution_kernels_all.pixel.fwhm = [output.resolution_kernels_all.pixel.fwhm; output.resolution_kernels.(['curve' num2str(n)]).pixel.fwhm(temp)];
%                     output.resolution_kernels_all.pixel.shape = output.resolution_kernels.(['curve' num2str(n)]).pixel.shape;
%                     output.resolution_kernels_all.binning.fwhm = [output.resolution_kernels_all.binning.fwhm; output.resolution_kernels.(['curve' num2str(n)]).binning.fwhm(temp)];
%                     output.resolution_kernels_all.binning.shape = output.resolution_kernels.(['curve' num2str(n)]).binning.shape;
%                     output.resolution_kernels_all.aperture.fwhm = [output.resolution_kernels_all.aperture.fwhm; output.resolution_kernels.(['curve' num2str(n)]).binning.fwhm(temp)];
%                     output.resolution_kernels_all.aperture.shape = output.resolution_kernels.(['curve' num2str(n)]).aperture.shape;
%                     output.resolution_kernels_all.classic_res.fwhm = [output.resolution_kernels_all.classic_res.fwhm; output.resolution_kernels.(['curve' num2str(n)]).classic_res.fwhm(temp)];
%                     output.resolution_kernels_all.classic_res.shape = output.resolution_kernels.(['curve' num2str(n)]).classic_res.shape;
%                     output.resolution_kernels_all.cm = output.resolution_kernels.(['curve' num2str(n)]).cm;
%                     output.resolution_kernels_all.fwhmwidth = output.resolution_kernels.(['curve' num2str(n)]).fwhmwidth;
%                     output.resolution_kernels_all.finesse = output.resolution_kernels.(['curve' num2str(n)]).finesse;
                end
            else
                output.resolution_kernels_all = [];
            end
        end
        
        output.n_curves = n_curves;
        output.curve_number = curve;
        
        
    case 'fit_it'
        if isempty(curve_handles); return; end
        
        message_handle = grasp_message(['Working: 1D Curve Fitting'],1,'sub');
        
        fitdata = grasp_plot_fit_callbacks_2('get_fit_data'); %get fit data and curve number
        
        %Fit it!. - Using M_Fit
        fun_name = 'pseudo_fn';
        
        start_params = status_flags.fitter.function_info_1d.values;
        %vary_params = start_params;
        %vary_params(find(vary_params==0)) = 1; %i.e. Otherwise if the actual starting param= 0, then vary param would = 0 and parameter would not move in the fit.
        %Only vary parameters that are not fixed
        vary_mask = double(not(status_flags.fitter.function_info_1d.fix));
        %vary_params = vary_params.*vary_mask;
        
        %NOTE:  Not quite sure what is going on with the mf_lsqr routine
        %in terms of what value of the vary parameters should be.  mf_lsqr
        %does not converge and gives Nan's for the errors if this parameter
        %is not correct.  The thing that seems to work best is if it is set
        %to 0.1.  This is then multiplied again inside mf_lsqr.  Until I
        %really understand the mf_lsqr program I can't really work out what
        %is going on.  Maybe look for another least squares minimising
        %routine
        
        
        vary_params = 0.1*vary_mask;
        vary_params = vary_params/100; %randomly chosen to get the curve fit to work properly
        
        %Check all parameters are not fixed
        if sum(vary_params)~=0;
            
            %Run the fit once
            [fit_params, fit_params_err,chi2] = mf_lsqr_grasp(fitdata.xdat_all,fitdata.ydat_all,fitdata.edat_all,start_params,vary_params,fun_name,fitdata);
            
            %transfer chi2 to permanent var
            mf_fitter.temp = chi2;
            %In old Grasp I used to do a second run though the minimiser
            %holding all but one parameter constant at a time to get a
            %better co-varience checked value for the error.  After
            %problems in getting the fitter to work at all this MIGHT BE
            %commented out.  In any case, it seems that mf_lsqr does do a
            %covariance check anyway so might never have been necessary in
            %the first place.
            
            %***** Covariance checking *****
            disp(' ');
            disp('Covariance Checking');
            fit_params_err_cov = zeros(size(fit_params_err)); %only store the error
            params_to_vary = find(vary_params~=0);
            for n = 1:length(params_to_vary);
                temp_vary_params = zeros(size(vary_params));
                temp_vary_params(params_to_vary(n)) = vary_params(params_to_vary(n));
                %Recal the fitting function with the temporary vary_params
                [temp_fit_params, temp_fit_params_err] = mf_lsqr_grasp(fitdata.xdat_all,fitdata.ydat_all,fitdata.edat_all,fit_params,temp_vary_params,fun_name,fitdata);
                fit_params_err_cov(params_to_vary(n)) = temp_fit_params_err(params_to_vary(n));
            end
            
            %Store fit parameters and error in the status_flags
            status_flags.fitter.function_info_1d.values = fit_params';
            status_flags.fitter.function_info_1d.err_values = fit_params_err_cov';
            
            %make final copy of any grouped parameters into theircopy positions
            param_number = 1;
            for fn_multiplex = 1:status_flags.fitter.number1d;
                for variable_loop = 1:status_flags.fitter.function_info_1d.no_parameters
                    if status_flags.fitter.function_info_1d.group(param_number) == 1;
                        status_flags.fitter.function_info_1d.values(param_number) = status_flags.fitter.function_info_1d.values(variable_loop);
                        status_flags.fitter.function_info_1d.err_values(param_number) = status_flags.fitter.function_info_1d.err_values(variable_loop);
                    end
                    param_number = param_number+1;
                end
            end
            
            %Display Covariance corrected fit Parameters and Names
            disp(' ')
            disp('Covariance Corrected Fit Params');
            l = length(status_flags.fitter.function_info_1d.long_names);
            for n = 1:l(1)
                disp([status_flags.fitter.function_info_1d.long_names{n} ' : ' num2str(status_flags.fitter.function_info_1d.values(n)) '  err: ' num2str(status_flags.fitter.function_info_1d.err_values(n))]);
            end
            disp(['Chi^2 = ' num2str(chi2,'%8.3f')]);
            disp(' ')
            
            %Dump the fit paramters, errors and parameter names into a structure that can be picked up elsewhere
            %fit_parameters = struct('name',fun_name,'pnames',pnames,'params',fit_params,'error',fit_params_err_cov,'error_non_cov_check',fit_params_err);
            fit_parameters = status_flags.fitter.function_info_1d;
            fit_parameters.number_functions = status_flags.fitter.number1d;
            
            %Add Fit History to the displayimage structure
            fit_history = {['Curve Fit:' status_flags.fitter.function_info_1d.name]};
            l = length(status_flags.fitter.function_info_1d.variable_names);
            for n = 1:l(1)
                fit_history(length(fit_history)+1,:) = {[status_flags.fitter.function_info_1d.variable_names{n} ' : ' num2str(fit_parameters.values(n)) ' err: ' num2str(fit_parameters.err_values(n))]};
            end
            fit_history(length(fit_history)+1,:) = {['Chi^2 = ' num2str(chi2,'%8.3f')]};
            
            
            
            %Add resolution history
            if not(isempty(fitdata.resolution_kernels_all));
            for n = 1:length(fitdata.resolution_kernels_all.history);
                    fit_history(length(fit_history)+1,:) = fitdata.resolution_kernels_all.history(n);
            end
            end
            
            %Place the history and fit parameters in each curves plot_info
            for n = 1:length(fitdata.curve_number)
                plot_info = get(curve_handles{fitdata.curve_number(n)}(1),'userdata');
                plot_info.fit1d_history = fit_history;
                plot_info.fit1d_parameters = fit_parameters;
                set(curve_handles{fitdata.curve_number(n)}(1),'userdata',plot_info);
            end
        end
        
        %Delete message
        delete(message_handle);
        
        %Draw the fitted function
        grasp_plot_fit_callbacks_2('draw_fn',fit_color); %update the function
        
        %Update the curve fit window with the new parameters
        grasp_plot_fit_callbacks_2('update_curve_fit_window');
end


