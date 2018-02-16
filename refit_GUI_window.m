function [ done ] = refit_GUI_window( to_do )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% v 7.4
% 2/4/15 MFF Liz


global grasp_env
global grasp_handles
global status_flags
global mf_fitter

if ishandle(mf_fitter.handles.refit); delete(mf_fitter.handles.refit); end
%grasp_handles.window_modules.sector.inner_radius

%% Initialize
%mf_GUI_callbacks('initialize');
mf_fitter.refit.continue = 0;

%% Callbacks
switch to_do
    
    case 'draw_window'
        fig_position = [-1250, 580, 310, 310];
        mf_fitter.handles.refit = figure(....
            'units','pixels',....
            'Position',fig_position,....
            'Name','Refit Window' ,....
            'NumberTitle', 'off',....
            'Tag','refit_GUI_window',...
            'color',grasp_env.background_color,....
            'menubar','none',....
            'resize','off')
            %'closerequestfcn','closereq');

        handle = mf_fitter.handles.refit;

        % Move on Button
        uicontrol('units','normalized','Position',[0.65 0.18 0.25 0.06],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','pushbutton','string','Move On','callBack','refit_GUI_window(''continue'')');

        % Break Point
        uicontrol(handle,'units','normalized','Position',[0.1 0.85 0.15 0.12],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','center','Style','text','String','Breakpoint:','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
        mf_fitter.handles.refit.bp = uicontrol(handle,'units','normalized','Position',[0.1 0.82 0.15 0.06],'tooltip','Enter break point file indices separated by a comma','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','edit','String',num2str(mf_fitter.refit.bp),'HorizontalAlignment','left','Tag','angular_bin','Visible','on','Callback', 'refit_GUI_window(''change_bps'')');

        % Center1
        uicontrol(handle,'units','normalized','Position',[0.3 0.85 0.15 0.12],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','center','Style','text','String','Center 1 Value:','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
        mf_fitter.handles.refit.x1 = uicontrol(handle,'units','normalized','Position',[0.3 0.82 0.15 0.06],'tooltip','Enter break point file indices separated by a comma','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','edit','String','-','HorizontalAlignment','left','Tag','angular_bin','Visible','on','Callback', 'refit_GUI_window(''change_x1'')');

        % Center1
        uicontrol(handle,'units','normalized','Position',[0.3 0.70 0.15 0.12],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','center','Style','text','String','Center 3 Value:','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
        mf_fitter.handles.refit.x3 = uicontrol(handle,'units','normalized','Position',[0.3 0.67 0.15 0.06],'tooltip','Enter value for center','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','edit','String','-','HorizontalAlignment','left','Tag','angular_bin','Visible','on','Callback', 'refit_GUI_window(''change_x3'')');

        % Fix?
        uicontrol(handle,'units','normalized','Position',[0.2 0.5 0.35 0.05],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','left','Style','text','String','fix?','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
        mf_fitter.handles.fixcenters = uicontrol(handle,'units','normalized','Position',[0.3 0.55 0.048 0.05],'tooltip','If checked, values will be fixed instead of just used as a guess','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','checkbox','HorizontalAlignment','right','Tag','fix_checkbox','Visible','on','Value',1);

        % Try it Button
        uicontrol('units','normalized','Position',[0.65 0.3 0.25 0.06],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','pushbutton','string','Try it!','callBack','refit_GUI_window(''try_it'')');

        
    case 'initialize'
        mf_fitter.refit.continue = 0;
        mf_fitter.refit.bp = 0;
        mf_fitter.fix.center1 = [];
        mf_fitter.fix.center3 = [];
        %start(refit_timer);
        refit_GUI_window('draw_window');
        
    case 'continue'
        mf_fitter.refit.continue = 1;
        %done = 1;
        %return
        %stop(refit_timer);
        
    case 'change_bps'
        mf_fitter.refit.bp = str2num(get(mf_fitter.handles.refit.bp,'String'));
        
    case 'change_x1'
        mf_fitter.refit.x1_val = str2num(get(mf_fitter.handles.refit.x1,'String'));
            
    case 'change_x3'
        mf_fitter.refit.x3_val = str2num(get(mf_fitter.handles.refit.x3,'String'));
        
    case 'try_it'
        bp = mf_fitter.refit.bp
        x1_val = mf_fitter.refit.x1_val
        x3_val = mf_fitter.refit.x3_val
        
        if (length(x1_val)==1) mf_fitter.fit_data.center1(1:bp) = mf_fitter.fit_data.center1(x1_val)
        else mf_fitter.fit_data.center1(1:bp) = mean(mf_fitter.fit_data.center1(x1_val))
        end
        if (length(x3_val)==1) mf_fitter.fit_data.center3(1:bp) = mf_fitter.fit_data.center3(x1_val)
        else mf_fitter.fit_data.center3(1:bp) = mean(mf_fitter.fit_data.center3(x1_val))
        end
        for i =1:bp
             mf_fitter_callbacks('matlab_fit4',3,i)
             close(mf_fitter.handles.plot_handle);
        end
        
        mf_fitter_table;
        mf_fitter_callbacks('center_store');
        mf_fitter_callbacks('center_plot')
        if( get(mf_fitter.handles.fixcenters,'Value') == 1)
            if (length(x1_val)==1) mf_fitter.fix.center1(1:bp) = mf_fitter.fit_data.center1(x1_val)
            else mf_fitter.fix.center1(1:bp) = mean(mf_fitter.fit_data.center1(x1_val))
            end
            %mf_fitter.fix.center1(1:num_refit) = mean(mf_fitter.fit_data.center1(num_refit:(num_refit-1)))
            mf_fitter.fix.center1 = mf_fitter.fix.center1';
            if (length(x3_val)==1) mf_fitter.fit_data.center3(1:bp) = mf_fitter.fit_data.center3(x1_val)
            else mf_fitter.fix.center3(1:bp) = mean(mf_fitter.fit_data.center3(x1_val))
            end
            %mf_fitter.fix.center3(1:num_refit) = mean(mf_fitter.fit_data.center3(num_refit:(num_refit-1))) 
            mf_fitter.fix.center3 = mf_fitter.fix.center3';
        end



end
%Peak Percentage Cutoff
%uicontrol(handle,'units','normalized','Position',[0.6 0.85 0.15 0.12],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','center','Style','text','String','Move on','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
%mf_fitter.handles.cutoff =
%uicontrol(handle,'units','normalized','Position',[0.6 0.82 0.15 0.06],'tooltip','Ready to move on?','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','button','String',num2str(mf_fitter.cutoff),'HorizontalAlignment','left','Tag','angular_bin','Visible','on','Callback', 'mf_GUI_callbacks(''change_cutoff'')');

%uicontrol(handle,'units','normalized','Position',[0.8 0.85 0.15 0.12],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','center','Style','text','String','Colormap Limits:','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
%mf_fitter.handles.climits = uicontrol(handle,'units','normalized','Position',[0.8 0.82 0.15 0.06],'tooltip','Enter TWO values for the colormap limits','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','edit','String',num2str(mf_fitter.climits),'HorizontalAlignment','left','Tag','angular_bin','Visible','on','Callback', 'mf_GUI_callbacks(''change_climits'')');

%1 and 2 pk files
%uicontrol(handle,'units','normalized','Position',[0.6 0.67 0.15 0.12],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','center','Style','text','String','1-Peak File:','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
%mf_fitter.handles.good_inner = uicontrol(handle,'units','normalized','Position',[0.6 0.65 0.15 0.06],'tooltip','Select a file with ONE peak','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','edit','String',num2str(mf_fitter.good_peaks.inner),'HorizontalAlignment','left','Tag','angular_bin','Visible','on','Callback', 'mf_GUI_callbacks(''good_peaks'')');

%uicontrol(handle,'units','normalized','Position',[0.8 0.67 0.15 0.12],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','center','Style','text','String','2-Peak File:','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
%mf_fitter.handles.good_outers = uicontrol(handle,'units','normalized','Position',[0.8 0.65 0.15 0.06],'tooltip','Select a file with TWO peaks','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','edit','String',num2str(mf_fitter.good_peaks.outers),'HorizontalAlignment','left','Tag','angular_bin','Visible','on','Callback', 'mf_GUI_callbacks(''good_peaks'')');

end

