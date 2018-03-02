function [ ] = mf_save( to_do, opt1, FigHandle )
% v 9.2 2/28/2018 E R Louden    

% This function opens the save options GUI window in order to set the folder 
% extension, control parameters, & other save preferences.
% It also exports the figures and data tables according to those preferences.
% The callbacks include:
%   change_save_option
%   save
%   save_files
%   window



%% Initialize

% global data structures
global grasp_env
global grasp_handles
global status_flags
global mf_fitter

% error code in case no figure handle is required
if (nargin<2)
    opt1 = 'void';
    FigHandle = 'null';
end
if (nargin<3)
    FigHandle = 'null';
end


%% Callbacks - window & save

switch to_do
        
    
    case 'change_save_options'
        if(strcmp(opt1,'none'))
           disp('please specifiy which save_option needs to be changed using the opt1 variable') 
        else
            % check if passing more than just one field
            if(strfind(opt1,'.')) 
                opt1a = opt1(1:(strfind(opt1,'.')-1));
                opt1b = opt1((strfind(opt1,'.')+1):length(opt1));
                temp = get(mf_fitter.handles.(opt1a).(opt1b),'string');
                if not(isempty(temp)) 
                    mf_fitter.(char(opt1a)).(char(opt1b)) = temp;
                end
            % only one field    
            else
                temp = get(mf_fitter.handles.(char(opt1)),'string');
                if not(isempty(temp)) 
                    mf_fitter.(char(opt1)) = temp;
                end
            end
        end

        
 
    case 'save'
        % Save Script
        
        dir = [mf_fitter.save_options.extension mf_fitter.save_options.folder '_MFFv9.2/'];
        mkdir(dir);
               
        if(mf_fitter.save_options.file_type.fig)
            figfile = [dir opt1 '.fig'];
            saveas(FigHandle,figfile);
        end
        if(mf_fitter.save_options.file_type.jpg)
            jpgfile = [dir opt1 '.jpg'];
            saveas(FigHandle,jpgfile);
        end
        if(mf_fitter.save_options.file_type.eps)
            epsfile = [dir opt1 '.eps'];
            saveas(FigHandle,epsfile);
        end
        if(mf_fitter.save_options.file_type.pdf)
            pdffile = [dir opt1 '.pdf'];
            saveas(FigHandle,pdffile);
        end
        if(mf_fitter.save_options.file_type.fig + mf_fitter.save_options.file_type.jpg + mf_fitter.save_options.file_type.eps + mf_fitter.save_options.file_type.pdf == 0)
            disp('no file types selected - check save window')
        end    
        
        
        
    case 'save_files'
        % checks desired save file types from check boxes on GUI
        
        mf_fitter.save_options.file_type.fig = get(mf_fitter.handles.save_options.figFile,'Value');
        mf_fitter.save_options.file_type.jpg = get(mf_fitter.handles.save_options.jpgFile,'Value');
        mf_fitter.save_options.file_type.eps = get(mf_fitter.handles.save_options.epsFile,'Value');
        mf_fitter.save_options.file_type.pdf = get(mf_fitter.handles.save_options.pdfFile,'Value');
        
    
        
    case 'window'
        
        % Draw Window
        fig_position = [800, 590, 500, 250];

        grasp_handles.window_modules.save.window = figure(....
            'units','pixels',....
            'Position',fig_position,....
            'Name','Experiment & Save Options' ,....
            'NumberTitle', 'off',....
            'Tag','mf_save_window',...
            'color',grasp_env.background_color,....
            'menubar','none',....
            'resize','off'); 

        handle = grasp_handles.window_modules.save.window;
        mf_fitter.handles.save_GUI = handle;

        % Some initial Values
        if(isempty(mf_fitter.save_options.extension))
            extString = 'enter extension string here';
        else
            extString = mf_fitter.save_options.extension;
        end

        if(isempty(mf_fitter.user_inputs.control_parameter))
            conParamString = 'enter comma separated values here';
        else
            conParamString = num2str(mf_fitter.user_inputs.control_parameter);
        end
        
        if(isempty(mf_fitter.save_options.folder))
            folderString = 'enter folder string here';
        else
            folderString = mf_fitter.save_options.folder;
        end

        
        % GUI Inputs 
        
        %Experiment Extension 
        uicontrol(handle,'units','normalized','Position',[0.1 0.8 0.5 0.12],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','center','Style','text','String','Experiment Extension','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
        mf_fitter.handles.save_options.extension = uicontrol(handle,'units','normalized','Position',[0.1 0.75 0.8 0.1],'tooltip','Sets the extension for where to save files','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','edit','String', extString,'HorizontalAlignment','left','Tag','angular_bin','Visible','on','Callback','mf_save(''change_save_options'',''save_options.extension'')'); %'mf_GUI_callbacks(''change_ext'')');

        %Experiment Folder
        uicontrol(handle,'units','normalized','Position',[0.1 0.6 0.5 0.12],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','center','Style','text','String','Data Run Folder ','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
        mf_fitter.handles.save_options.folder = uicontrol(handle,'units','normalized','Position',[0.1 0.55 0.8 0.1],'tooltip','Particular data run folder','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','edit','String',folderString,'HorizontalAlignment','left','Tag','angular_bin','Visible','on','Callback','mf_save(''change_save_options'',''save_options.folder'')'); %'mf_GUI_callbacks(''change_folder'')');
        
        %Experiment Cycles
        uicontrol(handle,'units','normalized','Position',[0.1 0.4 0.5 0.12],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','center','Style','text','String','Control Parameter: (i.e. Applied AC Cyles)','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
        mf_fitter.handles.user_inputs.control_parameter = uicontrol(handle,'units','normalized','Position',[0.1 0.35 0.8 0.1],'tooltip','Stores the control parameter values for plots','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','edit','String',conParamString,'HorizontalAlignment','left','Tag','angular_bin','Visible','on','Callback','mf_save(''change_save_options'',''user_inputs.control_parameter'')'); %'mf_GUI_callbacks(''change_cycs'')');

        %File Type Check Boxes
        %fig
        val = mf_fitter.save_options.file_type.fig;
        uicontrol(handle,'units','normalized','Position',[0.2 0.2 0.3 0.05],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','left','Style','text','String','.fig?','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
        mf_fitter.handles.save_options.figFile = uicontrol(handle,'units','normalized','Position',[0.2 0.1 0.04 0.08],'tooltip','Check to save as a matlab figure','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','checkbox','HorizontalAlignment','right','Tag','angular_bin','Visible','on','Value', val,'Callback','mf_save(''save_files'')'); %'mf_GUI_callbacks(''save_files'')');
        
        %jpeg
        val = mf_fitter.save_options.file_type.jpg;
        uicontrol(handle,'units','normalized','Position',[0.4 0.2 0.3 0.05],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','left','Style','text','String','.jpg?','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
        mf_fitter.handles.save_options.jpgFile = uicontrol(handle,'units','normalized','Position',[0.4 0.1 0.04 0.08],'tooltip','Check to save as a jpeg','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','checkbox','HorizontalAlignment','right','Tag','angular_bin','Visible','on','Value', val,'Callback','mf_save(''save_files'')');%'mf_GUI_callbacks(''save_files'')');
        
        %eps
        val = mf_fitter.save_options.file_type.eps;
        uicontrol(handle,'units','normalized','Position',[0.6 0.2 0.3 0.05],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','left','Style','text','String','.eps?','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
        mf_fitter.handles.save_options.epsFile = uicontrol(handle,'units','normalized','Position',[0.6 0.1 0.04 0.08],'tooltip','Check to save as an eps','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','checkbox','HorizontalAlignment','right','Tag','angular_bin','Visible','on','Value', val,'Callback','mf_save(''save_files'')');%'mf_GUI_callbacks(''save_files'')');
        
        %pdf
        val = mf_fitter.save_options.file_type.pdf;
        uicontrol(handle,'units','normalized','Position',[0.8 0.2 0.3 0.05],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','left','Style','text','String','.pdf?','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
        mf_fitter.handles.save_options.pdfFile = uicontrol(handle,'units','normalized','Position',[0.8 0.1 0.04 0.08],'tooltip','Check to save as a pdf','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','checkbox','HorizontalAlignment','right','Tag','angular_bin','Visible','on','Value', val,'Callback','mf_save(''save_files'')'); %'mf_GUI_callbacks(''save_files'')');
        
        % BUG - Button to save all open files


        % BUG - export case - this does not need to be its own function
end

