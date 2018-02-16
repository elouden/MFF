function [ ] = mf_fitter_save( to_do, name, FigHandle )
%Opens GUI to set folder extension & other save options

% v 9.0
% 11/12/2016 MFF Liz

global grasp_env
global grasp_handles
global status_flags
global mf_fitter

if (nargin<3)
    name = 'void';
    FigHandle = 'null';
end

switch to_do
    
    case 'window'
        %% Draw Window
        % monitor:  
        %fig_position = [-1250, 580, 310, 310];
        % laptop:
        fig_position = [800, 590, 500, 250];

        grasp_handles.window_modules.save.window = figure(....
            'units','pixels',....
            'Position',fig_position,....
            'Name','Experiment & Save Options' ,....
            'NumberTitle', 'off',....
            'Tag','mf_save_window',...
            'color',grasp_env.background_color,....
            'menubar','none',....
            'resize','off'); %,....
            %'closerequestfcn','sector_callbacks2(''close'');closereq');

        handle = grasp_handles.window_modules.save.window;

        %panel = uipanel(handle, 'Title', 'Sector','Units', 'normalized','FontSize',10,'BackgroundColor',grasp_env.background_color,'Position',[.05 .05 .5 .9], 'ForegroundColor', [1 1 1]);

        %% Some initial Values
        if(isempty(mf_fitter.extension))
            extString = 'enter extension string here';
        else
            extString = mf_fitter.extension;
        end

        if(isempty(mf_fitter.fit_data.cycles))
            cycString = 'enter comma separated values here';
        else
            cycString = num2str(mf_fitter.fit_data.cycles);
        end
        
        if(isempty(mf_fitter.fit_data.cycles))
            folderString = 'enter folder string here';
        else
            folderString = mf_fitter.folder;
        end

        %% GUI Inputs 

        %Experiment Extension 
        uicontrol(handle,'units','normalized','Position',[0.1 0.8 0.5 0.12],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','center','Style','text','String','Experiment Extension','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
        mf_fitter.handles.extension = uicontrol(handle,'units','normalized','Position',[0.1 0.75 0.8 0.1],'tooltip','Sets the extension for where to save files','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','edit','String', extString,'HorizontalAlignment','left','Tag','angular_bin','Visible','on','Callback', 'mf_GUI_callbacks(''change_ext'')');

        %Experiment Folder
        uicontrol(handle,'units','normalized','Position',[0.1 0.6 0.5 0.12],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','center','Style','text','String','Data Run Folder ','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
        mf_fitter.handles.folder = uicontrol(handle,'units','normalized','Position',[0.1 0.55 0.8 0.1],'tooltip','Particular data run folder','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','edit','String',folderString,'HorizontalAlignment','left','Tag','angular_bin','Visible','on','Callback', 'mf_GUI_callbacks(''change_folder'')');
        
        %Experiment Cycles
        uicontrol(handle,'units','normalized','Position',[0.1 0.4 0.5 0.12],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','center','Style','text','String','Control Parameter: (i.e. Applied AC Cyles)','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
        mf_fitter.handles.cycles = uicontrol(handle,'units','normalized','Position',[0.1 0.35 0.8 0.1],'tooltip','Stores the control parameter values for plots','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','edit','String',cycString,'HorizontalAlignment','left','Tag','angular_bin','Visible','on','Callback', 'mf_GUI_callbacks(''change_cycs'')');

        %File Type Check Boxes
        %fig
        val = mf_fitter.save.fig;
        uicontrol(handle,'units','normalized','Position',[0.2 0.2 0.3 0.05],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','left','Style','text','String','.fig?','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
        mf_fitter.handles.figFile = uicontrol(handle,'units','normalized','Position',[0.2 0.1 0.04 0.08],'tooltip','Check to save as a matlab figure','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','checkbox','HorizontalAlignment','right','Tag','angular_bin','Visible','on','Value', val,'Callback','mf_GUI_callbacks(''save_files'')');
        
        %jpeg
        val = mf_fitter.save.jpg;
        uicontrol(handle,'units','normalized','Position',[0.4 0.2 0.3 0.05],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','left','Style','text','String','.jpg?','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
        mf_fitter.handles.jpgFile = uicontrol(handle,'units','normalized','Position',[0.4 0.1 0.04 0.08],'tooltip','Check to save as a jpeg','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','checkbox','HorizontalAlignment','right','Tag','angular_bin','Visible','on','Value', val,'Callback','mf_GUI_callbacks(''save_files'')');
        
        %eps
        val = mf_fitter.save.eps;
        uicontrol(handle,'units','normalized','Position',[0.6 0.2 0.3 0.05],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','left','Style','text','String','.eps?','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
        mf_fitter.handles.epsFile = uicontrol(handle,'units','normalized','Position',[0.6 0.1 0.04 0.08],'tooltip','Check to save as an eps','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','checkbox','HorizontalAlignment','right','Tag','angular_bin','Visible','on','Value', val,'Callback','mf_GUI_callbacks(''save_files'')');
        
        %pdf
        val = mf_fitter.save.pdf;
        uicontrol(handle,'units','normalized','Position',[0.8 0.2 0.3 0.05],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','left','Style','text','String','.pdf?','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
        mf_fitter.handles.pdfFile = uicontrol(handle,'units','normalized','Position',[0.8 0.1 0.04 0.08],'tooltip','Check to save as a pdf','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','checkbox','HorizontalAlignment','right','Tag','angular_bin','Visible','on','Value', val,'Callback','mf_GUI_callbacks(''save_files'')');
        
% Button to save all open files

    case 'save'
%% Save Script
        
        %name = [mf_fitter.file ' AB ' num2str(current_ab)];
        dir = [mf_fitter.extension mf_fitter.folder '_Buzz/'];
        mkdir(dir);
               
        if(mf_fitter.save.fig)
            figfile = [dir name '.fig'];
            saveas(FigHandle,figfile);
        end
        if(mf_fitter.save.jpg)
            jpgfile = [dir name '.jpg'];
            saveas(FigHandle,jpgfile);
        end
        if(mf_fitter.save.eps)
            epsfile = [dir name '.eps'];
            saveas(FigHandle,epsfile);
        end
        if(mf_fitter.save.pdf)
            pdffile = [dir name '.pdf'];
            saveas(FigHandle,pdffile);
        end
        if(mf_fitter.save.fig + mf_fitter.save.jpg + mf_fitter.save.eps + mf_fitter.save.pdf == 0)
            disp('no file types selected - check save window')
        end        
    
end

