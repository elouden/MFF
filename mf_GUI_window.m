function mf_GUI_window
% Draws main GUI for MFF, pretends to be Sector Window so grasp plays nice

% v 8.0
% 11/11/16 MFF Liz


global grasp_env
global grasp_handles
global status_flags
global mf_fitter

if ishandle(grasp_handles.window_modules.sector.window); delete(grasp_handles.window_modules.sector.window); end
grasp_handles.window_modules.sector.inner_radius

%% Initialize
mf_GUI_callbacks('initialize');

%% Draw Window
% monitor:  
%fig_position = [-1250, 580, 310, 310];
% laptop:
fig_position = [1000, 580, 350, 310];

grasp_handles.window_modules.sector.window = figure(....
    'units','pixels',....
    'Position',fig_position,....
    'Name','Multi-file Fitter' ,....
    'NumberTitle', 'off',....
    'Tag','mf_GUI_window',...
    'color',grasp_env.background_color,....
    'menubar','none',....
    'resize','off',....
    'closerequestfcn','sector_callbacks2(''close'');closereq');

handle = grasp_handles.window_modules.sector.window;

panel = uipanel(handle, 'Title', 'Sector', 'Units', 'normalized','FontSize',10,'BackgroundColor',grasp_env.background_color,'Position',[.05 .27 .5 .7], 'ForegroundColor', [1 1 1]);
panel2 = uipanel(handle, 'Title', 'Smoothing', 'Units', 'normalized','FontSize',10,'BackgroundColor',grasp_env.background_color,'Position',[.05 .05 .5 .2], 'ForegroundColor', [1 1 1]);

%% Sector Box & Binning
%Inner radius
uicontrol(panel,'units','normalized','Position',[0.05,0.85,0.5,0.06],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','center','Style','text','String','Inner Radius:','BackgroundColor',grasp_env.background_color, 'ForegroundColor', [1 1 1]);
grasp_handles.window_modules.sector.inner_radius = uicontrol(panel,'units','normalized','Position',[0.6,0.85,0.35,0.06],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','edit','String',num2str(status_flags.analysis_modules.sectors.inner_radius),'HorizontalAlignment','left','Visible','on','CallBack','sector_callbacks2(''inner_radius'');');

%Outer Radius
uicontrol(panel,'units','normalized','Position',[0.05,0.73,0.51,0.06],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','center','Style','text','String','Outer Radius:','BackgroundColor',grasp_env.background_color, 'ForegroundColor', [1 1 1]);
grasp_handles.window_modules.sector.outer_radius = uicontrol(panel,'units','normalized','Position',[0.6,0.73,0.35,0.06],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','edit','String',num2str(status_flags.analysis_modules.sectors.outer_radius),'HorizontalAlignment','left','Visible','on','CallBack','sector_callbacks2(''outer_radius'');');

%Theta
uicontrol(panel,'units','normalized','Position',[0.05,0.54,0.5,0.06],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','center','Style','text','String','Theta:','BackgroundColor',grasp_env.background_color, 'ForegroundColor', [1 1 1]);
grasp_handles.window_modules.sector.theta = uicontrol(panel,'units','normalized','Position',[0.6,0.54,0.35,0.06],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','edit','String',num2str(status_flags.analysis_modules.sectors.theta),'HorizontalAlignment','left','Visible','on','CallBack','sector_callbacks2(''theta'');');

%Delta Theta
uicontrol(panel,'units','normalized','Position',[0.05,0.42,0.5,0.06],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','center','Style','text','String','Delta Theta:','BackgroundColor',grasp_env.background_color, 'ForegroundColor', [1 1 1]);
grasp_handles.window_modules.sector.delta_theta = uicontrol(panel,'units','normalized','Position',[0.6,0.42,0.35,0.06],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','edit','String',num2str(status_flags.analysis_modules.sectors.delta_theta),'HorizontalAlignment','left','Tag','sector_deltatheta','Visible','on','CallBack','sector_callbacks2(''delta_theta'');');

%Angle Binning
uicontrol(panel,'units','normalized','Position',[0.05 0.23 0.5 0.06],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','center','Style','text','String','Angle Bin:','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
grasp_handles.window_modules.radial_average.azimuth_bin = uicontrol(panel,'units','normalized','Position',[0.6 0.23 0.35 0.06],'tooltip','Re-Bin Size in Angle (degrees).  Note: Min Value Depends on the Sector Minimum Radius','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','edit','HorizontalAlignment','left','String', '1', 'Tag','angular_bin','Visible','on','callback','radial_average_callbacks2(''azimuth_bin'');');

% Default Values
%set(grasp_handles.window_modules.sector.inner_radius,'String',38);
%set(grasp_handles.window_modules.sector.outer_radius,'String',51);
%set(grasp_handles.window_modules.sector.theta,'String',304.5);
%set(grasp_handles.window_modules.sector.delta_theta,'String',50);
%set(grasp_handles.window_modules.radial_average.azimuth_bin,'String',1.2);


%APPLY SECTOR
uicontrol(panel,'units','normalized','Position',[0.3, 0.07, 0.4, 0.08],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','pushbutton','String','Apply','HorizontalAlignment','center','Tag','boxfox_button','Visible','on','CallBack','mf_GUI_callbacks(''apply_sector'')');

%% Smoothing

% smoothing fwhm
uicontrol(panel2,'units','normalized','Position',[0.05 0.65 0.25 0.25],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','left','Style','text','String','FWHM:','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
mf_fitter.handles.smoothing.fwhm = uicontrol(panel2,'units','normalized','Position',[0.05 0.25 0.25 0.3],'tooltip','FWHM for smoothing algorithm','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','edit','String',num2str(mf_fitter.smoothing.fwhm),'HorizontalAlignment','left','Tag','angular_bin','Visible','on','Callback', 'mf_GUI_callbacks(''smoothing_fwhm'')');

% smoothing step size
uicontrol(panel2,'units','normalized','Position',[0.39 0.65 0.35 0.25],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','left','Style','text','String','Step Size:','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
mf_fitter.handles.smoothing.step = uicontrol(panel2,'units','normalized','Position',[0.4 0.25 0.25 0.3],'tooltip','Phi step size for plotting','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','edit','String',num2str(mf_fitter.smoothing.step),'HorizontalAlignment','left','Tag','angular_bin','Visible','on','Callback', 'mf_GUI_callbacks(''smoothing_stepsize'')');

%smoothing on?
uicontrol(panel2,'units','normalized','Position',[0.7 0.65 0.3 0.25],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','left','Style','text','String','Smooth?','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
mf_fitter.handles.smoothing.switch = uicontrol(panel2,'units','normalized','Position',[0.77 0.27 0.1 0.3],'tooltip','Check if you would like to perform data smoothing','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','checkbox','HorizontalAlignment','right','Tag','angular_bin','Visible','on','Value', 1);


%% Other Options

% FITTING PARADIGMS
%uicontrol(handle,'units','normalized','Position',[0.62 0.8 0.35, 0.05],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','left','Style','text','String','peaks move?','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
mf_fitter.handles.fit_dropdown = uicontrol('units','normalized','Position',[0.6 0.9, 0.35 0.075],'tooltip','Select a fitting routine','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','popup','string',{'Fitting Paradigms','Peaks Fixed','2K Peaks Grow','14K Peaks Rotate','Raster Scan', 'Angle Binning Scan', 'I vs Q'}, 'Callback', 'mf_GUI_callbacks(''fit_dropdown'')');
%mf_fitter.handles.dropdown = uicontrol('units','normalized','Position',[0.60 0.03 0.3 0.06],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','popup','string',{'Options','Fraction Plot','Centers Plot','Intensity Colormap','Peak Decay','Angle Binning Scan','Raster Scan'}, 'Callback', 'mf_GUI_callbacks(''dropdown'')');

% PLOT
%FRACTION, COLORMAP DROPDOWN MENU
%mf_fitter.handles.dropdown = uicontrol('units','normalized','Position',[0.60 0.03 0.35 0.06],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','popup','string',{'Options','Fraction Plot','Centers Plot','Intensity Colormap','Peak Decay','Angle Binning Scan','Raster Scan'}, 'Callback', 'mf_GUI_callbacks(''dropdown'')');
mf_fitter.handles.plot_dropdown = uicontrol('units','normalized','Position',[0.60 0.8 0.35 0.075],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','popup','string',{'Plot Options', 'FMS','Peak Centers', 'Peak Separation', 'Intensity Colormap', 'Peak Decay'}, 'Callback', 'mf_GUI_callbacks(''plot_dropdown'')');



%1 and 2 pk files
uicontrol(handle,'units','normalized','Position',[0.6 0.66 0.17 0.12],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','left','Style','text','String','1-Peak:','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
mf_fitter.handles.good_inner = uicontrol(handle,'units','normalized','Position',[0.6 0.66 0.12 0.06],'tooltip','Select a file with ONE peak','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','edit','String',num2str(mf_fitter.good_peaks.inner),'HorizontalAlignment','left','Tag','angular_bin','Visible','on','Callback', 'mf_GUI_callbacks(''good_peaks'')');

uicontrol(handle,'units','normalized','Position',[0.75 0.66 0.17 0.12],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','left','Style','text','String','2-Peak:','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
mf_fitter.handles.good_outers = uicontrol(handle,'units','normalized','Position',[0.75 0.66 0.12 0.06],'tooltip','Select a file with TWO peaks','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','edit','String',num2str(mf_fitter.good_peaks.outers),'HorizontalAlignment','left','Tag','angular_bin','Visible','on','Callback', 'mf_GUI_callbacks(''good_peaks'')');

%3-peak
uicontrol(handle,'units','normalized','Position',[0.88 0.73 0.3 0.05],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','left','Style','text','String','3-Peak?','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
mf_fitter.handles.three_peak = uicontrol(handle,'units','normalized','Position',[0.9 0.67 0.048 0.05],'tooltip','Check if there is no good 1-peak reference','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','checkbox','HorizontalAlignment','right','Tag','angular_bin','Visible','on','Value', 0);



%Peak Percentage Cutoff
uicontrol(handle,'units','normalized','Position',[0.6 0.5 0.15 0.12],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','center','Style','text','String','Intensity Cutoff:','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
mf_fitter.handles.int_cutoff = uicontrol(handle,'units','normalized','Position',[0.6 0.47 0.15 0.06],'tooltip','Sets minimum fraction of total intensity to be considered a peak','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','edit','String',num2str(mf_fitter.int_cutoff),'HorizontalAlignment','left','Tag','angular_bin','Visible','on','Callback', 'mf_GUI_callbacks(''change_int_cutoff'')');

uicontrol(handle,'units','normalized','Position',[0.8 0.5 0.15 0.12],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','center','Style','text','String','Peak Factor:','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
mf_fitter.handles.dphi_cutoff = uicontrol(handle,'units','normalized','Position',[0.8 0.47 0.15 0.06],'tooltip','Multiplicative factor - how much larger than FWHM to be considered for avereages','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','edit','String',num2str(mf_fitter.dphi_cutoff),'HorizontalAlignment','left','Tag','angular_bin','Visible','on','Callback', 'mf_GUI_callbacks(''change_dphi_cutoff'')');



%FIT TABLE BUTTON
uicontrol('units','normalized','Position',[0.625 0.32 0.3 0.1],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','pushbutton','string','View Params', 'Callback', 'mf_fitter_table');

% Experiment & Save Options
uicontrol('units','normalized','Position',[0.625 0.2 0.2 0.1],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','pushbutton','tooltip','Click to enter the save folder information and applied cycles','string','Save Options','Callback','mf_fitter_save(''window'')');
uicontrol(handle,'units','normalized','Position',[0.84 0.26 0.35 0.05],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','left','Style','text','String','Save?','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
mf_fitter.handles.save = uicontrol(handle,'units','normalized','Position',[0.85 0.215 0.048 0.05],'tooltip','Check if you want to auto save files based on the save options','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','checkbox','HorizontalAlignment','right','Tag','angular_bin','Visible','on','Value', 0);
    

%LOOP BUTTON
uicontrol('units','normalized','Position',[0.675 0.03 0.2 0.15],'FontName',grasp_env.font,'FontSize',1.25*grasp_env.fontsize,'Style','pushbutton','string','GO!','Callback','mf_GUI_callbacks(''go'')');
   



%% Old options 

%FWHM fix box
%uicontrol(handle,'units','normalized','Position',[0.62 0.45 0.35 0.05],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','left','Style','text','String','Fix FWHM?','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
%mf_fitter.handles.fix_fwhm = uicontrol(handle,'units','normalized','Position',[0.87 0.45 0.048 0.05],'tooltip','Would you like the FWHM to be fixed in the final fit?','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','checkbox','HorizontalAlignment','right','Tag','angular_bin','Visible','on','Value', 1);

%awesome?
%uicontrol(handle,'units','normalized','Position',[0.66 0.35 0.35 0.05],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','left','Style','text','String','EIA?','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
%mf_fitter.handles.awesome = uicontrol(handle,'units','normalized','Position',[0.87 0.35 0.048 0.05],'tooltip','???','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','checkbox','HorizontalAlignment','right','Tag','angular_bin','Visible','on','Value', 0);


 












%Garbage which is hidden, but kept to play nice with grasp_update
%Anisotropy & Anisotropy angle
enable = 'off';
uicontrol(handle,'units','normalized','Position',[0.07,0.26,0.25,0.06],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','right','Style','text','String','Anisotropy:','BackgroundColor',grasp_env.background_color, 'ForegroundColor', [1 1 1], 'Visible', 'off');
grasp_handles.window_modules.sector.anisotropy = uicontrol(handle,'units','normalized','Position',[0.35,0.26,0.16,0.06],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','edit','String',num2str(status_flags.analysis_modules.sectors.anisotropy),'HorizontalAlignment','left','Visible','off','CallBack','sector_callbacks(''anisotropy'');','enable',enable);
uicontrol(handle,'units','normalized','Position',[0.53,0.26,0.15,0.06],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','right','Style','text','String','Angle:','BackgroundColor',grasp_env.background_color, 'ForegroundColor', [1 1 1], 'Visible', 'off');
grasp_handles.window_modules.sector.anisotropy_angle = uicontrol(handle,'units','normalized','Position',[0.7,0.26,0.16,0.06],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','edit','String',num2str(status_flags.analysis_modules.sectors.anisotropy_angle),'HorizontalAlignment','left','Visible','off','CallBack','sector_callbacks(''anisotropy_angle'');','enable',enable);


%Quick Angles +30 +45 +60 +90
uicontrol(handle,'units','normalized','Position',[0.55,0.58,0.1,0.06],'string','+30','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','pushbutton','Tag','sector_theta_plus30','Visible','off','Value',0,...
    'CallBack','sector_callbacks(''angle_plus'',30);');
uicontrol(handle,'units','normalized','Position',[0.67,0.58,0.1,0.06],'string','+45','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','pushbutton','Tag','sector_theta_plus45','Visible','off','Value',0,...
    'CallBack','sector_callbacks(''angle_plus'',45);');
uicontrol(handle,'units','normalized','Position',[0.55,0.50,0.1,0.06],'string','+60','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','pushbutton','Tag','sector_theta_plus60','Visible','off','Value',0,...
    'CallBack','sector_callbacks(''angle_plus'',60);');
uicontrol(handle,'units','normalized','Position',[0.67,0.50,0.1,0.06],'string','+90','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','pushbutton','Tag','sector_theta_plus90','Visible','off','Value',0,...
    'CallBack','sector_callbacks(''angle_plus'',90);');

%Mirror Sectors
number_of_sectors_string = '1';
for n = 2:status_flags.analysis_modules.sectors.mirror_sectors_max
    number_of_sectors_string = [number_of_sectors_string '|' num2str(n)];
end
uicontrol(handle,'units','normalized','Position',[0.65, 0.9,0.3,0.06],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','right','Style','text','String','Mirror Sectors:','BackgroundColor',grasp_env.background_color, 'ForegroundColor', [1 1 1], 'Visible', 'off');
grasp_handles.window_modules.sector.mirror_sectors = uicontrol(handle,'units','normalized','Position',[0.7,0.85,0.25,0.06],'HorizontalAlignment','center','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','Popup','Tag','number_of_sectors',...
    'String',number_of_sectors_string,'CallBack','sector_callbacks(''mirror_sectors'');','Value',status_flags.analysis_modules.sectors.mirror_sectors, 'Visible', 'off');


%Refresh & draw the sectors
sector_callbacks2;
