function [] = mf_fitter_table()

%Creates window showing fit parameters, chi^2, and allows user to export to txt. 

% v 9.1
% 5/8/2017 MFF Liz

global grasp_env;
global mf_fitter;

% stop the pile-up of fit parameter tables!
if(isfield(mf_fitter.handles,'table'))
    if(ishandle(mf_fitter.handles.table))
        close(mf_fitter.handles.table)
    end
end

%create figure window.
f = figure('Name', 'Please check fit parameters', 'Position', [ -1260, -100, 1016, 500], 'Color', grasp_env.background_color);


%EXPORT BUTTON
uicontrol(f,'units','normalized','Position',[0.88,0.91,0.1,0.08],...
    'Style','pushbutton','String','Export to .txt','HorizontalAlignment',.... 
    'center','Visible','on','CallBack','mf_export');

%Refit file
uicontrol(f,'units','normalized','Position',[0.02,0.88,0.1,0.085],...
    'Style','text','String','Examine Fit:','tooltip', 'Enter the index of a fit you would like to examine',...
    'HorizontalAlignment','center','Visible','on', 'BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
 
mf_fitter.handles.refit = uicontrol(f,'units','normalized','Position',[0.12,0.925,0.1,0.05],...
                          'Style','edit','tooltip', 'Enter the index of a fit you would like to examine',... 
                          'HorizontalAlignment', 'left','Visible','on', 'Callback', 'mf_fitter_NEWcallbacks(''refit'',  3, str2num( get(mf_fitter.handles.refit, ''String'')))');

%Update table
uicontrol(f,'units','normalized','Position',[0.32,0.91,0.1,0.08],...
    'Style','pushbutton','String','Update Table','tooltip', 'Update table after custom fit',...
    'HorizontalAlignment','center','Visible','on', 'Callback', 'mf_fitter_NEWcallbacks(''update'',  3 , str2num( get(mf_fitter.handles.refit, ''String'')))');
    
%Import data to more clearly named files
Numors = mf_fitter.fit_data.names;
Background = mf_fitter.fit_data.background;
FWHM = mf_fitter.fit_data.fwhm(:,1);
I1 = mf_fitter.fit_data.intensity1;
X1 = mf_fitter.fit_data.center1;
I2 = mf_fitter.fit_data.intensity2;
X2 = mf_fitter.fit_data.center2;
I3 = mf_fitter.fit_data.intensity3;
X3 = mf_fitter.fit_data.center3;
Chi2 = mf_fitter.fit_data.chi2;


% Check for negative intensities
negative_numors = '';

for i = 1:length(mf_fitter.fit_data.names)
    if (I1(i,1) < 0 || I2(i,1) < 0 || I3(i,1) < 0)
        negative_numors = [negative_numors num2str(Numors(i,1)) ', '];       
    end
end
        
% Negative intensities 
uicontrol(f,'units','normalized','Position',[0.48,0.88,0.1,0.085],...
    'Style','text','String','Negative Intensities:',...
    'HorizontalAlignment','center','Visible','on', 'BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
      
uicontrol(f,'units','normalized','Position',[0.6,0.91,0.2,0.08],...
    'Style','text','String',negative_numors, 'ForegroundColor', [1 0 0],...
    'HorizontalAlignment','center', 'Visible','on');

        if(isempty(mf_fitter.fit_data.cycles))
            Cycle = zeros(mf_fitter.depth,1); 
        else
            cycSize = size(mf_fitter.fit_data.cycles);
            if(cycSize(1) == 1)
                Cycle = mf_fitter.fit_data.cycles';
            else
                Cycle = mf_fitter.fit_data.cycles;
            end
        end
        if(isempty(Chi2))
            s = size(I1);
            Chi2 = zeros(s);
            Chi2 = Chi2(:,1);
        end

% Titles at top of table
Titles = {'Numors', 'Cycles,','Chi-Squared', 'FWHM', 'I1', 'I1 Error','X1','X1 Error','I2', 'I2 Error', 'X2', 'X2 Error', 'I3','I3 Error', 'X3','X3 Error'};

% Creating checkboxes at end
ColumnType = {[], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []};

% Which columns can be edited (only checkboxes)
EditableNum= [false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false];

% Concatenating all data
DataTable = [Numors, Cycle, Chi2, FWHM, I1, X1, I2, X2, I3, X3];

DataCell = num2cell(DataTable);

        if(isempty(mf_fitter.fit_data.fms)== 0)
            FMS = mf_fitter.fit_data.fms';
            FMS2 = mf_fitter.error_prop.sig_frac;
            FGS = mf_fitter.fit_data.fgs';
            mf_fitter.handles.ExportTable = table(Numors, Cycle, Chi2, Background, FWHM, I1, X1, I2, X2, I3, X3, FMS, FMS2, FGS);
        else
            %Creating table
            mf_fitter.handles.ExportTable = table(Numors, Cycle, Chi2, Background, FWHM, I1, X1, I2, X2, I3, X3);
        end

%drawing table
mf_fitter.handles.test_table = uitable('Parent', f, 'Data', DataCell, 'ColumnName', Titles,...
                         'ColumnFormat', ColumnType, 'ColumnEditable', EditableNum,...
                         'Units', 'normalized', 'Position', [0,0,1,.9], 'ColumnWidth', {70});
                     
mf_fitter.handles.table = gcf;
                     
end