function [] = mf_table()
% v 9.2 3/01/2018 E R Louden    

% This function creates a window showing fit parameters, chi^2, and allows 
% the user to `refit' bad fits and export paramters to a .txt file

%% Initialize

% global data structures
global grasp_env;
global mf_fitter;


%% Declutter
% Remove open table window

if(isfield(mf_fitter.handles,'table'))
    if(ishandle(mf_fitter.handles.table))
        close(mf_fitter.handles.table)
    end
end


%% Figure Window

% create figure window.
f = figure('Name', 'Please check fit parameters', 'Position', [ -1260, -100, 1016, 500], 'Color', grasp_env.background_color);


% EXPORT BUTTON
uicontrol(f,'units','normalized','Position',[0.88,0.91,0.1,0.08],...
    'Style','pushbutton','String','Export to .txt','HorizontalAlignment',.... 
    'center','Visible','on','CallBack','mf_export');

% Refit files
uicontrol(f,'units','normalized','Position',[0.02,0.88,0.1,0.085],...
    'Style','text','String','Examine Fit:','tooltip', 'Enter the index of a fit you would like to examine',...
    'HorizontalAlignment','center','Visible','on', 'BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
 
mf_fitter.handles.refit = uicontrol(f,'units','normalized','Position',[0.12,0.925,0.1,0.05],...
                          'Style','edit','tooltip', 'Enter the index of a fit you would like to examine',... 
                          'HorizontalAlignment', 'left','Visible','on', 'Callback', 'mf_fitter_NEWcallbacks(''refit'',  mf_fitter.algorithm_options.num_peaks, str2num( get(mf_fitter.handles.refit, ''String'')))');

%% Update table

uicontrol(f,'units','normalized','Position',[0.32,0.91,0.1,0.08],...
    'Style','pushbutton','String','Update Table','tooltip', 'Update table after custom fit',...
    'HorizontalAlignment','center','Visible','on', 'Callback', 'mf_fitter_NEWcallbacks(''update'',  mf_fitter.algorithm_options.num_peaks, str2num( get(mf_fitter.handles.refit, ''String'')))');
    
% Import data to more clearly named files
cycName = mf_fitter.algorithm_options.current_cycle; 

Numors = mf_fitter.numors;

y0 = mf_fitter.fit_data.(cycName).background;
FWHM = mf_fitter.fit_data.(cycName).fwhm;
I1 = mf_fitter.fit_data.(cycName).intensity1;
X1 = mf_fitter.fit_data.(cycName).center1;
I2 = mf_fitter.fit_data.(cycName).intensity2;
X2 = mf_fitter.fit_data.(cycName).center2;
I3 = mf_fitter.fit_data.(cycName).intensity3;
X3 = mf_fitter.fit_data.(cycName).center3;

% BUG - should include Itot and I1 ?

% Check if the control parameter and Chi2 are avaiable for inclusion
if(isempty(mf_fitter.user_inputs.control_parameter))
    conParam = zeros(mf_fitter.depth,1); 
else
    conParam = str2num(mf_fitter.user_inputs.control_parameter);
end

if(isempty(mf_fitter.fit_data.(cycName).chi2))
    Chi2 = zeros(mf_fitter.depth,1); 
else
    Chi2 = mf_fitter.fit_data.(cycName).chi2;
end


% Check everything for size
corSize = size(Numors);

% find which dimension is the actual length
nS = find(corSize ~= 1); % this will reflect the number of data files
nO = find(corSize == 1); % this will be either 1 or 2 depending on variable

variables = {'conParam','Chi2','y0','FWHM','I1','X1','I2','X2','I3','X3'};
for v = 1:length(variables)
    varSize = size(eval(variables{v}));
    
    if(varSize(nS) == corSize(nS))
        % everything is good to go
    elseif(varSize(nO) == corSize(nS))
        % the variable simply needs to be flipped
        eval([variables{v} '= (transpose(eval(variables{v})))'])
    else
        % there is a dimensional mismatch
        disp(['error: variable ' variables{v} ' does not have the correct dimensions'])
    end
    
end


%% Nonphysical fits
% Check for negative intensities - these are nonphysical fits
negative_numors = '';

for i = 1:length(mf_fitter.numors)
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


%% Build the table

% Titles at top of table
Titles = {'Numors', 'Control Parameter,','Chi-Squared', 'y0', 'y0 Error', 'FWHM', 'FWHM Error', 'I1', 'I1 Error','X1','X1 Error','I2', 'I2 Error', 'X2', 'X2 Error', 'I3','I3 Error', 'X3','X3 Error'};

% Creating checkboxes at end
ColumnType = {[], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []};

% Which columns can be edited (only checkboxes)
EditableNum= [false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false];

% Concatenating all data
DataTable = [Numors, conParam, Chi2, y0, FWHM, I1, X1, I2, X2, I3, X3];

DataCell = num2cell(DataTable);

% BUG - find a better way to check if other quantities (such as fms, delta
% phi) have been calculated to insert into table
%         if(isempty(mf_fitter.fit_data.(cycName).fms)== 0)
%             FMS = mf_fitter.fit_data.(cycName).fms';
%             % BUG - do I still propagate the errors at some point?
%             %FMS2 = mf_fitter.error_prop.sig_frac;
%             FGS = mf_fitter.fit_data.(cycName).fgs';
%             mf_fitter.handles.ExportTable = table(Numors, conParam, Chi2, y0, FWHM, I1, X1, I2, X2, I3, X3, FMS, FMS2, FGS);
%         else
            mf_fitter.handles.ExportTable = table(Numors, conParam, Chi2, y0, FWHM, I1, X1, I2, X2, I3, X3);
%         end

        
%% Draw the Table

mf_fitter.handles.test_table = uitable('Parent', f, 'Data', DataCell, 'ColumnName', Titles,...
                         'ColumnFormat', ColumnType, 'ColumnEditable', EditableNum,...
                         'Units', 'normalized', 'Position', [0,0,1,.9], 'ColumnWidth', {70});
                     
mf_fitter.handles.table = gcf;
                     
end