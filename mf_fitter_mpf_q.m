function [ ] = mf_fitter_mpf_q( )   
%   Fits the I vs Q      

% v 9.2
% 1/8/2017 MFF Liz

%%
T = 14;
offset = 272;
%9.3G, 14K A peakpos1 = [-8.45187517700000;-7.68915335500000;-7.47352154800000;-6.87633342300000;-6.69900075100000;-6.37595843600000;-5.66311063200000;-5.26966137300000;-4.77753474300000;-4.57352948500000;-3.93995958600000;-3.24127178800000;-2.63441652700000;-2.08772010300000;-2.16719148300000;-2.29879820400000;-2.33026628200000];
%9.3G, 14K A peakpos2 = [5.90762311800000;4.92568230300000;4.39117319500000;4.30777653900000;4.41996344500000;3.93916621700000;3.25245582100000;2.68130154000000;2.53264448400000;2.11011047600000;2.34221284000000;2.19088092000000;1.65485484800000;1.29943095100000;1.18069352600000;1.20828193600000;1.29487810600000];
%9.3G, 14K B peakpos1 = [-6.26000000000000;-5.50000000000000;-3.85000000000000;-3.25000000000000;-1.52000000000000;-1.12000000000000]
peakpos = [6.91000000000000;6.42000000000000;5.55000000000000;5.19000000000000;4.04000000000000;2.10000000000000];

%% INITIALIZE global vars and function fitting
global mf_fitter;
global status_flags;
global grasp_handles;

% determine number of files ('depth') for use elsewhere
mf_fitter.depth = status_flags.selector.fdpth_max - 1;


% set the curve fit general parameters that are always the same
status_flags.fitter.function_info_1d.name = 'Gaussian';
status_flags.fitter.function_info_1d.variable_names = {'y0','i0', 'xc','fwhm'};
status_flags.fitter.function_info_1d.long_names = {'Background','Integrated Intensity','Centre','FWHM'}
status_flags.fitter.function_info_1d.math_code = {'y = y0 + (i0 / (fwhm*sqrt(pi/2) / sqrt(log(4)))) .* exp((-2*(x-xc).^2) / (fwhm.^2 / log(4)));'};
status_flags.fitter.function_info_1d.auto_guess_code = {'[amp i_amp] = max(y);    %Peak Intensity','xc = x(i_amp);         %Centre Position','y0 = min(y);         %Background','x_diff = diff(x); %Differences between adjacent x''s.','ytemp = y(1:length(y)-1);','i0 = abs((sum((ytemp-y0).*x_diff)));','sigma = i0/(amp*sqrt(2*pi));','fwhm = 2*sqrt(2*log(2))*sigma; %The sqrt is because this is a 2D gaussian','guess_values = [y0, i0, xc, fwhm];'};
status_flags.fitter.function_info_1d.point_click_code = {'text_handle = grasp_message([''Click on Background''],1,''sub''); %Grasp Message','[x y0]=ginput(1); %Mouse input','delete(text_handle);pause(0.1);','text_handle = grasp_message([''Click on Peak''],1,''sub'');  %Grasp Message','[xc amp]=ginput(1); %Mouse input, Centre and Peak intensity','delete(text_handle);pause(0.1);','text_handle = grasp_message([''Click on Peak Width''],1,''sub'');  %Grasp Message','[fwhm y]=ginput(1); %Mouse input','delete(text_handle);pause(0.1);','fwhm=abs(fwhm-xc);','amp=amp-y0;','i0 = amp*fwhm*sqrt(pi)/(sqrt(2)*sqrt(2*log(2))); %Convert Peak Intensity to IntInt','guess_values = [y0, i0, xc, fwhm];'}

% initialize to zeros
status_flags_names = {'group','fix','values','err_values'}
for i=1:length(status_flags_names)
    variable = char(status_flags_names(i));
   status_flags.fitter.function_info_1d.(variable) = zeros(1,4); 
end      
       
status_flags.fitter.function_info_1d.no_parameters = 4;
status_flags.fitter.fn1d = 1;


mf_fitter.exp = mf_fitter.folder;


%% Perform the I vs Q fit

for img_num = 1:mf_fitter.depth  

    
        % for 14 K change sector to make fitted peaks
        if(T == 14)
            val = peakpos(img_num)+offset
            set(grasp_handles.window_modules.sector.theta,'String',num2str(val))
            sector_callbacks2('theta');
        end
        
        status_flags.fitter.number1d = 1;
        
        % load desired depth file and update main grasp GUI
        % the +1 is necessary as file 1 represents a sum
        status_flags.selector.fd = img_num+1;
        grasp_update;

        % plot I vs q and create the current figure, store handle
        radial_average_callbacks('averaging_control','radial_q');
        mf_fitter.handles.plot_handle = gcf;
        grasp_update;
        
        % Perform auto guess (only for 1-peak fit) and fit functions
        grasp_plot_fit_callbacks('auto_guess','auto_guess');
        
        % having trouble auto-guessing power, start it at 5, fix it temporarily so auto-guess doesn't overwrite
        %status_flags.fitter.function_info_1d.fix(3) = 1;
        %status_flags.fitter.function_info_1d.values(3) = 5;
        %grasp_update;
       
        
        % turn off "fix" for power
        %status_flags.fitter.function_info_1d.fix(3) = 0;
        %grasp_update;
        
        % perform the fit
        set(grasp_handles.window_modules.curve_fit1d.curve_number,'value',1)
        grasp_plot_fit_callbacks('fit_it');
        
        % store the data - loop through values & error
         for j=1:2
             
            if(j==1) h = status_flags.fitter.function_info_1d.values; 
            else h = status_flags.fitter.function_info_1d.err_values;
            end
       
%             mf_fitter.IvsQ_fit.background(img_num,j) = h(1);
%             mf_fitter.IvsQ_fit.multiplier(img_num,j) = h(2);
%             mf_fitter.IvsQ_fit.power(img_num,j) = h(3);
%             mf_fitter.IvsQ_fit.intensity(img_num,j) = h(4);
%             mf_fitter.IvsQ_fit.center(img_num,j) = h(5);
%             mf_fitter.IvsQ_fit.fwhm(img_num,j) = h(6);

            mf_fitter.IvsQ_fit.background(img_num,j) = h(1);
            mf_fitter.IvsQ_fit.intensity(img_num,j) = h(2);
            mf_fitter.IvsQ_fit.center(img_num,j) = h(3);
            mf_fitter.IvsQ_fit.fwhm(img_num,j) = h(4);

        end

        
end  

%% Make a table to export the fit data

% Import data to more clearly named files
Numors = mf_fitter.fit_data.names;
Background = mf_fitter.IvsQ_fit.background;
%Mult = mf_fitter.IvsQ_fit.background;
%Power = mf_fitter.IvsQ_fit.power;
I = mf_fitter.IvsQ_fit.intensity;
Q = mf_fitter.IvsQ_fit.center;
FWHM = mf_fitter.IvsQ_fit.fwhm;

% Get Cycles
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

% Titles at top of table
%Titles = {'Numors', 'Cycles','Q','Q_err','y0','y0_err','m','m_err','n','n_err','I','I_err','FWHM','FWHM_err'};
 Titles = {'Numors', 'Cycles','Q','Q_err','y0','y0_err','I','I_err','FWHM','FWHM_err'};
   
% Concatenating all data
DataTable = [Numors, Cycle, Q, Background, I, FWHM];
DataCell = num2cell(DataTable);

mf_fitter.handles.Export_IvsQ_Table = table(Numors, Cycle, Background, Q, I, FWHM);
    
    % case 'data_storage'
        % set the handle based on whether grabbing the values or the errors
       
% mf_fitter_table;
% fileName = [mf_fitter.extension mf_fitter.folder '_Buzz/' 'data.txt'];
% writetable(mf_fitter.handles.ExportTable, fileName);
    

