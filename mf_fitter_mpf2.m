function [ ] = mf_fitter_mpf2( )
    
% User provides two, "good" images, provides guess values for initial fits

% This fitting algorithm operates under the assumption that the peak
% postions can potentially move.  This movement is not strictly enforced,
% thus in the case of the perpendicular AC Cycles, where the peaks are
% likely fixed, it is still able to perfrom an accurate fit.

% First cycle: y0, c2, fwhm fixed; guess values for i0, c1, c3 scale with the 'index' of the file
% Second cycle: prompts user to refit or fix peaks not confident in
% Third cycle: 4 step process, user supplies when peaks are fixed
%               groups y0, fhwm; never fixes y0 or i0
%               1. fixed c, fwhm to averages for c2 & fwhm, to gui supplied for c1 & c3
%               2. fixed c, determine fwhm
%               3. fixed fwhm to adv_avg fwhm, determines c2
%               4. fixed c2 to adv_avg c2, determines i0
% Fourth cycle: use avg c2 & fwhm, fix to the fit data for i0, y0; determine final GS c
% Final cycle: fit one last time with grasp, same as 3.1, determine final i0

% v 8.0
% 5/13/2015 MFF Liz

%number of fitting cycles
Nc = 6;
wait_position = [0.85 .2 .25 .08];
refit_check = 1;

%% INITIALIZE global vars
global mf_fitter;
global status_flags;
global grasp_handles;

mf_fitter_callbacks('initialize');


%% FIT reference peaks
notification = msgbox('Fitting reference peaks...');

mf_fitter.good_peaks.inner_num = [];
mf_fitter.good_peaks.outers_num = [];

%locates good files' index 
mf_fitter.good_peaks.inner_num = find(mf_fitter.fit_data.names == mf_fitter.good_peaks.inner);
mf_fitter.good_peaks.outers_num = find(mf_fitter.fit_data.names == mf_fitter.good_peaks.outers);

%Checks to see if both successfully found
if isempty(mf_fitter.good_peaks.inner_num)
    disp('One-peak file could not be found');
    return
end
if isempty(mf_fitter.good_peaks.outers_num)
    disp('Two-peak file could not be found');
    return
end

%Fitting first (one peak) reference file
mf_fitter.flag = 1;
mf_fitter_callbacks('set_free',1,mf_fitter.good_peaks.inner_num); 
mf_fitter_callbacks('fit',1,mf_fitter.good_peaks.inner_num);
pause(0.25);
close(mf_fitter.handles.plot_handle);
mf_fitter.flag = 0;

if ( get(mf_fitter.handles.three_peak,'Value') )
  mf_fitter_callbacks('set_free',3,mf_fitter.good_peaks.inner_num);
  mf_fitter_callbacks('fit',3,mf_fitter.good_peaks.outers_num);
  pause(0.25);
  close(mf_fitter.handles.plot_handle);
else
    %Fitting second (two peak) reference file
    mf_fitter_callbacks('set_free',2,mf_fitter.good_peaks.inner_num);
    mf_fitter_callbacks('fit',2,mf_fitter.good_peaks.outers_num);
    pause(0.25);
    close(mf_fitter.handles.plot_handle);
end

close(notification);


%% Check for EIA
%player = audioplayer(mf_fitter.awesome.y, mf_fitter.awesome.Fs);
%if get(mf_fitter.handles.awesome, 'Value')   
%    play(player)
%end


%% AVERAGE widths, background, & MS center
mf_fitter_callbacks('avg',0,0,'fwhm');
mf_fitter_callbacks('avg',0,0,'background');
mf_fitter_callbacks('avg',0,0,'center2');


%%  FITTING CYCLE 1
% FIT ALL w/ 3 peak with matlab fitter, fixed background, fwhm, & x2 to average values
% guess values for intensity and cycles scale with how far along the transition is

%Progress bar
name = ['Fitting Cycle: 1 of ', num2str(Nc), '...']
notification = waitbar(0, name, 'Name', 'Please wait','units', 'Normalized', 'Position', wait_position);

%Fitting cycle
for i = 1:mf_fitter.depth
    
    try
        %Perform Matlab fit
        mf_fitter_callbacks('matlab_fit1',3,i)

        %Array to check if all files fit properly
        mf_fitter.did_fit(i) = 1;
        mf_fitter.fit_data.chi2(i,1) = mf_fitter.temp;
     
    catch
         mf_fitter.did_fit(i) = 0;
         disp('fit error');
    
    end
    
    pause(0.25); 
    close(mf_fitter.handles.plot_handle);
    
    %Update wait bar each cycle
    waitbar(i/(Nc*mf_fitter.depth));

end

% Display fit params in table (optional)
mf_fitter_table;
close(notification);

mf_fitter_callbacks('center_store');
mf_fitter_callbacks('center_plot')

mf_fitter_callbacks('fraction')
pause(2)
close(gcf)


%%  FITTING CYCLE 2
% Prompts user to re-examine bad peaks
% Select the last peak you were confident in, attempt to refit once, then set what value to fix it to.

name = ['Fitting Cycle: 2 of ', num2str(Nc), '...']
notification = waitbar(0, name, 'Name', 'Please wait','units', 'Normalized', 'Position', wait_position);

mf_fitter.refit.continue = 0;
refit_GUI_window('initialize');
while(mf_fitter.refit.continue ~= 1)
   pause(2)
end

close(notification);
close(mf_fitter.handles.table);
close(mf_fitter.handles.cplot);

mf_fitter_table;

mf_fitter_callbacks('center_store');
mf_fitter_callbacks('center_plot')


%%  FITTING CYCLE 3
% FIT ALL w/ 3 peak with grasp fitter, fixed fwhm & x2

name = ['Fitting Cycle: 3 of ', num2str(Nc), '...']
notification = waitbar(0, name, 'Name', 'Please wait','units', 'Normalized', 'Position', wait_position);

% Uses the values from the refit gui window
num_refit = length(mf_fitter.fix.center1);
for i=(num_refit+1):mf_fitter.depth
    mf_fitter.fix.center1(i) = mf_fitter.fit_data.center1(i);
    mf_fitter.fix.center3(i) = mf_fitter.fit_data.center3(i);
end
    
for j=1:4
    for i = 1:mf_fitter.depth   
        try
            mf_fitter_callbacks('set_mfp1',3,i,j);
            mf_fitter_callbacks('fit',3,i);

            mf_fitter.did_fit(i) = 1;
            mf_fitter.fit_data.chi2(i,1) = mf_fitter.temp;

         catch
             mf_fitter.did_fit(i) = 0;
             disp('fit error');

        end

        pause(0.25);
        close(mf_fitter.handles.plot_handle);
        
        waitbar((2/Nc)+i/(j*Nc*mf_fitter.depth));
    end
    
    if(j==2)
        figure
        errorbar(mf_fitter.fit_data.names',mf_fitter.fit_data.fwhm(:,1),mf_fitter.fit_data.fwhm(:,2));
        xlabel('File Name');
        ylabel('FWHM');
        mf_fitter_callbacks('adv_avg',0,0,'fwhm');
    end
    
    if(j==3)
       mf_fitter_callbacks('adv_avg',0,0,'center2');
    end

    close(mf_fitter.handles.table);
    close(mf_fitter.handles.cplot);

    mf_fitter_table;
    mf_fitter_callbacks('center_store');
    mf_fitter_callbacks('center_plot')

end

close(notification)
close(mf_fitter.handles.table)
close(mf_fitter.handles.cplot)

mf_fitter_table;


%% AVERAGE widths, background, & MS center
mf_fitter_callbacks('avg',0,0,'fwhm');
mf_fitter_callbacks('avg',0,0,'background');
mf_fitter_callbacks('adv_avg',0,0,'center2');


 %%  FITTING CYCLE 4
% FIT ALL w/ 3 peak with matlab fitter, fixed fwhm & x2

name = ['Fitting Cycle: 4 of ', num2str(Nc), '...']
notification = waitbar(0, name, 'Name', 'Please wait','units', 'Normalized', 'Position', wait_position);


for i = 1:mf_fitter.depth
    
    try
        mf_fitter_callbacks('matlab_fit4',3,i)
    
        mf_fitter.did_fit(i) = 1;
        mf_fitter.fit_data.chi2(i,1) = mf_fitter.temp;
        
    catch
        mf_fitter.did_fit(i) = 0;
        disp('fit error');
    
    end
    
    pause(0.25);  
    close(mf_fitter.handles.plot_handle);

    waitbar((4/5)+i/(Nc*mf_fitter.depth));

end

close(notification);
close(mf_fitter.handles.table)

mf_fitter_table;
mf_fitter_callbacks('center_store');
mf_fitter_callbacks('center_plot');

% save final center plot
mf_fitter.handles.fig = gcf;
name = [mf_fitter.exp ' center_plot']
mf_fitter.file = name;
mf_fitter_callbacks('save',0,0,2)
close(gcf)

mf_fitter_callbacks('adv_avg',0,0,'center2');

%%  Fitting Cycle 5

name = ['Fitting Cycle: 5 of ', num2str(Nc), '...'];
notification = waitbar(0, name, 'Name', 'Please wait','units', 'Normalized', 'Position', wait_position);

A = isnan(mf_fitter.fit_data.center1(:,2));
B = isnan(mf_fitter.fit_data.center3(:,2));
       
goodA = [];
goodA = find(A==0);
bothgood = 0;
i=1;
j=1;

% Determine where both center values have a non-infinite error bar
for i=1:length(goodA)
        if(B(goodA(i))==0) 
            bothgood(j) = goodA(i);
            j = j+1;
        end
        i = i+1;
end
    
num_refit = bothgood(refit_check);

mf_fitter.fix.center1(1:num_refit) = mean(mf_fitter.fit_data.center1(num_refit:(num_refit+1))) 
mf_fitter.fix.center1 = mf_fitter.fix.center1';
mf_fitter.fix.center3(1:num_refit) = mean(mf_fitter.fit_data.center3(num_refit:(num_refit+1))) 
mf_fitter.fix.center3 = mf_fitter.fix.center3';

for i=num_refit:mf_fitter.depth
    mf_fitter.fix.center1(i) = mf_fitter.fit_data.center1(i)
    mf_fitter.fix.center3(i) = mf_fitter.fit_data.center3(i)
end

mf_fitter_callbacks('adv_avg',0,0,'center2');

% create the directory to store the individual data files
dir = [mf_fitter.extension mf_fitter.folder '/Individual Fits'];
mkdir(dir);

for i = 1:mf_fitter.depth
    
    try
        mf_fitter_callbacks('set_mfp1',3,i,1); 
        mf_fitter_callbacks('fit',3,i);
        
        close(mf_fitter.handles.plot_handle)
        mf_fitter_callbacks('plot_fit',3,i,0);
        
        mf_fitter.did_fit(i) = 1;
        mf_fitter.fit_data.chi2(i,1) = mf_fitter.temp;     
    
        mf_fitter.handles.fig = gcf;
        plot_template()
        name = ['individual Fits/File_' num2str(mf_fitter.fit_data.names(i))];
        mf_fitter.file = name;
        mf_fitter_callbacks('save',0,0,1);
        pause(0.25);
        close(mf_fitter.handles.plot_handle);

    catch
         mf_fitter.did_fit(i) = 0;
         disp('fit error');
         
    end
    
    waitbar((4/Nc)+i/(Nc*mf_fitter.depth));
end

close(notification);
close(mf_fitter.handles.table)

%% save fraction and colormap plots
save_names = {'fraction', 'colormap'};
for i=1:length(save_names)
    variable = char(save_names(i));
    mf_fitter_callbacks(variable);
    pause(3)
    mf_fitter.handles.fig = gcf;
    name = [mf_fitter.exp ' ' variable];
    mf_fitter.file = name;
    mf_fitter_callbacks('save',0,0,2);
    close(gcf);
end

mf_fitter_table;
fileName = [mf_fitter.extension mf_fitter.folder '/data.txt'];
writetable(mf_fitter.handles.ExportTable, fileName);

end