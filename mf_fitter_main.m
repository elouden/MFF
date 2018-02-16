function [ ] = mf_fitter_main( )

%% MFF Multi-File Fitting Algorithm

% User provides two, "good" images, provides centers for initial fits

% First cycle: centers fixed, 3-peak fit fixed, FWHM allowed to vary.
% Second cycle: FWHM fixed to avg, centers allowed to vary, num of peaks
%               determined by fractional intensity cutoff.
% Final cycle: centers fixed to avg, FWHM fixed, 3-peak fixed.

% v 8.0
% 5/13/2015 MFF Liz

global mf_fitter;
global status_flags;


% Experiment Details
%mf_fitter.extension = '/Users/edewaard/Documents/Eskildsen Research/MgB2 Buzz Analysis/May 2015 Figures & Fits/';
%mf_fitter.extension = '/Users/edewaard/Dropbox/Eskildsen Research/2014_11 MgB2 ORNL/'
%mf_fitter.extension = '/Users/edewaard/Dropbox/Eskildsen Research/2016_09 MgB2 ILL D33/pure Analysis/'
%mf_fitter.folder = 'MFF_Buzz w_peakrot_strexp_RCRO/'
%mf_fitter.folder = 'MFF_Buzz/'
%mf_fitter.exp = 'ORNL CG2 11_2014 RCRO'
mf_fitter.exp = 'ILL D33 9_2016'

% NEW FITTING PARADIGM
if ( get(mf_fitter.handles.move_peaks,'Value') )
    disp('Fitting routine based on the MS peaks rotating out from the central peak has been selected')
    %mf_fitter_mpf2()
    % new routine, make drop down menu later
    mf_fitter_mpf3()
else

%% INITIALIZE global vars

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

player = audioplayer(mf_fitter.awesome.y, mf_fitter.awesome.Fs);

if get(mf_fitter.handles.awesome, 'Value')
   
    play(player)

end

%% AVERAGE widths & background

mf_fitter_callbacks('avg',0,0,'fwhm');
mf_fitter_callbacks('avg',0,0,'background');

%% FIT ALL w/ 3 peak, fixed centers, free FWHM

%Progress bar
wait_position = [-.525 .115 .25 .08];
notification = waitbar(0,'Fitting cycle: 1 of 3...', 'Name', 'Please wait','units', 'Normalized', 'Position', wait_position);

%Fitting cycle
for i = 1:mf_fitter.depth
    
    try
        %Setting parameters for fit
        mf_fitter_callbacks('set_centers',3,i); 
        
        %Fitting
        mf_fitter_callbacks('fit',3,i);
        
        %Array to check if all files fit properly
        mf_fitter.did_fit(i) = 1;
        mf_fitter.fit_data.chi2(i,1) = mf_fitter.temp;
    
    catch
        mf_fitter.did_fit(i) = 0;
        disp('fuck you mf_lsqr');
    end
    
    pause(0.25);
    close(mf_fitter.handles.plot_handle);
    
    %Update wait bar each cycle
    waitbar(i/(3*mf_fitter.depth));

end

%Display fit params in table (optional)
mf_fitter_table;
gcf

%% Check for failed fits and display
failed_fits = {};

for i = 1:mf_fitter.depth
    
    if not(mf_fitter.did_fit(i))
           disp('Fit failed: ')
           disp(mf_fitter.fit_data.names(i,1))
 %       failed_fits = failed_fits + num2str(mf_fitter.fit_data.names(i,1)) + ' ';
    failed_fits = 1;
    end
end

if not(isempty(failed_fits))
    
%    notification = msgbox('The following files did not fit:' + failed_fits, 'Icon', 'warn');
    mf_fitter_table;
    
end
%     

%% AVERAGE Background & FWHM
mf_fitter_callbacks('avg',0,0,'background')
mf_fitter_callbacks('adv_avg',0,0,'fwhm')
%a = median(mf_fitter.fit_data.fwhm);
%mf_fitter.averages.fwhm = a(1,1);

%% DISCRIMINATE to find peaks
mf_fitter_callbacks(0,'discriminate');


%% FIT ALL w/ mixed peaks, loose centers

%Update progress bar
waitbar(0.333,notification,'Fitting cycle: 2 of 3...');

for i = 1:mf_fitter.depth
    
    %Calculate num of peaks from discriminate
    num_peaks = mf_fitter.inner_peak(i,1) + 2*mf_fitter.outer_peaks(i,1);
      
    try
        mf_fitter_callbacks('set_fwhm', num_peaks,i); 
        mf_fitter_callbacks('fit',3,i);
        mf_fitter.did_fit(i) = 1;
    
    catch
        mf_fitter.did_fit(i) = 0;
        disp('fuck you mf_lsqr');
    end
    
    pause(0.25);
    
    close(mf_fitter.handles.plot_handle);
    waitbar((i+mf_fitter.depth)/(3*mf_fitter.depth));

end

%Display params (optional)
mf_fitter_table;
mf_fitter_callbacks('center_store')
mf_fitter_callbacks('center_plot')


%% AVERAGE centers
mf_fitter_callbacks('adv_avg',0,0,'center1');
mf_fitter_callbacks('adv_avg',0,0,'center2');
mf_fitter_callbacks('adv_avg',0,0,'center3');

   
%% FIT ALL w/ three peaks, fixed centers

waitbar(0.667,notification,'Fitting cycle: 3 of 3...');

if get(mf_fitter.handles.fix_fwhm, 'Value')
    
    %FIT WITH FIXED CENTERS AND FIXED FWHM

     for i = 1:mf_fitter.depth
    
        try
            mf_fitter_callbacks('set_both',3,i); 
            mf_fitter_callbacks('fit',3,i);
            mf_fitter.did_fit(i) = 1;
    
        catch
            mf_fitter.did_fit(i) = 0;
            disp('fuck you mf_lsqr');
        end
    
        pause(0.25);
        close(mf_fitter.handles.plot_handle);
        waitbar((i+2*mf_fitter.depth)/(3*mf_fitter.depth));

    end
else 
    
    %FIT JUST WITH FIXED CENTERS 
    
    for i = 1:mf_fitter.depth
    
        try
            mf_fitter_callbacks('set_centers',3,i); 
            mf_fitter_callbacks('fit',3,i);
            mf_fitter.did_fit(i) = 1;
    
        catch
            mf_fitter.did_fit(i) = 0;
            disp('fuck you mf_lsqr');
        end
    
        pause(0.25);
        close(mf_fitter.handles.plot_handle);
        waitbar((i+2*mf_fitter.depth)/(3*mf_fitter.depth));

    end
end

close(notification);

%% DISPLAY final table
mf_fitter_table;

%% CALCULATE Fraction in MS and GS
%mf_fitter_callbacks('fraction');

for i = 1:3;

beep

pause(2)

end

stop(player)

end
end

