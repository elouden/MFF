function [ ] = mf_fitter_mpf2K( )
% v 9.2 2/19/2018 E R Louden    

% This function performs the actual fitting cycle.  This fitting algorithm
% is ideal for the supercooled (2K) VL transition data.  While all peaks
% are fit, initial guess values and order of fitting  have been tailored
% for when the ES peaks nucleate at their final orientations and grow.
% This discontinuous transitions necessitates 3 Bragg peaks.

% Fitting Cycle Order-
%       Reference Files:    fits a 1-peak and a 2-peak file, used to get initial guess values (GRASP)
%       Second Cycle:       Fix centers & fwhm, determine intesnities (MATLAB)
%       Third Cycle:        Fix intensiites & fwhm, determine centers  (Matlab)
%       Fourth Cycle:       Fix centers, intensities free, determine fwhm from files with I/Itot > 10%  (MATLAB)
%       Fifth Cycle:        Fix fwhm, centers & intensities free (MATLAB)


%% INITIALIZE 

% global data structures
global mf_fitter;
global status_flags;
global grasp_handles;

%mf_fitter_callbacks('initialize');
mf_fitter_NEWcallbacks('initialize');


% number of fitting cycles
mf_fitter.algorithm_options.total_cycles = 6;
mf_fitter.algorithm_options.fitter = {'G','M','M','M','M','M'};
mf_fitter.algorithm_options.fitter_cycle = [];
mf_fitter.algorithm_options.num_peaks = 3;
mf_fitter.algorithm_options.algorithm_name = '2K';


% wait_position = [0.85 .2 .25 .08];
% refit_check = 1;

% get file save name
fileName = [mf_fitter.save_options.extension mf_fitter.save_options.folder '_MFFv9.2/'];



%% Check which Fitting Cycles are to be viewed
view1 = mf_fitter.algorithm_options.cycle_view.view1;
view2 = mf_fitter.algorithm_options.cycle_view.view2;
view3 = mf_fitter.algorithm_options.cycle_view.view3;
view4 = mf_fitter.algorithm_options.cycle_view.view4;
view5 = mf_fitter.algorithm_options.cycle_view.view5;


%% Reference Files - Cycle 1
% For simplicity, fit reference peaks in GRASP 

mf_fitter.algorithm_options.current_cycle = 'cycle1';
notification = msgbox('Fitting reference peaks...');
mf_fitter.algorithm_options.fitter_cycle = 1;

%Checks to see if both successfully found
if isempty(mf_fitter.user_inputs.reference_files.one)
    disp('One-peak file could not be found'); % later base this off check boxes
    return
end

if isempty(mf_fitter.user_inputs.reference_files.two)
    disp('Two-peak file could not be found');
    return
end

n1 = find(mf_fitter.numors == mf_fitter.user_inputs.reference_files.one);
n2 = find(mf_fitter.numors == mf_fitter.user_inputs.reference_files.two);

% Fitting first (one peak) reference file
mf_fitter_NEWcallbacks('set_free',1,n1); 
mf_fitter_NEWcallbacks('fit_GRASP',1,n1);
mf_fitter_NEWcallbacks('data_storage',1,n1);
pause(0.25);
close(mf_fitter.handles.plot_handle);

if ( get(mf_fitter.handles.three_peak,'Value') )
  % Fitting second (three peak) reference file
  mf_fitter_NEWcallbacks('set_free',3,n2);
  mf_fitter_NEWcallbacks('fit_GRASP',3,n2);
  mf_fitter_NEWcallbacks('data_storage',3,n2);
  pause(0.25);
  close(mf_fitter.handles.plot_handle);
else
    % Fitting second (two peak) reference file
    mf_fitter_NEWcallbacks('set_free',2,n2);
    mf_fitter_NEWcallbacks('fit_GRASP',2,n2);
    mf_fitter_NEWcallbacks('data_storage',2,n2);
    pause(0.25);
    close(mf_fitter.handles.plot_handle);
end

close(notification);

% old code to play "everything is awesome" song
% player = audioplayer(mf_fitter.awesome.y, mf_fitter.awesome.Fs);
% if get(mf_fitter.handles.awesome, 'Value') play(player); end


%% Perform Data Smoothing

if(get(mf_fitter.handles.smoothing.switch,'Value')  & isempty(mf_fitter.data.smoothed) )
    % smooth_prep formats the necessary inputs and then calls the smooth function
    mf_fitter_NEWcallbacks('smooth_prep')
elseif(mf_fitter.handles.smoothing.switch == 0)
    disp('Smoothing Switch is turned off');
elseif(isempty(mf_fitter.SmoothedData) == 0)
    disp('Data has already been smoothed with these parameters')
else
   disp('No smoothing has been performed');
end

fig_h = [];


%% Load Smoothed Data
% for use in every fitting cycle

phi = mf_fitter.data.smoothed.phi; 
Int = mf_fitter.data.smoothed.Int;
Int_err = mf_fitter.data.smoothed.Int_err;

% BUG - write table with smoothed data to export

%% Second Fitting Cycle
% Fix centers & fwhm, determine intensities

mf_fitter.algorithm_options.current_cycle = 'cycle2';
mf_fitter.algorithm_options.fitter_cycles = 2;
curCycName = 'cycle2';
prevCycName = 'cycle1';

% Grab fit data from reference files
ref1 = find(mf_fitter.numors == mf_fitter.user_inputs.reference_files.one);
ref2 = find(mf_fitter.numors == mf_fitter.user_inputs.reference_files.two);

% Fixed Values 
FWHM = (mf_fitter.fit_data.(prevCycName).fwhm(ref1)+mf_fitter.fit_data.(prevCycName).fwhm(ref2))/2;
XC1 = mf_fitter.fit_data.(prevCycName).center1(ref2) - status_flags.analysis_modules.sectors.theta;
XC2 = mf_fitter.fit_data.(prevCycName).center2(ref1) - status_flags.analysis_modules.sectors.theta;
XC3 = mf_fitter.fit_data.(prevCycName).center3(ref2) - status_flags.analysis_modules.sectors.theta;

% Guess Values (remember data has been smoothed & centered)
Y0 = (mf_fitter.fit_data.(prevCycName).background(ref1)+mf_fitter.fit_data.(prevCycName).background(ref2))/2;  % GRASP does one y0 per guassian peak
I01f = mf_fitter.fit_data.(prevCycName).intensity1(ref2);
I02i = mf_fitter.fit_data.(prevCycName).intensity2(ref1);
I03f = mf_fitter.fit_data.(prevCycName).intensity3(ref2);

for n=1:mf_fitter.depth    
    
    % change intensity guess as a function of n
    I01 = (n/mf_fitter.depth)*I01f;
    I02 = (mf_fitter.depth-n+1)/(mf_fitter.depth)*I02i; % +1 to avoid lower and upper bound being 0
    I03 = (n/mf_fitter.depth)*I03f;
    
    % order: y0, fwhm, i01, xc1, i02, xc2, i03, xc3
    % matlabFit(phi, int, interr, problemVarNames, problemVarValues, start, lower, upper)
    [F1 CI1 fig_h(n) chi2] = matlabFit(phi, Int(n,:), Int_err(n,:), {'fwhm','xc1', 'xc2', 'xc3'}, {FWHM, XC1, XC2, XC3}, [I01, I02, I03, Y0], [0, 0, 0, Y0], 1.25*[I01f, I02i, I03f, Y0]);
         
    % intensities and background free
    mf_fitter.fit_data.(curCycName).background(n,1) = F1.y0;
    mf_fitter.fit_data.(curCycName).background(n,2) = CI1(2,4) - F1.y0;
    mf_fitter.fit_data.(curCycName).intensity1(n,1) = F1.i01;
    mf_fitter.fit_data.(curCycName).intensity1(n,2) = CI1(2,1) - F1.i01;
    mf_fitter.fit_data.(curCycName).intensity2(n,1) = F1.i02;
    mf_fitter.fit_data.(curCycName).intensity2(n,2) = CI1(2,2) - F1.i02;
    mf_fitter.fit_data.(curCycName).intensity3(n,1) = F1.i03;
    mf_fitter.fit_data.(curCycName).intensity3(n,2) = CI1(2,3) - F1.i03;
    
    % centers and FWHM fixed, error should stay 0
    mf_fitter.fit_data.(curCycName).center1(n,1) = XC1;
    mf_fitter.fit_data.(curCycName).center2(n,1) = XC2;
    mf_fitter.fit_data.(curCycName).center3(n,1) = XC3; 
    mf_fitter.fit_data.(curCycName).fwhm(n,1) = FWHM;
    
    % chi2
    mf_fitter.fit_data.(curCycName).chi2(n) = chi2;
    
end

% unless user wants to view full cycle, close figures
if(view2) 
else close(fig_h) 
end
fig_h = [];
    
% Verify fits - BUG
mf_fitter_NEWcallbacks('FitCheck')

% export fits to folder
writetable(mf_fitter.handles.ExportTable, [fileName curCycName '_data.txt']);

try
    close(mf_fitter.handles.table) 
end 


%% Third Fitting Cycle
% Fix intensiites & fwhm, determine centers 
% no values from GRASP, no other corrections necessary
% values change with file number

mf_fitter.algorithm_options.current_cycle = 'cycle3';
mf_fitter.algorithm_options.fitter_cycles = 3;
curCycName = 'cycle3';
prevCycName = 'cycle2';

for n=1:mf_fitter.depth
    
    % Fixed Values
    FWHM = mf_fitter.fit_data.(prevCycName).fwhm(1);
    I01 = mf_fitter.fit_data.(prevCycName).intensity1(n);
    I02 = mf_fitter.fit_data.(prevCycName).intensity2(n);
    I03 = mf_fitter.fit_data.(prevCycName).intensity3(n);
    
    % Guess values
    Y0 = mf_fitter.fit_data.(prevCycName).background(n);      
    XC1 = mf_fitter.fit_data.(prevCycName).center1(n);
    XC2 = mf_fitter.fit_data.(prevCycName).center2(n);
    XC3 = mf_fitter.fit_data.(prevCycName).center3(n);

    % order: y0, fwhm, i01, xc1, i02, xc2, i03, xc3
    % matlabFit(phi, int, interr, problemVarNames, problemVarValues, start, lower, upper)
    [F2 CI2 fig_h(n) chi2] = matlabFit(phi, Int(n,:), Int_err(n,:), { 'fwhm', 'i01', 'i02', 'i03'}, {FWHM, I01, I02, I03}, [XC1, XC2, XC3, Y0], [1.5*XC1, -2, 0.5*XC3, 0.5*Y0], [0.5*XC1, 2, 1.5*XC3, 1.5*Y0])    
     
    % centers and background free
    mf_fitter.fit_data.(curCycName).background(n,1) = F2.y0;
    mf_fitter.fit_data.(curCycName).background(n,2) = CI2(2,4) - F2.y0;
    mf_fitter.fit_data.(curCycName).center1(n,1) = F2.xc1;
    mf_fitter.fit_data.(curCycName).center1(n,2) = CI2(2,1) - F2.xc1;
    mf_fitter.fit_data.(curCycName).center2(n,1) = XC2;
    mf_fitter.fit_data.(curCycName).center2(n,2) = CI2(2,2) - F2.xc1;
    mf_fitter.fit_data.(curCycName).center3(n,1) = F2.xc3; 
    mf_fitter.fit_data.(curCycName).center3(n,2) = CI2(2,3) - F2.xc3;
    
    % intensities and FWHM fixed, error should stay 0
    mf_fitter.fit_data.(curCycName).fwhm(n,1) = FWHM;
    mf_fitter.fit_data.(curCycName).intensity1(n,1) = I01;
    mf_fitter.fit_data.(curCycName).intensity2(n,1) = I02;
    mf_fitter.fit_data.intensity3(n,1) = I03;
    
    % chi2
    mf_fitter.fit_data.(curCycName).chi2(n) = chi2;
    
end

% unless user wants to view full cycle, close figures
if(view3) 
else close(fig_h) 
end
fig_h = [];

% Verify fits - BUG
mf_fitter_NEWcallbacks('FitCheck')
try
    close(mf_fitter.handles.table) 
end

%% Fourth Fitting Cycle
% Fix centers, intensities free, determine fwhm 
% But only peaks that meet minimum intensity cut-off
% no values from GRASP, no other corrections necessary
% values change with file number

mf_fitter.algorithm_options.current_cycle = 'cycle4';
mf_fitter.algorithm_options.fitter_cycles = 4;
curCycName = 'cycle4';
prevCycName = 'cycle3';


% BUG
% for extracting covariances when all data must be fit
% N = 1:length(mf_fitter.fit_data.cycles);

for n = 1:mf_fitter.depth
    
    % Fixed Values
    XC1 = mf_fitter.fit_data.(prevCycName).center1(n);
    XC2 = mf_fitter.fit_data.(prevCycName).center2(n);
    XC3 = mf_fitter.fit_data.(prevCycName).center3(n);

    % Guess values
    FWHM = mf_fitter.fit_data.(prevCycName).fwhm(n);
    Y0 = mf_fitter.fit_data.(prevCycName).background(n);  
    I01 = mf_fitter.fit_data.(prevCycName).intensity1(n);
    I02 = mf_fitter.fit_data.(prevCycName).intensity2(n);
    I03 = mf_fitter.fit_data.(prevCycName).intensity3(n);
    
    % order: y0, fwhm, i01, xc1, i02, xc2, i03, xc3
    % matlabFit(phi, int, interr, problemVarNames, problemVarValues, start, lower, upper)
    
    % original fit with centers fixed, for final algorithm
    [F3 CI3 fig_h(n) chi2] = matlabFit(phi, Int(N(n),:), Int_err(N(n),:), { 'xc1', 'xc2', 'xc3'}, {XC1, XC2, XC3}, [FWHM, I01, I02, I03, Y0], [0.4*FWHM, 0.5*I01, 0.5*I02, 0.5*I03, 0.8*Y0], [1.6*FWHM, 1.25*I01, 1.25*I02, 1.25*I03, 1.2*Y0])    
      
    %uncomment for convariance -  do not fix xc1 and xc3, just give initial values
    %[F3 CI3 fig_h(n)] = matlabFit(phi, Int(N(n),:), Int_err(N(n),:),{},{},[FWHM, I01, I02, I03, XC1, XC2, XC3, Y0], [0.4*FWHM, 0.4*I01, 0.4*I02, 0.4*I03, XC1-3, XC2-0.5, XC3-3, 0.2*Y0], [1.6*FWHM, 1.6*I01, 1.6*I02, 1.6*I03, XC1+3, XC2+0.5, XC3+3, 1.8*Y0])    
     
    % uncomment for completely free fit
    % [F3 CI3 fig_h(n)] = matlabFit(phi, Int(N(n),:), Int_err(N(n),:),{},{},[FWHM, FWHM, FWHM, I01, I02, I03, XC1, XC2, XC3, Y0], [0.4*FWHM, 0.4*FWHM, 0.4*FWHM, 0.4*I01, 0.4*I02, 0.4*I03, XC1-3, XC2-0.5, XC3-3, -Inf], [1.6*FWHM, 1.6*FWHM, 1.6*FWHM, 1.6*I01, 1.6*I02, 1.6*I03, XC1+3, XC2+0.5, XC3+3, Inf], 2)    
    
    % fwhm, intensities, and background free
    mf_fitter.fit_data.(curCycName).fwhm(n,1) = F3.fwhm;
    mf_fitter.fit_data.(curCycName).fwhm(n,2) = CI3(2,1) - F3.fwhm;
    mf_fitter.fit_data.(curCycName).intensity1(n,1) = F3.i01;
	mf_fitter.fit_data.(curCycName).intensity1(n,2) = CI3(2,2) - F3.i01;
	mf_fitter.fit_data.(curCycName).intensity2(n,1) = F3.i02;
    mf_fitter.fit_data.(curCycName).intensity2(n,2) = CI3(2,3) - F3.i02;
    mf_fitter.fit_data.(curCycName).intensity3(n,1) = F3.i03;
    mf_fitter.fit_data.(curCycName).intensity3(n,2) = CI3(2,4) - F3.i03;
    mf_fitter.fit_data.(curCycName).background(n,1) = F3.y0;
    mf_fitter.fit_data.(curCycName).background(n,2) = CI3(2,5) - F3.y0;
    
    % centers fixed, error should stay 0
    mf_fitter.fit_data.(curCycName).center1(n,1) = XC1;
    mf_fitter.fit_data.(curCycName).center2(n,1) = XC2;
    mf_fitter.fit_data.(curCycName).center3(n,1) = XC3;
    
     % uncomment when extracting covariance
%      mf_fitter.fit_data.center1(n,1) = F3.xc1;
%      mf_fitter.fit_data.center1(n,2) = CI3(2,5) - F3.xc1;
%      mf_fitter.fit_data.center2(n,1) = F3.xc2;
%      mf_fitter.fit_data.center2(n,2) = CI3(2,6) - F3.xc2;
%      mf_fitter.fit_data.center3(n,1) = F3.xc3;
%      mf_fitter.fit_data.center3(n,2) = CI3(2,7) - F3.xc3;
%      mf_fitter.fit_data.background(n,1) = F3.y0;
%      mf_fitter.fit_data.background(n,2) = CI3(2,8) - F3.y0;
     
     % uncomment for completely free 
%     mf_fitter.fit_data.fwhm1(n,1) = F3.fwhm1;
%     mf_fitter.fit_data.fwhm1(n,2) = CI3(2,1) - F3.fwhm1;
%     mf_fitter.fit_data.fwhm2(n,1) = F3.fwhm2;
%     mf_fitter.fit_data.fwhm2(n,2) = CI3(2,2) - F3.fwhm2;
%     mf_fitter.fit_data.fwhm3(n,1) = F3.fwhm3;
%     mf_fitter.fit_data.fwhm3(n,2) = CI3(2,3) - F3.fwhm3;
%     
%     mf_fitter.fit_data.intensity1(n,1) = F3.i01;
%     mf_fitter.fit_data.intensity1(n,2) = CI3(2,4) - F3.i01;
%     mf_fitter.fit_data.intensity2(n,1) = F3.i02;
%     mf_fitter.fit_data.intensity2(n,2) = CI3(2,5) - F3.i02;
%     mf_fitter.fit_data.intensity3(n,1) = F3.i03;
%     mf_fitter.fit_data.intensity3(n,2) = CI3(2,6) - F3.i03;
%     
%      mf_fitter.fit_data.center1(n,1) = F3.xc1;
%      mf_fitter.fit_data.center1(n,2) = CI3(2,7) - F3.xc1;
%      mf_fitter.fit_data.center2(n,1) = F3.xc2;
%      mf_fitter.fit_data.center2(n,2) = CI3(2,8) - F3.xc2;
%      mf_fitter.fit_data.center3(n,1) = F3.xc3;
%      mf_fitter.fit_data.center3(n,2) = CI3(2,9) - F3.xc3;
%      mf_fitter.fit_data.background(n,1) = F3.y0;
%      mf_fitter.fit_data.background(n,2) = CI3(2,10) - F3.y0;
     
    % chi2
    mf_fitter.fit_data.(curCycName).chi2(n) = chi2;
    
end

% unless user wants to view full cycle, close figures
if(view4) 
else close(fig_h) 
end
fig_h = [];

% Verify fits
mf_fitter_NEWcallbacks('FitCheck')   
try
    close(mf_fitter.handles.table) 
end  
   
    
%% Fifth Fitting Cycle
% Double Check centers with new fwhm
% no values from GRASP, no other corrections necessary
% values change with file number

mf_fitter.algorithm_options.current_cycle = 'cycle5';
mf_fitter.algorithm_options.fitter_cycles = 5;
curCycName = 'cycle5';
prevCycName = 'cycle4';

% Fixed Values
% BUG - later add code to descriminate against outliers

% determine fms to see which peaks meet the intensity criteria
fms = mf_fitter.fit_data.(prevCycName).intensity2(:,1) ./ (mf_fitter.fit_data.(prevCycName).intensity1(:,1) + mf_fitter.fit_data.(prevCycName).intensity2(:,1) + mf_fitter.fit_data.(prevCycName).intensity3(:,1));

% fms must be greater than int_cutoff (MS peak meets critera)
% fms must be less than 1 - 2*int_cutoff (ES peaks meet criteria)
int_cutoff = mf_fitter.user_inputs.int_cutoff;
FWHM = mean(mf_fitter.fit_data.(prevCycName).fwhm(find(fms > int_cutoff & fms < (1-2*int_cutoff)),:));
FWHM(2) = std(mf_fitter.fit_data.(prevCycName).fwhm(find(fms > int_cutoff & fms < (1-2*int_cutoff)),:));

for n=1:mf_fitter.depth

    % Guess Values
    Y0 = mf_fitter.fit_data.(prevCycName).background(n);  
    I01 = mf_fitter.fit_data.(prevCycName).intensity1(n,1);
    I02 = mf_fitter.fit_data.(prevCycName).intensity2(n,1);
    I03 = mf_fitter.fit_data.(prevCycName).intensity3(n,1);
    XC1 = mf_fitter.fit_data.(prevCycName).center1(n,1);
    XC2 = mf_fitter.fit_data.(prevCycName).center2(n,1);
    XC3 = mf_fitter.fit_data.(prevCycName).center3(n,1);

    % order: y0, fwhm, i01, xc1, i02, xc2, i03, xc3
    % matlabFit(phi, int, interr, problemVarNames, problemVarValues, start, lower, upper)
    [F4 CI4 fig_h(n) chi2] = matlabFit(phi, Int(n,:), Int_err(n,:), {'fwhm','i02', 'xc2'}, {FWHM(1), I02, XC2}, [I01, XC1, I03, XC3, Y0], [0.75*I01,  0.75*I03, 1.5*XC1, 0.5*XC3, 0.75*Y0], [1.4*I01, 1.4*I03, 0.5*XC1, 1.5*XC3, 1.25*Y0])    %xc1 is neg, lower bound needs to be more neg
         
    
    % BUG - do I want I02, XC2 fixed here?
    
    % background, intensities, and centers free
    mf_fitter.fit_data.(curCycName).background(n,1) = F4.y0;
    mf_fitter.fit_data.(curCycName).background(n,2) = CI4(2,5) - F4.y0;

    mf_fitter.fit_data.(curCycName).intensity1(n,1) = F4.i01;
    mf_fitter.fit_data.(curCycName).intensity1(n,2) = CI4(2,1) - F4.i01;
    mf_fitter.fit_data.(curCycName).intensity3(n,1) = F4.i03;
    mf_fitter.fit_data.(curCycName).intensity3(n,2) = CI4(2,2) - F4.i03;
    mf_fitter.fit_data.(curCycName).center1(n,1) = F4.xc1;
    mf_fitter.fit_data.(curCycName).center1(n,2) = CI4(2,3) - F4.xc1; 
    mf_fitter.fit_data.(curCycName).center3(n,1) = F4.xc3; 
    mf_fitter.fit_data.(curCycName).center3(n,2) = CI4(2,4) - F4.xc3;
    
    % FWHM fixed, error should be standard deviation
    mf_fitter.fit_data.(curCycName).fwhm(n,1) = FWHM(1);
    mf_fitter.fit_data.(curCycName).fwhm(n,2) = FWHM(2);
    mf_fitter.fit_data.(curCycName).intensity2(n,1) = I02;
    mf_fitter.fit_data.(curCycName).center2(n,1) = XC2;
    
    % chi2
    mf_fitter.fit_data.(curCycName).chi2(n) = chi2;
    
%     save = get(mf_fitter.handles.save,'Value')
%     if(save)
%        numor = mf_fitter.fit_data.names(n); 
%        mf_save('save',['I vs Xi - ' num2str(numor)],fig_h(n)) 
%     end
    
end

% unless user wants to view full cycle, close figures
if(view4) 
else close(fig_h) 
end
fig_h = [];

% Verify fits
mf_fitter_NEWcallbacks('FitCheck')   
try
    close(mf_fitter.handles.table) 
end  

%% Sixth Fitting Cycle
% Final fit
% values change with file number

mf_fitter.algorithm_options.current_cycle = 'cycle6';
mf_fitter.algorithm_options.fitter_cycles = 6;
curCycName = 'cycle6';
prevCycName = 'cycle5';

for r = 1:2
for n=1:mf_fitter.depth
    
    % Guess values (remember data has been smoothed & centered)  
    Y0 = mf_fitter.fit_data.(prevCycName).background(n,1);  
    I01 = mf_fitter.fit_data.(prevCycName).intensity1(n,1);
    I02 = mf_fitter.fit_data.(prevCycName).intensity2(n,1);
    I03 = mf_fitter.fit_data.(prevCycName).intensity3(n,1);
    XC1 = mf_fitter.fit_data.(prevCycName).center1(n,1);
    XC2 = mf_fitter.fit_data.(prevCycName).center2(n,1);
    XC3 = mf_fitter.fit_data.(prevCycName).center3(n,1);
    
    % Fixed values
    FWHM(1) = mf_fitter.fit_data.(prevCycName).fwhm(n,1);
    FWHM(2) = mf_fitter.fit_data.(prevCycName).fwhm(n,2);
    
    % order: y0, fwhm, i01, xc1, i02, xc2, i03, xc3
    % matlabFit(phi, int, interr, problemVarNames, problemVarValues, start, lower, upper)
    [F5 CI5 fig_h(n) chi2] = matlabFit(phi, Int(n,:), Int_err(n,:), {'fwhm'}, {FWHM(1)}, [I01, XC1, I02, XC2, I03, XC3, Y0], [0.5*I01, 0.5*I02, 0.5*I03, 1.5*XC1, XC2-1, 0.5*XC3, 0.75*Y0], [1.25*I01, 1.25*I02, 1.25*I03, 0.5*XC1, XC2+1, 1.5*XC3, 1.25*Y0])    %xc1 is neg, lower bound needs to be more neg
         
    % background, intensities, and centers free
    mf_fitter.fit_data.(curCycName).background(n,1) = F5.y0;
    mf_fitter.fit_data.(curCycName).background(n,2) = CI5(2,5) - F5.y0;
    
    mf_fitter.fit_data.(curCycName).intensity1(n,1) = F5.i01;
    mf_fitter.fit_data.(curCycName).intensity1(n,2) = CI5(2,1) - F5.i01;
    mf_fitter.fit_data.(curCycName).intensity2(n,1) = F5.i02;
    mf_fitter.fit_data.(curCycName).intensity1(n,2) = CI5(2,2) - F5.i02;
    mf_fitter.fit_data.(curCycName).intensity3(n,1) = F5.i03;
    mf_fitter.fit_data.(curCycName).intensity3(n,2) = CI5(2,3) - F5.i03;
    
    mf_fitter.fit_data.(curCycName).center1(n,1) = F5.xc1;
    mf_fitter.fit_data.(curCycName).center1(n,2) = CI5(2,4) - F5.xc1;
    mf_fitter.fit_data.(curCycName).center2(n,1) = F5.xc2;
    mf_fitter.fit_data.(curCycName).center2(n,2) = CI5(2,5) - F5.xc2;
    mf_fitter.fit_data.(curCycName).center3(n,1) = F5.xc3; 
    mf_fitter.fit_data.(curCycName).center3(n,2) = CI5(2,6) - F5.xc3;
    
    % FWHM fixed, error should be standard deviation
    mf_fitter.fit_data.fwhm(n,1) = FWHM(1);
    mf_fitter.fit_data.fwhm(n,2) = FWHM(2);
    
    % chi2
    mf_fitter.fit_data.(curCycName).chi2(n) = chi2;
    
    if(r == 1)
        close(fig_h(n))
    end
    
    save = get(mf_fitter.handles.save,'Value');
    if(save & r == 2)
       numor = mf_fitter.fit_data.names(n); 
       mf_save('save',['I vs Xi - ' num2str(numor)],fig_h(n));
    end
    
end    
end

mf_table;

    
%% Extra plots - check if user wants to make these
    save = get(mf_fitter.handles.save,'Value')
    % BUG - move to plot_Maker
    
    if(save)
        phi = mf_fitter.data.smoothed.phi; 
        Int = mf_fitter.data.smoothed.Int;
        Int_err = mf_fitter.data.smoothed.Int_err;
        
        
        %% Make Colormap
            cm_han = CM(phi, Int, cyc);
            mf_save('save','Colormap',cm_han)

            
        %% Peak Separation
            mf_fitter_callbacks('center_separation');
            y  = abs(mf_fitter.fit_data.(curCycName).center3(:,1) - mf_fitter.fit_data.(curCycName).center1(:,1));
            y_err = sqrt((mf_fitter.fit_data.(curCycName).center3(:,2)).^2+(mf_fitter.fit_data.(curCycName).center1(:,2)).^2);
            ps_han = figure;
            plotData(cyc, y, y_err, 'Applied || AC Cycles', '\Delta \phi (degrees)', ['Peak Separation - ' mf_fitter.folder])
            set(gca,'xScale','log')
            
            mf_save('save','PeakSeparation',ps_han)
        
        
        %% Intensity
            int_han = figure
            errorbar(cyc,mf_fitter.fit_data.(curCycName).intensity1(:,1),mf_fitter.fit_data.(curCycName).intensity1(:,2),'ok')
            hold on
            errorbar(cyc,mf_fitter.fit_data.(curCycName).intensity3(:,1),mf_fitter.fit_data.(curCycName).intensity3(:,2),'sk')
            title('Intensity')
            xlabel('Applied || AC Cycles')
            ylabel('Peak Intensity (arb. units)')
            set(gca,'xscale','log')
            legend('peak 1', 'peak 2')
            plot_template
            
            mf_save('save','PeakIntensity',int_han)
        
            
        %% Angle Decay
            set(mf_fitter.handles.plot_dropdown,'Value',6) 
            mf_GUI_callbacks( 'plot_dropdown' )
            h = gcf; 
            
            mf_save('save','AngleDecay',h)
            
            
    end
    
end
%
