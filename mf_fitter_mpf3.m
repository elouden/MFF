function [ ] = mf_fitter_mpf3( )
%   Reference Files:  fits a 1-peak and a 2-peak file, used to get initial guess values (GRASP)
%   First Cycle: Fix intensities & fwhm, determine centers (MATLAB)
%   Second Cycle:  Fix centers & fwhm, determine intensities (MATLAB)
%   Third Cycle:   Fix centers, intensities free, determine fwhm from files with dxc > 2*fwhm  (MATLAB)
%   Fourth Cycle:  Fix fwhm, centers & intensities free (MATLAB)
%                   repeated 3x to converge on final fit

%       **Note:  so far in all cycles, i02 == 0; xc2 ==0

% v 9.1
% 5/7/2017 MFF Liz


%% Initialize Global Variables & Functions script
global mf_fitter;
global status_flags;
global grasp_handles;

mf_fitter_callbacks('initialize');

mf_fitter.fitter.type = {'G','M','M','M','M','M'};
mf_fitter.fitter.cycles = [];

%% Check which Fitting Cycles are to be viewed
view1 = mf_fitter.cycleview.view1;
view2 = mf_fitter.cycleview.view2;
view3 = mf_fitter.cycleview.view3;
view4 = mf_fitter.cycleview.view4;
view5 = mf_fitter.cycleview.view5;


%% For simplicity, fit reference peaks in GRASP 

notification = msgbox('Fitting reference peaks...');
mf_fitter.fitter.cycles = 1;

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
mf_fitter_callbacks('data_storage',1,mf_fitter.good_peaks.inner_num);
pause(0.25);
close(mf_fitter.handles.plot_handle);
mf_fitter.flag = 0;


if ( get(mf_fitter.handles.three_peak,'Value') )
  mf_fitter_callbacks('set_free',3,mf_fitter.good_peaks.inner_num);
  mf_fitter_callbacks('fit',3,mf_fitter.good_peaks.outers_num);
  mf_fitter_callbacks('data_storage',3,mf_fitter.good_peaks.outers_num);
  pause(0.25);
  close(mf_fitter.handles.plot_handle);
else
    %Fitting second (two peak) reference file
    mf_fitter_callbacks('set_free',2,mf_fitter.good_peaks.inner_num);
    mf_fitter_callbacks('fit',2,mf_fitter.good_peaks.outers_num);
    mf_fitter_callbacks('data_storage',2,mf_fitter.good_peaks.outers_num);
    pause(0.25);
    close(mf_fitter.handles.plot_handle);
end

close(notification);


%% Perform Data Smoothing

if(mf_fitter.handles.smoothing.switch  & isempty(mf_fitter.SmoothedData) )
    disp(['Perform MRE smoothing with FWHM: ' num2str(mf_fitter.smoothing.fwhm)]);
    for i = 1:length(mf_fitter.fit_data.names)
        [phi, Int(i,:), Int_err(i,:)] = getData(i);
        
        
        phisym = status_flags.analysis_modules.sectors.theta;
        dTheta = status_flags.analysis_modules.sectors.delta_theta;
        ss = mf_fitter.smoothing.step;
        
        %phirng = [phisym - dTheta/2; phisym + dTheta/2; ss]   % because we are centering, do not need the phisym value
        phirng = [-dTheta/2; dTheta/2; ss]
        [phiS, IntS(i,:), Int_errS(i,:)] = smooth(phi, Int(i,:), Int_err(i,:), i, phisym, 0, phirng, mf_fitter.smoothing.fwhm)
        
        fig_h = plotData(phiS, IntS(i,:), Int_errS(i,:), '\phi  - \phi_0', 'Intensity (arb. units)', ['Smoothed Data Numor - ' num2str(i)]);
        if(view1)
            pause(1)
        end
        close(fig_h)
    end
    
    if(mf_fitter.smoothing.fwhm == 1)
        mf_fitter.SmoothedData.phi = phi - phisym;
        mf_fitter.SmoothedData.phi = mf_fitter.SmoothedData.phi'
        mf_fitter.SmoothedData.Int = Int;
        mf_fitter.SmoothedData.Int_err = Int_err;
        disp('data symmetrized but not smoothed')
    else   
        mf_fitter.SmoothedData.phi = phiS;
        mf_fitter.SmoothedData.Int = IntS;
        mf_fitter.SmoothedData.Int_err = Int_errS;
    end
    
    disp('Smoothing Complete');
elseif(mf_fitter.handles.smoothing.switch == 0)
    disp('Smoothing Switch is turned off');
elseif(isempty(mf_fitter.SmoothedData) == 0)
    disp('Data has already been smoothed with these parameters')
else
   disp('No smoothing has been performed');
end

%% Load Smoothed Data
% for use in every fitting cycle

phi = mf_fitter.SmoothedData.phi;
Int = mf_fitter.SmoothedData.Int;
Int_err = mf_fitter.SmoothedData.Int_err;


%% First Fitting Cycle
% Fix intensities & fwhm, determine centers
mf_fitter.fitter.cycles = 2;

% Use guess values from 2-peak reference file for most
ref1 = mf_fitter.good_peaks.outers_num;
ref2 = mf_fitter.good_peaks.inner_num;

% Fixed Values
FWHM = mf_fitter.fit_data.fwhm(ref1);
I01 = mf_fitter.fit_data.intensity1(ref1);
I02 = 0; %mf_fitter.fit_data.intensity2(1);
I03 = mf_fitter.fit_data.intensity3(ref1);

% Guess values (remember data has been smoothed & centered)
Y0 = 2 * mf_fitter.fit_data.background(ref1);  % GRASP does one y0 per guassian peak
XC1 = mf_fitter.fit_data.center1(ref1) - status_flags.analysis_modules.sectors.theta;
XC2 = mf_fitter.fit_data.center2(ref1);
XC3 = mf_fitter.fit_data.center3(ref1) - status_flags.analysis_modules.sectors.theta;

N = length(mf_fitter.fit_data.names);
for n=1:N
    % guess how close peaks end up 
    fwhm_f = mf_fitter.fit_data.fwhm(ref2);
    
    dxc = fwhm_f - FWHM;
    XC1 = XC1*((N+1)-n)/(N+1) + dxc/2;
    XC3 = XC3*((N+1)-n)/(N+1) + dxc/2;
    
    % order: y0, fwhm, i01, xc1, i02, xc2, i03, xc3
    % matlabFit(phi, int, interr, problemVarNames, problemVarValues, start, lower, upper)
    [F1 CI1 fig_h(n)] = matlabFit(phi, Int(n,:), Int_err(n,:), {'fwhm','i01', 'i02', 'xc2', 'i03'}, {FWHM, I01, I02, XC2, I03}, [XC1, XC3, Y0], [-15, 0, 0], [0, 15, 5*Y0]);
         
    mf_fitter.fit_data.background(n,1) = F1.y0;
    mf_fitter.fit_data.background(n,2) = CI1(2,3) - F1.y0;
    mf_fitter.fit_data.fwhm(n,1) = FWHM;
    mf_fitter.fit_data.intensity1(n,1) = I01;
    mf_fitter.fit_data.intensity2(n,1) = I02;
    mf_fitter.fit_data.intensity3(n,1) = I03;
    mf_fitter.fit_data.center1(n,1) = F1.xc1;
    mf_fitter.fit_data.center1(n,2) = CI1(2,1) - F1.xc1;
    mf_fitter.fit_data.center2(n,1) = XC2;
    mf_fitter.fit_data.center3(n,1) = F1.xc3; 
    mf_fitter.fit_data.center3(n,2) = CI1(2,2) - F1.xc3;
    
end

if(view2) 
else close(fig_h) 
end
fig_h = [];

% Verify fits
mf_fitter_NEWcallbacks('FitCheck')
try
    close(mf_fitter.handles.table) 
end 

%% Second Fitting Cycle
% Fix centers & fwhm, determine intensities (still keeping I02 = 0)
% no values from GRASP, no other corrections necessary
% values change with file number
mf_fitter.fitter.cycles = 3;

N = length(mf_fitter.fit_data.names);
for n=1:N
    % Fixed Values
    FWHM = mf_fitter.fit_data.fwhm(1);
    I02 = 0; %mf_fitter.fit_data.intensity2(1);
    XC1 = mf_fitter.fit_data.center1(n);
    XC2 = mf_fitter.fit_data.center2(1);
    XC3 = mf_fitter.fit_data.center3(n);

    % Guess values
    Y0 = mf_fitter.fit_data.background(1);  
    I01 = mf_fitter.fit_data.intensity1(1);
    I03 = mf_fitter.fit_data.intensity3(1);

    % order: y0, fwhm, i01, xc1, i02, xc2, i03, xc3
    % matlabFit(phi, int, interr, problemVarNames, problemVarValues, start, lower, upper)
    [F2 CI2 fig_h(n)] = matlabFit(phi, Int(n,:), Int_err(n,:), { 'fwhm', 'xc1', 'i02', 'xc2', 'xc3'}, {FWHM, XC1, I02, XC2, XC3}, [I01, I03, Y0], 0.25*[I01, I03, 0], 2*[I01, I03, (5/2)*Y0])    
         
    mf_fitter.fit_data.background(n,1) = F2.y0;
    mf_fitter.fit_data.background(n,2) = CI2(2,3) - F2.y0
    mf_fitter.fit_data.intensity1(n,1) = F2.i01;
    mf_fitter.fit_data.intensity1(n,2) = CI2(2,1) - F2.i01;
    mf_fitter.fit_data.intensity3(n,1) = F2.i03;
    mf_fitter.fit_data.intensity3(n,2) = CI2(2,2) - F2.i03;
end

    if(view3) 
    else close(fig_h) 
    end
    fig_h = [];

    
% Verify fits
mf_fitter_NEWcallbacks('FitCheck')
try
    close(mf_fitter.handles.table) 
end
    
%% Third Fitting Cycle
% Fix centers, intensities free, determine fwhm (still keeping I02 = 0)
% But only peaks separated by more than 2xfwhm 
% no values from GRASP, no other corrections necessary
% values change with file number
mf_fitter.fitter.cycles = 4;

dxc_c = 1.75 * mf_fitter.fit_data.fwhm(1,1);

% for final algorithm when above conditional is required
% N = find(abs(mf_fitter.fit_data.center1 - mf_fitter.fit_data.center3) > dxc_c);

% for extracting covariances when all data must be fit
 N = 1:length(mf_fitter.fit_data.cycles);

for n = 1:length(N)
    % Fixed Values
    I02 = 0; %mf_fitter.fit_data.intensity2(1);
    XC1 = mf_fitter.fit_data.center1(N(n));
    XC2 = mf_fitter.fit_data.center2(1);
    XC3 = mf_fitter.fit_data.center3(N(n));

    % Guess values
    FWHM = mf_fitter.fit_data.fwhm(1);
    Y0 = mf_fitter.fit_data.background(N(n));  
    I01 = mf_fitter.fit_data.intensity1(N(n));
    I03 = mf_fitter.fit_data.intensity3(N(n));
    
    
    % order: y0, fwhm, i01, xc1, i02, xc2, i03, xc3
    % matlabFit(phi, int, interr, problemVarNames, problemVarValues, start, lower, upper)
    
    % original fit with centers fixed, for final algorithm
    % [F3 CI3 fig_h(n)] = matlabFit(phi, Int(N(n),:), Int_err(N(n),:), { 'xc1', 'i02', 'xc2', 'xc3'}, {XC1, I02, XC2, XC3}, [FWHM, I01, I03, Y0], 0.4*[FWHM, I01, I03, Y0], 1.6*[FWHM, I01, I03, Y0])    
     
    %uncomment for convariance -  do not fix xc1 and xc3, just give initial values
    % [F3 CI3 fig_h(n)] = matlabFit(phi, Int(N(n),:), Int_err(N(n),:), { 'i02', 'xc2'}, {I02, XC2}, [FWHM, I01, I03, XC1, XC2, Y0], [0.4*FWHM, 0.4*I01, 0.4*I03, XC1-3, XC3-3, 0.4*Y0], [1.6*FWHM, 1.6*I01, 1.6*I03, XC1+3, XC3+3, 1.6*Y0])    
     
    % uncomment for completely free fit
    [F3 CI3 fig_h(n)] = matlabFit(phi, Int(N(n),:), Int_err(N(n),:),{},{},[FWHM, FWHM, I01/(I03+I01), I03+I01, XC1, XC3, Y0], [0.3*FWHM, 0.3*FWHM, 0.4*I01/(I03+I01), 0.4*(I03+I01), XC1-3, XC3-3, 0], [2*FWHM, 2*FWHM, 1.6*I01/(I03+I01), 1.6*(I03+I01), XC1+3, XC3+3, 4*Y0], 3);    
    
    
    
%     mf_fitter.fit_data.fwhm(n,1) = F3.fwhm;
%     mf_fitter.fit_data.fwhm(n,2) = CI3(2,1) - F3.fwhm;
%     mf_fitter.fit_data.intensity1(n,1) = F3.i01;
%     mf_fitter.fit_data.intensity1(n,2) = CI3(2,2) - F3.i01;
%     mf_fitter.fit_data.intensity3(n,1) = F3.i03;
%     mf_fitter.fit_data.intensity3(n,2) = CI3(2,3) - F3.i03;
    
    % uncomment when extracting covariance
%     mf_fitter.fit_data.center1(n,1) = F3.xc1;
%     mf_fitter.fit_data.center1(n,2) = CI3(2,4) - F3.xc1;
%     mf_fitter.fit_data.center2(n,1) = F3.xc2;
%     mf_fitter.fit_data.center2(n,2) = CI3(2,5) - F3.xc2;
%     mf_fitter.fit_data.background(n,2) = CI3(2,6) - F3.y0;
%     mf_fitter.fit_data.fwhm(n,1) = F3.fwhm;

    % comment when extracting covariance/open fits
%     mf_fitter.fit_data.background(n,2) = CI3(2,4) - F3.y0;
%     mf_fitter.fit_data.fwhm(n,1) = F3.fwhm;

    %uncomment for open fit
    mf_fitter.fit_data.fwhm1(n,1) = F3.fwhm1;
    mf_fitter.fit_data.fwhm1(n,2) = CI3(2,1) - F3.fwhm1;
    mf_fitter.fit_data.fwhm2(n,1) = F3.fwhm3;
    mf_fitter.fit_data.fwhm2(n,2) = CI3(2,2) - F3.fwhm3;
    
    mf_fitter.fit_data.intensity1(n,1) = F3.i1*F3.itot;
    mf_fitter.fit_data.intensity1(n,2) = sqrt((CI3(2,3)*F3.itot)^2 + (CI3(2,4)*F3.i1)^2);
    mf_fitter.fit_data.intensity3(n,1) = (1-F3.i1)*F3.itot;
    mf_fitter.fit_data.intensity3(n,2) =  sqrt((CI3(2,3)*F3.itot)^2 + (CI3(2,4)*(1-F3.i1))^2);
    
    mf_fitter.fit_data.center1(n,1) = F3.xc1;
    mf_fitter.fit_data.center1(n,2) = CI3(2,5) - F3.xc1;    
    mf_fitter.fit_data.center3(n,1) = F3.xc3;
    mf_fitter.fit_data.center3(n,2) = CI3(2,6) - F3.xc3;
    
    mf_fitter.fit_data.background(n,1) = F3.y0;
    mf_fitter.fit_data.background(n,2) = CI3(2,7) - F3.y0;
end

    if(view4) 
    else close(fig_h) 
    end
    fig_h = [];

    
% Verify fits
mf_fitter_NEWcallbacks('FitCheck')   
try
    close(mf_fitter.handles.table) 
end  
   
disp('pause to save fwhm - under fit_data.fwhm')

%% Fourth Fitting Cycle
% Double Check centers with new fwhm
% no values from GRASP, no other corrections necessary
% values change with file number
mf_fitter.fitter.cycles = 5;

% Fixed Values
% later add code to descriminate against outliers
nonZero = find(mf_fitter.fit_data.fwhm(:,2) > 0); % indicates that is where a fit was performed (it has an associated error)
FWHM = mean(mf_fitter.fit_data.fwhm(nonZero,:));
I02 = 0;
XC2 = mf_fitter.fit_data.center2(1);

N = length(mf_fitter.fit_data.names);
Nc = max(find(abs(mf_fitter.fit_data.center1 - mf_fitter.fit_data.center3) > dxc_c));

for n=1:N
    % Guess values  
        Y0 = mf_fitter.fit_data.background(n);  
    if(n<Nc)
        I01 = mf_fitter.fit_data.intensity1(n,1);
        I03 = mf_fitter.fit_data.intensity3(n,1);
        XC1 = mf_fitter.fit_data.center1(n,1);
        XC3 = mf_fitter.fit_data.center3(n,1);
    else
        I01 = mf_fitter.fit_data.intensity1(Nc,1);
        I03 = mf_fitter.fit_data.intensity3(Nc,1);
        XC1 = mf_fitter.fit_data.center1(Nc,1);
        XC3 = mf_fitter.fit_data.center3(Nc,1);
    end
    
    % use new fit type (opt = 1) to better determine the total intenisty
    ITOT = I01 + I03;
    I1 = I01 / (ITOT);

     % ORDER
     % i1, itot,  xc1, xc3, y0
     %matlabFit(phi, int, interr, problemVarNames, problemVarValues, start, lower, upper, opt)
    [F4 CI4 fig_h(n)] = matlabFit(phi, Int(n,:), Int_err(n,:), {'fwhm'}, {FWHM(1)}, [I1, ITOT, XC1, XC3, Y0], [ 0.75*I1, 0.8*ITOT, 1.5*XC1, 0.5*XC3, 0.5*Y0], [1.25*I1, 1.5*ITOT, 0.5*XC1, 1.5*XC3, 1.5*Y0], 1)    %xc1 is neg, lower bound needs to be more neg
         
     mf_fitter.fit_data.background(n,1) = F4.y0;
     mf_fitter.fit_data.background(n,2) = CI4(2,5) - F4.y0;
     mf_fitter.fit_data.fwhm(n,1) = FWHM(1);
     mf_fitter.fit_data.fwhm(n,2) = FWHM(2);
     mf_fitter.fit_data.center1(n,1) = F4.xc1;
     mf_fitter.fit_data.center1(n,2) = CI4(2,3) - F4.xc1;
     mf_fitter.fit_data.center2(n,1) = XC2;
     mf_fitter.fit_data.center3(n,1) = F4.xc3; 
     mf_fitter.fit_data.center3(n,2) = CI4(2,4) - F4.xc3;
     
     % intensity calculations are a little more complicated, also store itot and i1
     itot = F4.itot;
     itot_err = CI4(2,2) - itot;
     mf_fitter.fit_data.I_tot(n,1) = itot;
     mf_fitter.fit_data.I_tot(n,2) = itot_err;
     
     i1 = F4.i1;
     i1_err = CI4(2,1) - i1;
     mf_fitter.fit_data.I1(n,1) = i1;
     mf_fitter.fit_data.I1(n,2) = i1_err;
     
     mf_fitter.fit_data.intensity1(n,1) = itot*i1;
     mf_fitter.fit_data.intensity1(n,2) = sqrt((i1*itot_err)^2+(itot*i1_err)^2);
     mf_fitter.fit_data.intensity3(n,1) = itot*(1-i1);
     mf_fitter.fit_data.intensity3(n,2) = sqrt(((1-i1)*itot_err)^2+(itot*(-1)*i1_err)^2);

end

    
    if(view5) 
    else close(fig_h) 
    end
    fig_h = [];

    
% Verify fits
mf_fitter_NEWcallbacks('FitCheck')
try
    close(mf_fitter.handles.table) 
end 
    
%% Fifth Fitting Cycle
% Final fit (essentially the same as the fourth, just iterated 2x more)
% values change with file number
mf_fitter.fitter.cycles = 6;

phi = mf_fitter.SmoothedData.phi;
Int = mf_fitter.SmoothedData.Int;
Int_err = mf_fitter.SmoothedData.Int_err;

% Fixed Values
% later add code to descriminate against outliers
nonZero = find(mf_fitter.fit_data.fwhm(:,2) > 0) % indicates that is where a fit was performed (it has an associated error)
FWHM = mean(mf_fitter.fit_data.fwhm(nonZero,:))
I02 = 0;
XC2 = mf_fitter.fit_data.center2(1);

N = length(mf_fitter.fit_data.names);

for r = 1:2
for n=1:N
    % Guess values (remember data has been smoothed & centered) 
        Y0 = mf_fitter.fit_data.background(n,1); 
        XC1 = mf_fitter.fit_data.center1(n,1);
        XC3 = mf_fitter.fit_data.center3(n,1);
        ITOT = mf_fitter.fit_data.I_tot(n,1);
        I1 = mf_fitter.fit_data.I1(n,1);
    

     % ORDER
     % i1, itot,  xc1, xc3, y0
     %matlabFit(phi, int, interr, problemVarNames, problemVarValues, start, lower, upper, opt)
     [F5 CI5 fig_h(n)] = matlabFit(phi, Int(n,:), Int_err(n,:), {'fwhm'}, {FWHM(1)}, [I1, ITOT, XC1, XC3, Y0], [ 0.75*I1, 0.75*ITOT, 1.5*XC1, 0.5*XC3, 0.5*Y0], [1.25*I1, 1.5*ITOT, 0.5*XC1, 1.5*XC3, 1.5*Y0], 1)    %xc1 is neg, lower bound needs to be more neg
         
     mf_fitter.fit_data.background(n,1) = F5.y0;
     mf_fitter.fit_data.background(n,2) = CI5(2,5) - F5.y0;
     mf_fitter.fit_data.fwhm(n,1) = FWHM(1);
     mf_fitter.fit_data.fwhm(n,2) = FWHM(2);
     mf_fitter.fit_data.center1(n,1) = F5.xc1;
     mf_fitter.fit_data.center1(n,2) = CI5(2,3) - F5.xc1;
     mf_fitter.fit_data.center2(n,1) = XC2;
     mf_fitter.fit_data.center3(n,1) = F5.xc3; 
     mf_fitter.fit_data.center3(n,2) = CI5(2,4) - F5.xc3;

     itot = F5.itot;
     itot_err = CI5(2,2) - itot;
     mf_fitter.fit_data.I_tot(n,1) = itot;
     mf_fitter.fit_data.I_tot(n,2) = itot_err;
     
     i1 = F5.i1;
     i1_err = CI5(2,1) - i1;
     mf_fitter.fit_data.I1(n,1) = i1;
     mf_fitter.fit_data.I1(n,2) = i1_err;
     
     mf_fitter.fit_data.intensity1(n,1) = itot*i1;
     mf_fitter.fit_data.intensity1(n,2) = sqrt((i1*itot_err)^2+(itot*i1_err)^2);
     mf_fitter.fit_data.intensity3(n,1) = itot*(1-i1);
     mf_fitter.fit_data.intensity3(n,2) = sqrt(((1-i1)*itot_err)^2+(itot*(-1)*i1_err)^2);

    
    if(r == 1)
        close(fig_h(n))
    end
    
    save = get(mf_fitter.handles.save,'Value');
    if(save & r == 2)
       numor = mf_fitter.fit_data.names(n); 
       mf_fitter_save('save',['I vs Xi - ' num2str(numor)],fig_h(n));
    end
    
end    
end


% Verify fits
mf_fitter_NEWcallbacks('FitCheck')
    
    
if(save)
    fileName = [mf_fitter.extension mf_fitter.folder '_Buzz/' 'data.txt'];
    writetable(mf_fitter.handles.ExportTable, fileName);
end 

%% Extra plots - check if user wants to make these

    save = get(mf_fitter.handles.save,'Value')
    
    if(save)
        % mf_fitter_plotMaker()
        
        phi = mf_fitter.SmoothedData.phi;
        Int = mf_fitter.SmoothedData.Int;
        cyc = mf_fitter.fit_data.cycles;

        %% Make Colormap
            cm_han = CM(phi, Int, cyc);
            mf_fitter_save('save','Colormap',cm_han)

        %% Peak Separation
            mf_fitter_callbacks('center_separation');
            y  = abs(mf_fitter.fit_data.center3(:,1) - mf_fitter.fit_data.center1(:,1));
            y_err = sqrt((mf_fitter.fit_data.center3(:,2)).^2+(mf_fitter.fit_data.center1(:,2)).^2);
            ps_han = figure;
            plotData(cyc, y, y_err, 'Applied || AC Cycles', '\Delta \phi (degrees)', ['Peak Separation - ' mf_fitter.folder])
            set(gca,'xScale','log')
            
            mf_fitter_save('save','PeakSeparation',ps_han)
        
        
        %% Intensity
            int_han = figure
            errorbar(cyc,mf_fitter.fit_data.intensity1(:,1),mf_fitter.fit_data.intensity1(:,2),'ok')
            hold on
            errorbar(cyc,mf_fitter.fit_data.intensity3(:,1),mf_fitter.fit_data.intensity3(:,2),'sk')
            title('Intensity')
            xlabel('Applied || AC Cycles')
            ylabel('Peak Intensity (arb. units)')
            set(gca,'xscale','log')
            legend('peak 1', 'peak 2')
            plot_template
            
            mf_fitter_save('save','PeakIntensity',int_han)
            
            int_tot_han = figure
            errorbar(cyc,mf_fitter.fit_data.I_tot(:,1),mf_fitter.fit_data.I_tot(:,2),'ok')
            %hold on
            %errorbar(cyc,mf_fitter.fit_data.intensity3(:,1),mf_fitter.fit_data.intensity3(:,2),'sk')
            title('Total Intensity')
            xlabel('Applied || AC Cycles')
            ylabel('Total Intensity (arb. units)')
            set(gca,'xscale','log')
            %legend('peak 1', 'peak 2')
            plot_template
            
            mf_fitter_save('save','TotalIntensity',int_tot_han)
        
        %% Angle Decay
            set(mf_fitter.handles.plot_dropdown,'Value',6) 
            mf_GUI_callbacks( 'plot_dropdown' )
            h = gcf; 
            
            mf_fitter_save('save','AngleDecay',h)
    else
    end
    
end

% 
% % %% Functions to be Called
% % 
% function [phi, Int, Int_err] = getData(img_num)
%         global status_flags
%         global mf_fitter
%         global plot_info
%         global grasp_handles
%     
%         curAB = get(grasp_handles.window_modules.radial_average.azimuth_bin,'String');
%        
%         % load desired depth file and update main grasp GUI
%         % the +1 is necessary as file 1 represents a sum
%         status_flags.selector.fd = img_num+1;
%         set(grasp_handles.window_modules.radial_average.azimuth_bin,'String',0.1);
%         radial_average_callbacks2('azimuth_bin');
%         grasp_update;
% 
%         % plot I vs xi and create the current figure, store handle
%         radial_average_callbacks('averaging_control','azimuthal');
%         mf_fitter.handles.plot_handle = gcf;
%         close(gcf)
%         grasp_update;
%             
%         % save data into appropriate variables
%         phi = plot_info.export_data(:,1);
%         Int = plot_info.export_data(:,2);
%         Int_err = plot_info.export_data(:,3);
%         
%         % return to original selected angle binning
%         set(grasp_handles.window_modules.radial_average.azimuth_bin,'String',curAB)
%         radial_average_callbacks2('azimuth_bin');
%         grasp_update;
% 
% end
% 
% function [h] = plotData(x,y,z, xName, yName, titleName)
%     h = figure
%     errorbar(x,y,z,'ok')
%     
%     xlabel(xName)
%     ylabel(yName)
%     title(titleName)
%     
%     plot_template
% end
% 
% function [] = plot_template()
%     scale = 1.5;
%     fig_size = [7.3 5.7];
%     fig_position = [1.25 1.25];
%     position = scale*[fig_position fig_size];
% 
%     set(gca,...
%         'Color',[1 1 1],...
%         'FontName','Arial','FontSize',scale*8,...
%         'Xcolor',[0 0 0],'Ycolor',[0 0 0],...
%         'XMinorTick','on','YMinorTick','on',...
%         'Box','on');
%         %'XLim',xlim,'XMinorTick', 'on',...
%         %'YLim',ylim,'YMinorTick', 'on',...
%         %'Units','centimeters','Position', position,...
% 
% hold off
% 
% end
% 
% function [fitobject ci fig_h] = matlabFit(phi, int, interr, problemVarNames, problemVarValues, start, lower, upper)
%         % ORDER
%         % y0, fwhm, i01,  xc1, i02, xc2, i03, xc3
% 
% % 
%          fo = fitoptions('Method','NonlinearLeastSquares','StartPoint',start,'Lower',lower,'Upper',upper);
% % 
%          ft = fittype('y0 + (i01/(fwhm * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(x-xc1)^2/(fwhm^2/log(4))) + (i02/(fwhm * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(x-xc2)^2/(fwhm^2/log(4))) + (i03/(fwhm * sqrt(pi/2) / sqrt(log(4))))*exp(-2*(x-xc3)^2/(fwhm^2/log(4)))',...
%                          'options',fo,'problem',problemVarNames);        
% 
% %        
%          [fitobject, gof] = fit(phi', int', ft, 'problem', problemVarValues)
%          ci = confint(fitobject);
%          fig_h = plotData(phi, int, interr, '\phi - \phi_0', 'Intensity (arb. units)', 'Fit Cycle 1')
%          hold on
%          plot(fitobject)
%          hold off
% %   
% end
% 
% function [] = plotXC(cm_han, cyc, xc1, xc1_err, xc3, xc3_err)
% figure(cm_han)
% hold on
% 
% errorbar(cyc, xc1, xc1_err, 'ow')
% errorbar(cyc, xc3, xc3_err, 'ow')
% end
% 
% function [h] = CM(phi, Int, cm_x)
% global mf_fitter
% 
%   indicesToDelete = (find(abs(phi) > 31)); % change to reference phi ring later
%   sizeInt = size(Int)
%   for i = 1:length(indicesToDelete)
%       for j = 1:sizeInt(1)
%       Int(j,indicesToDelete(i)) = min(Int(j,:))
%       end
%   end
%     
%    % data already smoothed and stored in 2D array above
%    % just swap so that x corresponds to rows
%    Int = Int';
%    
%    % offset phi to account for the way pcolor positions pixels
%     %global phi_shift
%     phi_shift = phi;
%     for i=1:(length(phi)-1)
%         delta_phi(i) = (phi(i+1) - phi(i))/2;
%         phi_shift(i) = phi(i) - delta_phi(i);
%     end  
%     
%     % perform the background subtraction and intensity integration (=1) 
%     ang = numel(phi);
%     
%     N = length(mf_fitter.fit_data.names);
%     
%     for i=1:N
%         % Find the average background 
%             %bg =  sum(Int(1:3,i) + Int((ang-2):ang,i)) / 6;
%             %backsub_int(:,i) = Int(:,i) - bg;
%             
%             %new BG sub method
%             % grab the (25% of total) minimum values and average them for BG
%             temp_int = sort(Int(:,i));
%             BG = temp_int(1:0.25*numel(temp_int));
%             BG = sum(BG)/numel(BG);
%             
%             backsub_int(:,i) = Int(:,i) - BG;
% 
%         % Make the integrated intensity equal to one across every column
%             %s = sum(backsub_int(:,i));
%             %integrated_int(:,i) = backsub_int(:,i)/s;
%             % check
%             %sum(integrated_int(:,i))
%         
%         % Normalize to the maximum in each column
%             m = max(backsub_int(:,i));
%             norm_int(:,i) = backsub_int(:,i)/m;
%             
%     end
% 
%     
%     % For Logarithmic intensity scale    
%         %log_i = log(integrated_int-min(min(integrated_int)));
%         log_i = log(norm_int-min(min(norm_int)));
%         %log_i = log(norm_int)
%         
%         
%     % Add junk x and phi value (based on the way matlab command pcolor works)
%         cm_x(N+1) = 2*cm_x(N);
%         phi_shift(ang+1) = phi_shift(ang) + delta_phi(ang-1);
%         phi(ang+1) = phi(ang) + delta_phi(ang-1);
%         
%         
%    % plot the colormap
%    intName = {'Int','Smoothed Intensity','log_i','Logarithmic Intensity Scale','norm_int','Normalized to Max Per Field'}
%    for i = 2:2
%    %for i = 3
%        h = figure
%        
%        % get intensity data
%        cData = eval(char(intName(2*i-1)));
%        
%        % Add junk Int column, row
%        %size(cData)
%        %ang
%        %length(phi)
%        %length(cm_x)
%        cData(2:(ang+1),:) = cData;
%        %cData(1,:) = 0; % shows up as red
%        %cData(ang,:) = 0; % shows up as red
%        cData(ang+1,:) = 0;
%        %cData(:,2) = 0; % shows up as red
%        cData(:,N+1) = 0;
%        %cData(:,N) = 0; % shows up as red
%        
%        % make the colormap and label axes
%        pcolor(cm_x,phi_shift,cData)
%        title(intName(2*i))
%        xlabel('H(T)')
%        ylabel('\phi  - \phi_0  (degrees)')
%        
%        CM_lim = get(gca,'CLim')
%        sort_int = sort(cData(:))
%        l = ceil(length(sort_int(:))/10)
%        CM_lim(1) = mean(sort_int(2:l)) % always one value that is - inf, want to ignore this one
%        %CM_lim = [min(min())]
%        CM_lim(2) = max(max(cData))
%        set(gca,'XScale','log','CLim',CM_lim)
%        
%        plot_template
%        
%    end
%    
% end
% 
% 
% function [phiS, intS, intErrS] = smooth( nphi, nInt, nIntErr, fn, phisym, symtrz, phirng, gwidth) %, intrng, savfig, clsfig, hl, fig_h, cm, cm_x, path)   
% %% Allocating memory for processed data
% N = numel(fn);
% phi = phirng(1):phirng(3):phirng(2);
% Int = zeros([N length(phi)]);
% IntErr = Int;
%       
% for n = 1:N
% % center phi range around zero (subtract 360 from values > 180)
% %     data = readtable([pathname char(fn(n))],'HeaderLines',hl(n),'Delimiter','\t'); 
% %     nphi = data{:,{'Az_Angle'}};
% %     nInt = data{:,{'I'}};
% %     nIntErr = data{:,{'Err_I'}};
% %     
% %     %center about I-phase position (stored in phisym)
% 
% %      nphi = nphi - phisym
% %       
%      phiabove = find(nphi > 180);
% %     phisym
%      if(phisym(n) > 180)
%          j = 1;
%         if(length(phisym) == 1) pnum = 1; else pnum = n;   end
%          nphi = nphi - phisym(pnum)
%      else
%          nphi(phiabove) = nphi(phiabove) - 360;
% %     
%     %center about I-phase position (stored in phisym)
%      j = 0;
% % %     if(phisym > 180)
% % %         j = 1;
% % %     end
%          if(length(phisym) == 1) pnum = 1; else pnum = n;   end
% %         phisym(pnum)
%          max(nphi)
%          min(nphi)
%          nphi = nphi - phisym(pnum) + j*360;
%      end
% %     
%     % Symmetrize
%     if(symtrz);
%         nphi = [nphi-phisym(n); -(nphi-phisym(n))];
%         nInt = [nInt; nInt];
%         nIntErr = [nIntErr; nIntErr];
%     end
%     
%      % Smooth
%      WhtInt = Int(1,:);
%      WhtIntErrSq = WhtInt;
%      Wht = WhtInt;
% %     if(phisym > 180)
% %         disp('break')
% %     end
%      for i = 1:length(phi)
%          for j = 1:length(nphi)
%              WhtInt(j) = exp(-(phi(i)-nphi(j))^2/gwidth^2)*nInt(j);
%              WhtIntErrSq(j) = (exp(-(phi(i)-nphi(j))^2/gwidth^2)*nIntErr(j))^2;
%              Wht(j) = exp(-(phi(i)-nphi(j))^2/gwidth^2);
%          end
%          length(fn)
%          phisym
%          if(phisym(n) > 180 & length(fn)>1)
%              disp('break')
%          end
%          Int(n,i) = sum(WhtInt)/sum(Wht);
%          IntErr(n,i) = sqrt(sum(WhtIntErrSq))/sum(Wht);
%      end
% %     
%        phiS = phi;
%        intS = Int;
%        intErrS = IntErr;
% %     %Store in data structure
% %     name = ['Field_' char(fh(n))];
% %     ILL_Dat.(name).phi = phi;
% %     ILL_Dat.(name).Int = Int(n,:);
% %     ILL_Dat.(name).IntErr = IntErr(n,:);
% %     
%  end
% 
% end