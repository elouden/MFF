function [ ] = mf_fitter_mpf2K( )
    
%   Reference Files:  fits a 1-peak and a 2-peak file, used to get initial guess values (GRASP)
%   First Cycle: Fix centers & fwhm, determine intesnities (MATLAB)
%   Second Cycle:  % Fix intensiites & fwhm, determine centers  (Matlab)
%   Third Cycle:   Fix centers, intensities free, determine fwhm from files with dxc > 2*fwhm  (MATLAB)
%   Fourth Cycle:  Fix fwhm, centers & intensities free (MATLAB)

%       

% v 9.0
% 11/9/2016 MFF Liz

% v 9.0
% 11/18/2016 MFF Liz

%% INITIALIZE global vars
global mf_fitter;
global status_flags;
global grasp_handles;

mf_fitter_callbacks('initialize');


%number of fitting cycles
Nc = 6;
wait_position = [0.85 .2 .25 .08];
refit_check = 1;


%% Check which Fitting Cycles are to be viewed
view1 = mf_fitter.cycleview.view1;
view2 = mf_fitter.cycleview.view2;
view3 = mf_fitter.cycleview.view3;
view4 = mf_fitter.cycleview.view4;
view5 = mf_fitter.cycleview.view5;


%% For simplicity, fit reference peaks in GRASP 

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

% old code to play "everything is awesome" song
% player = audioplayer(mf_fitter.awesome.y, mf_fitter.awesome.Fs);
% if get(mf_fitter.handles.awesome, 'Value') play(player); end


%% Perform Data Smoothing

if(get(mf_fitter.handles.smoothing.switch,'Value')  & isempty(mf_fitter.SmoothedData) )
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

fig_h = []

%% First Fitting Cycle
% Fix centers & fwhm, determine intensities


N = length(mf_fitter.fit_data.names);

phi = mf_fitter.SmoothedData.phi;
Int = mf_fitter.SmoothedData.Int;
Int_err = mf_fitter.SmoothedData.Int_err;

% Fixed Values
FWHM = mf_fitter.fit_data.fwhm(1);
XC1 = mf_fitter.fit_data.center1(N) - status_flags.analysis_modules.sectors.theta;
XC2 = mf_fitter.fit_data.center2(1) - status_flags.analysis_modules.sectors.theta;
XC3 = mf_fitter.fit_data.center3(N) - status_flags.analysis_modules.sectors.theta;


% Guess values (remember data has been smoothed & centered)
Y0 = 2 * mf_fitter.fit_data.background(1);  % GRASP does one y0 per guassian peak
I01 = mf_fitter.fit_data.intensity1(N);
I02 = mf_fitter.fit_data.intensity2(1);
I03 = mf_fitter.fit_data.intensity3(N);

for n=1:N
    % guess how close peaks end up 
    %fwhm_f = mf_fitter.fit_data.fwhm(N);
    
    %dxc = fwhm_f - FWHM;
    %XC1 = XC1*((N+1)-n)/(N+1) + dxc/2;
    %XC3 = XC3*((N+1)-n)/(N+1) + dxc/2;
    
    % order: y0, fwhm, i01, xc1, i02, xc2, i03, xc3
    % matlabFit(phi, int, interr, problemVarNames, problemVarValues, start, lower, upper)
    n
    [F1 CI1 fig_h(n)] = matlabFit(phi, Int(n,:), Int_err(n,:), {'fwhm','xc1', 'xc2', 'xc3'}, {FWHM, XC1, XC2, XC3}, [I01, I02, I03, Y0], 0.25*[I01, I02, I03, Y0], 1.75*[I01, I02, I03, Y0]);
         
    mf_fitter.fit_data.background(n,1) = F1.y0;
    mf_fitter.fit_data.background(n,2) = CI1(2,4) - F1.y0;
    mf_fitter.fit_data.fwhm(n,1) = FWHM;
    mf_fitter.fit_data.intensity1(n,1) = F1.i01;
    mf_fitter.fit_data.intensity1(n,2) = CI1(2,1) - F1.i01;
    mf_fitter.fit_data.intensity2(n,1) = F1.i02;
    mf_fitter.fit_data.intensity2(n,2) = CI1(2,2) - F1.i02;
    mf_fitter.fit_data.intensity3(n,1) = F1.i03;
    mf_fitter.fit_data.intensity3(n,2) = CI1(2,3) - F1.i03;
    mf_fitter.fit_data.center1(n,1) = XC1;
    %mf_fitter.fit_data.center1(n,2) = CI1(2,1) - F1.xc1;
    %mf_fitter.fit_data.center1(img_num,2) = F2;
    mf_fitter.fit_data.center2(n,1) = XC2;
    %mf_fitter.fit_data.center2(img_num,2) = %((fitobject.a-ci(1)) + (ci(2)-fitobject.a))/2; % + (((fitobject.b-ci(3))+(ci(4)-fitobject.b))/2)^2 );
    mf_fitter.fit_data.center3(n,1) = XC3; 
    %mf_fitter.fit_data.center3(n,2) = CI1(2,2) - F1.xc3;
    
end

    if(view2) 
    else close(fig_h) 
    end
    fig_h = [];
    
    
%% Second Fitting Cycle
% Fix intensiites & fwhm, determine centers 
% no values from GRASP, no other corrections necessary
% values change with file number

phi = mf_fitter.SmoothedData.phi;
Int = mf_fitter.SmoothedData.Int;
Int_err = mf_fitter.SmoothedData.Int_err;

N = length(mf_fitter.fit_data.names);
for n=1:N
    % Fixed Values
    FWHM = mf_fitter.fit_data.fwhm(1);
    I01 = mf_fitter.fit_data.intensity1(n);
    I02 = mf_fitter.fit_data.intensity2(n);
    I03 = mf_fitter.fit_data.intensity3(n);
    
    % Guess values
    Y0 = mf_fitter.fit_data.background(n);      
    XC1 = mf_fitter.fit_data.center1(n);
    XC2 = mf_fitter.fit_data.center2(n);
    XC3 = mf_fitter.fit_data.center3(n);



    % order: y0, fwhm, i01, xc1, i02, xc2, i03, xc3
    % matlabFit(phi, int, interr, problemVarNames, problemVarValues, start, lower, upper)
    [F2 CI2 fig_h(n)] = matlabFit(phi, Int(n,:), Int_err(n,:), { 'fwhm', 'i01', 'i02', 'i03'}, {FWHM, I01, I02, I03}, [XC1, XC2, XC3, Y0], [1.5*XC1, -2, 0.5*XC3, 0.5*Y0], [0.5*XC1, 2, 1.5*XC3, 1.5*Y0])    
         
    mf_fitter.fit_data.background(n,1) = F2.y0;
    mf_fitter.fit_data.background(n,2) = CI2(2,4) - F2.y0
    %mf_fitter.fit_data.fwhm(n,1) = FWHM;
    %mf_fitter.fit_data.intensity1(n,1) = F2.i01;
    %mf_fitter.fit_data.intensity1(n,2) = CI2(2,1) - F2.i01;
    %mf_fitter.fit_data.intensity3(n,1) = F2.i03;
    %mf_fitter.fit_data.intensity3(n,2) = CI2(2,2) - F2.i03;
    mf_fitter.fit_data.center1(n,1) = F2.xc1;
    mf_fitter.fit_data.center1(n,2) = CI2(2,1) - F2.xc1;
    %mf_fitter.fit_data.center1(img_num,2) = F2;
    mf_fitter.fit_data.center2(n,1) = XC2;
    mf_fitter.fit_data.center2(n,2) = CI2(2,2) - F2.xc1;
    %mf_fitter.fit_data.center2(img_num,2) = %((fitobject.a-ci(1)) + (ci(2)-fitobject.a))/2; % + (((fitobject.b-ci(3))+(ci(4)-fitobject.b))/2)^2 );
    mf_fitter.fit_data.center3(n,1) = F2.xc3; 
    mf_fitter.fit_data.center3(n,2) = CI2(2,3) - F2.xc3;
end

    if(view3) 
    else close(fig_h) 
    end
    fig_h = [];

%% Third Fitting Cycle
% Fix centers, intensities free, determine fwhm (still keeping I02 = 0)
% But only peaks separated by more than 2xfwhm 
% no values from GRASP, no other corrections necessary
% values change with file number

phi = mf_fitter.SmoothedData.phi;
Int = mf_fitter.SmoothedData.Int;
Int_err = mf_fitter.SmoothedData.Int_err;

% for final algorithm when above conditional is required
% dxc_c = 2 * mf_fitter.fit_data.fwhm(1,1);
% N = find(abs(mf_fitter.fit_data.center1 - mf_fitter.fit_data.center3) > dxc_c);

% for extracting covariances when all data must be fit
 N = 1:length(mf_fitter.fit_data.cycles);

for n = 1:length(N)
    % Fixed Values
    XC1 = mf_fitter.fit_data.center1(N(n));
    XC2 = mf_fitter.fit_data.center2(N(n));
    XC3 = mf_fitter.fit_data.center3(N(n));

    % Guess values
    FWHM = mf_fitter.fit_data.fwhm(1);
    Y0 = mf_fitter.fit_data.background(N(n));  
    I01 = mf_fitter.fit_data.intensity1(N(n));
    I02 = mf_fitter.fit_data.intensity2(N(n));
    I03 = mf_fitter.fit_data.intensity3(N(n));
    
    
    % order: y0, fwhm, i01, xc1, i02, xc2, i03, xc3
    % matlabFit(phi, int, interr, problemVarNames, problemVarValues, start, lower, upper)
    
    % original fit with centers fixed, for final algorithm
    % [F3 CI3 fig_h(n)] = matlabFit(phi, Int(N(n),:), Int_err(N(n),:), { 'xc1', 'xc2', 'xc3'}, {XC1, XC2, XC3}, [FWHM, I01, I02, I03, Y0], [0.4*FWHM, 0.5*I01, 0.5*I02, 0.5*I03, 0.8*Y0], [1.6*FWHM, 1.25*I01, 1.25*I02, 1.25*I03, 1.2*Y0])    
      
    %uncomment for convariance -  do not fix xc1 and xc3, just give initial values
    %[F3 CI3 fig_h(n)] = matlabFit(phi, Int(N(n),:), Int_err(N(n),:),{},{},[FWHM, I01, I02, I03, XC1, XC2, XC3, Y0], [0.4*FWHM, 0.4*I01, 0.4*I02, 0.4*I03, XC1-3, XC2-0.5, XC3-3, 0.2*Y0], [1.6*FWHM, 1.6*I01, 1.6*I02, 1.6*I03, XC1+3, XC2+0.5, XC3+3, 1.8*Y0])    
     
    % uncomment for completely free fit
    [F3 CI3 fig_h(n)] = matlabFit(phi, Int(N(n),:), Int_err(N(n),:),{},{},[FWHM, FWHM, FWHM, I01, I02, I03, XC1, XC2, XC3, Y0], [0.4*FWHM, 0.4*FWHM, 0.4*FWHM, 0.4*I01, 0.4*I02, 0.4*I03, XC1-3, XC2-0.5, XC3-3, -Inf], [1.6*FWHM, 1.6*FWHM, 1.6*FWHM, 1.6*I01, 1.6*I02, 1.6*I03, XC1+3, XC2+0.5, XC3+3, Inf], 2)    
    

    %mf_fitter.fit_data.fwhm(n,1) = F3.fwhm;
    %mf_fitter.fit_data.fwhm(n,2) = CI3(2,1) - F3.fwhm;
    % comment for completely free
%     mf_fitter.fit_data.intensity1(n,1) = F3.i01;
%     mf_fitter.fit_data.intensity1(n,2) = CI3(2,2) - F3.i01;
%     mf_fitter.fit_data.intensity2(n,1) = F3.i02;
%     mf_fitter.fit_data.intensity2(n,2) = CI3(2,3) - F3.i02;
%     mf_fitter.fit_data.intensity3(n,1) = F3.i03;
%     mf_fitter.fit_data.intensity3(n,2) = CI3(2,4) - F3.i03;
    
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
    mf_fitter.fit_data.fwhm1(n,1) = F3.fwhm1;
    mf_fitter.fit_data.fwhm1(n,2) = CI3(2,1) - F3.fwhm1;
    mf_fitter.fit_data.fwhm2(n,1) = F3.fwhm2;
    mf_fitter.fit_data.fwhm2(n,2) = CI3(2,2) - F3.fwhm2;
    mf_fitter.fit_data.fwhm3(n,1) = F3.fwhm3;
    mf_fitter.fit_data.fwhm3(n,2) = CI3(2,3) - F3.fwhm3;
    
    mf_fitter.fit_data.intensity1(n,1) = F3.i01;
    mf_fitter.fit_data.intensity1(n,2) = CI3(2,4) - F3.i01;
    mf_fitter.fit_data.intensity2(n,1) = F3.i02;
    mf_fitter.fit_data.intensity2(n,2) = CI3(2,5) - F3.i02;
    mf_fitter.fit_data.intensity3(n,1) = F3.i03;
    mf_fitter.fit_data.intensity3(n,2) = CI3(2,6) - F3.i03;
    
     mf_fitter.fit_data.center1(n,1) = F3.xc1;
     mf_fitter.fit_data.center1(n,2) = CI3(2,7) - F3.xc1;
     mf_fitter.fit_data.center2(n,1) = F3.xc2;
     mf_fitter.fit_data.center2(n,2) = CI3(2,8) - F3.xc2;
     mf_fitter.fit_data.center3(n,1) = F3.xc3;
     mf_fitter.fit_data.center3(n,2) = CI3(2,9) - F3.xc3;
     mf_fitter.fit_data.background(n,1) = F3.y0;
     mf_fitter.fit_data.background(n,2) = CI3(2,10) - F3.y0;
     
     % comment when extracting covariance
     % mf_fitter.fit_data.background(n,1) = F3.y0;
     % mf_fitter.fit_data.background(n,2) = CI3(2,5) - F3.y0;
end

    if(view4) 
    else close(fig_h) 
    end
    fig_h = [];

    
%% Fourth Fitting Cycle
% Double Check centers with new fwhm
% no values from GRASP, no other corrections necessary
% values change with file number

phi = mf_fitter.SmoothedData.phi;
Int = mf_fitter.SmoothedData.Int;
Int_err = mf_fitter.SmoothedData.Int_err;

% Fixed Values
% later add code to descriminate against outliers
%nonZero = find(mf_fitter.fit_data.fwhm(:,2) > 0) % indicates that is where a fit was performed (it has an associated error)
fms = mf_fitter.fit_data.intensity2(:,1) ./ (mf_fitter.fit_data.intensity1(:,1) + mf_fitter.fit_data.intensity2(:,1) + mf_fitter.fit_data.intensity3(:,1));
goodFits = find(fms > 0.1 & fms < 0.8)
FWHM = mean(mf_fitter.fit_data.fwhm(goodFits,:))
%I02 = 0;
%XC2 = mf_fitter.fit_data.center2(1);

N = length(mf_fitter.fit_data.names);
%Nc = max(find(abs(mf_fitter.fit_data.center1 - mf_fitter.fit_data.center3) > dxc_c));

for n=1:N
    % Guess values (remember data has been smoothed & centered)  
        Y0 = mf_fitter.fit_data.background(n);  
    %if(N<Nc)
        I01 = mf_fitter.fit_data.intensity1(n,1);
        I02 = mf_fitter.fit_data.intensity2(n,1);
        I03 = mf_fitter.fit_data.intensity3(n,1);
        XC1 = mf_fitter.fit_data.center1(n,1);
        XC2 = mf_fitter.fit_data.center2(n,1);
        XC3 = mf_fitter.fit_data.center3(n,1);
%     else
%         I01 = mf_fitter.fit_data.intensity1(Nc,1);
%         I02 = mf_fitter.fit_data.intensity2(Nc,1);
%         I03 = mf_fitter.fit_data.intensity3(Nc,1);
%         XC1 = mf_fitter.fit_data.center1(Nc,1);
%         XC2 = mf_fitter.fit_data.center2(Nc,1);
%         XC3 = mf_fitter.fit_data.center3(Nc,1);
%     end
    
    % order: y0, fwhm, i01, xc1, i02, xc2, i03, xc3
    % matlabFit(phi, int, interr, problemVarNames, problemVarValues, start, lower, upper)
    [F4 CI4 fig_h(n)] = matlabFit(phi, Int(n,:), Int_err(n,:), {'fwhm','i02', 'xc2'}, {FWHM(1), I02, XC2}, [I01, XC1, I03, XC3, Y0], [0.75*I01,  0.75*I03, 1.5*XC1, 0.5*XC3, 0.75*Y0], [1.4*I01, 1.4*I03, 0.5*XC1, 1.5*XC3, 1.25*Y0])    %xc1 is neg, lower bound needs to be more neg
         
    mf_fitter.fit_data.background(n,1) = F4.y0;
    mf_fitter.fit_data.background(n,2) = CI4(2,5) - F4.y0;
    mf_fitter.fit_data.fwhm(n,1) = FWHM(1);
    mf_fitter.fit_data.fwhm(n,2) = FWHM(2);
    mf_fitter.fit_data.intensity1(n,1) = F4.i01;
    mf_fitter.fit_data.intensity1(n,2) = CI4(2,1) - F4.i01;
    mf_fitter.fit_data.intensity2(n,1) = I02;
    mf_fitter.fit_data.intensity3(n,1) = F4.i03;
    mf_fitter.fit_data.intensity3(n,2) = CI4(2,2) - F4.i03;
    mf_fitter.fit_data.center1(n,1) = F4.xc1;
    mf_fitter.fit_data.center1(n,2) = CI4(2,3) - F4.xc1;
    mf_fitter.fit_data.center2(n,1) = XC2;
    mf_fitter.fit_data.center3(n,1) = F4.xc3; 
    mf_fitter.fit_data.center3(n,2) = CI4(2,4) - F4.xc3;
    
%     save = get(mf_fitter.handles.save,'Value')
%     if(save)
%        numor = mf_fitter.fit_data.names(n); 
%        mf_fitter_save('save',['I vs Xi - ' num2str(numor)],fig_h(n)) 
%     end
    
end

    
    if(view5) 
    else close(fig_h) 
    end
    fig_h = [];

%% Fifth Fitting Cycle
% Final fit
% values change with file number

phi = mf_fitter.SmoothedData.phi;
Int = mf_fitter.SmoothedData.Int;
Int_err = mf_fitter.SmoothedData.Int_err;

% Fixed Values
% later add code to descriminate against outliers
nonZero = find(mf_fitter.fit_data.fwhm(:,2) > 0) % indicates that is where a fit was performed (it has an associated error)
FWHM = mean(mf_fitter.fit_data.fwhm(nonZero,:))

N = length(mf_fitter.fit_data.names);
%Nc = max(find(abs(mf_fitter.fit_data.center1 - mf_fitter.fit_data.center3) > dxc_c));

for r = 1:2
for n=1:N
    % Guess values (remember data has been smoothed & centered)  
        Y0 = mf_fitter.fit_data.background(n,1);  
        I01 = mf_fitter.fit_data.intensity1(n,1);
        I02 = mf_fitter.fit_data.intensity2(n,1);
        I03 = mf_fitter.fit_data.intensity3(n,1);
        XC1 = mf_fitter.fit_data.center1(n,1);
        XC2 = mf_fitter.fit_data.center2(n,1);
        XC3 = mf_fitter.fit_data.center3(n,1);

    
    % order: y0, fwhm, i01, xc1, i02, xc2, i03, xc3
    % matlabFit(phi, int, interr, problemVarNames, problemVarValues, start, lower, upper)
    [F5 CI5 fig_h(n)] = matlabFit(phi, Int(n,:), Int_err(n,:), {'fwhm'}, {FWHM(1)}, [I01, XC1, I02, XC2, I03, XC3, Y0], [0.5*I01, 0.5*I02, 0.5*I03, 1.5*XC1, XC2-1, 0.5*XC3, 0.75*Y0], [1.25*I01, 1.25*I02, 1.25*I03, 0.5*XC1, XC2+1, 1.5*XC3, 1.25*Y0])    %xc1 is neg, lower bound needs to be more neg
         
    mf_fitter.fit_data.background(n,1) = F5.y0;
    mf_fitter.fit_data.background(n,2) = CI5(2,5) - F5.y0;
    mf_fitter.fit_data.fwhm(n,1) = FWHM(1);
    mf_fitter.fit_data.fwhm(n,2) = FWHM(2);
    mf_fitter.fit_data.intensity1(n,1) = F5.i01;
    mf_fitter.fit_data.intensity1(n,2) = CI5(2,1) - F5.i01;
    mf_fitter.fit_data.intensity2(n,1) = F5.i02;
    mf_fitter.fit_data.intensity1(n,2) = CI5(2,2) - F5.i02;
    mf_fitter.fit_data.intensity3(n,1) = F5.i03;
    mf_fitter.fit_data.intensity3(n,2) = CI5(2,3) - F5.i03;
    mf_fitter.fit_data.center1(n,1) = F5.xc1;
    mf_fitter.fit_data.center1(n,2) = CI5(2,4) - F5.xc1;
    mf_fitter.fit_data.center2(n,1) = F5.xc2;
    mf_fitter.fit_data.center2(n,2) = CI5(2,5) - F5.xc2;
    mf_fitter.fit_data.center3(n,1) = F5.xc3; 
    mf_fitter.fit_data.center3(n,2) = CI5(2,6) - F5.xc3;
    
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

mf_fitter_table;
fileName = [mf_fitter.extension mf_fitter.folder '_Buzz/' 'data.txt'];
writetable(mf_fitter.handles.ExportTable, fileName);
    
%% Extra plots - check if user wants to make these
    save = get(mf_fitter.handles.save,'Value')
    
    if(save)
        phi = mf_fitter.SmoothedData.phi;
        Int = mf_fitter.SmoothedData.Int;
        cyc = mf_fitter.fit_data.cycles;

        %% Make Colormap
            cm_han = CM(phi, Int, cyc);

            %xc1 = mf_fitter.fit_data.center1(:,1)+2;
            %xc1_err = mf_fitter.fit_data.center1(:,2);
            %xc3 = mf_fitter.fit_data.center3(:,1)+2;
            %xc3_err = mf_fitter.fit_data.center3(:,2);
            %plotXC(cm_han, cyc, xc1,  xc1_err, xc3, xc3_err);

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
        
        %% Angle Decay
            set(mf_fitter.handles.plot_dropdown,'Value',6) 
            mf_GUI_callbacks( 'plot_dropdown' )
            h = gcf; 
            
            mf_fitter_save('save','AngleDecay',h)
    else
    end
    
end
% 
% %% FIT reference peaks
% notification = msgbox('Fitting reference peaks...');
% 
% mf_fitter.good_peaks.inner_num = [];
% mf_fitter.good_peaks.outers_num = [];
% 
% %locates good files' index 
% mf_fitter.good_peaks.inner_num = find(mf_fitter.fit_data.names == mf_fitter.good_peaks.inner);
% mf_fitter.good_peaks.outers_num = find(mf_fitter.fit_data.names == mf_fitter.good_peaks.outers);
% 
% %Checks to see if both successfully found
% if isempty(mf_fitter.good_peaks.inner_num)
%     disp('One-peak file could not be found');
%     return
% end
% if isempty(mf_fitter.good_peaks.outers_num)
%     disp('Two-peak file could not be found');
%     return
% end
% 
% %Fitting first (one peak) reference file
% mf_fitter.flag = 1;
% mf_fitter_callbacks('set_free',1,mf_fitter.good_peaks.inner_num); 
% mf_fitter_callbacks('fit',1,mf_fitter.good_peaks.inner_num);
% pause(0.25);
% close(mf_fitter.handles.plot_handle);
% mf_fitter.flag = 0;
% 
% if ( get(mf_fitter.handles.three_peak,'Value') )
%   mf_fitter_callbacks('set_free',3,mf_fitter.good_peaks.inner_num);
%   mf_fitter_callbacks('fit',3,mf_fitter.good_peaks.outers_num);
%   pause(0.25);
%   close(mf_fitter.handles.plot_handle);
% else
%     %Fitting second (two peak) reference file
%     mf_fitter_callbacks('set_free',2,mf_fitter.good_peaks.inner_num);
%     mf_fitter_callbacks('fit',2,mf_fitter.good_peaks.outers_num);
%     pause(0.25);
%     close(mf_fitter.handles.plot_handle);
% end
% 
% close(notification);
% 
% 
% %% AVERAGE widths, background, & MS center
% mf_fitter_callbacks('avg',0,0,'fwhm');
% mf_fitter_callbacks('avg',0,0,'background');
% mf_fitter_callbacks('avg',0,0,'center2');
% 
% %% Perform Data Smoothing
% 
% if(mf_fitter.handles.smoothing.switch  & isempty(mf_fitter.SmoothedData) )
%     disp(['Perform MRE smoothing with FWHM: ' num2str(mf_fitter.smoothing.fwhm)]);
%     for i = 1:length(mf_fitter.fit_data.names)
%         [phi, Int(i,:), Int_err(i,:)] = getData(i);
%         
%         
%         phisym = status_flags.analysis_modules.sectors.theta;
%         dTheta = status_flags.analysis_modules.sectors.delta_theta;
%         ss = mf_fitter.smoothing.step;
%         
%         %phirng = [phisym - dTheta/2; phisym + dTheta/2; ss]   % because we are centering, do not need the phisym value
%         phirng = [-dTheta/2; dTheta/2; ss]
%         [phiS, IntS(i,:), Int_errS(i,:)] = smooth(phi, Int(i,:), Int_err(i,:), i, phisym, 0, phirng, mf_fitter.smoothing.fwhm)
%         
%         fig_h = plotData(phiS, IntS(i,:), Int_errS(i,:), '\phi  - \phi_0', 'Intensity (arb. units)', ['Smoothed Data Numor - ' num2str(i)]);
%         if(view1)
%             pause(1)
%         end
%         close(fig_h)
%     end
%     
%     mf_fitter.SmoothedData.phi = phiS;
%     mf_fitter.SmoothedData.Int = IntS;
%     mf_fitter.SmoothedData.Int_err = Int_errS;
%     
%     disp('Smoothing Complete');
% elseif(mf_fitter.handles.smoothing.switch == 0)
%     disp('Smoothing Switch is turned off');
% elseif(isempty(mf_fitter.SmoothedData) == 0)
%     disp('Data has already been smoothed with these parameters')
% else
%    disp('No smoothing has been performed');
% end
% 
% 
% %%  FITTING CYCLE 1
% % FIT ALL w/ 3 peak with matlab fitter, fixed background, fwhm, & x2 to average values
% % guess values for intensity and cycles scale with how far along the transition is
% 
% phi = mf_fitter.SmoothedData.phi;
% Int = mf_fitter.SmoothedData.Int;
% Int_err = mf_fitter.SmoothedData.Int_err;
% 
% 
% N = length(mf_fitter.fit_data.names);
% 
% % Fixed Values
% Y0 = mf_fitter.fit_data.background(1);
% FWHM = mf_fitter.fit_data.fwhm(1);
% XC2 = mf_fitter.fit_data.center2(1) - status_flags.analysis_modules.sectors.theta;
% 
% % Guess values (remember data has been smoothed & centered)
% %Y0 = 2 * mf_fitter.fit_data.background(1);  % GRASP does one y0 per guassian peak
% XC1 = mf_fitter.fit_data.center1(N) - status_flags.analysis_modules.sectors.theta;
% XC3 = mf_fitter.fit_data.center3(N) - status_flags.analysis_modules.sectors.theta;
% I01 = mf_fitter.fit_data.intensity1(N);
% I02 = mf_fitter.fit_data.intensity2(1);
% I03 = mf_fitter.fit_data.intensity3(N);
% 
% for n=1:N
%     % have central peak decrease & side peaks grow
%     % have starting guess values for centers be the final positions (but give more wiggle room)
% 
%     i01 = I01*(n/N);
%     i02 = I02*(N+1-n)/(N+1);
%     i03 = I03*(n/N);
%      
%     
%     % order: y0, fwhm, i01, xc1, i02, xc2, i03, xc3
%     % matlabFit(phi, int, interr, problemVarNames, problemVarValues, start, lower, upper)
%     [F1 CI1 fig_h(n)] = matlabFit(phi, Int(n,:), Int_err(n,:), {'y0','fwhm','xc2'}, {Y0, FWHM, XC2}, [i01, i02, i03, XC1, XC3], [0.75*i01, 0.75*i02, 0.75*i03, 1.5*XC1, 0.5*XC3], [1.25*i01, 1.25*i02, 1.25*i03, 0.5*XC1, 1.5*XC3]);
%     
%              
%     mf_fitter.fit_data.background(n,1) = Y0;
%     %mf_fitter.fit_data.background(n,2) = CI1(2,3) - F1.y0;
%     mf_fitter.fit_data.fwhm(n,1) = FWHM;
%     mf_fitter.fit_data.intensity1(n,1) = F1.i01;
%     mf_fitter.fit_data.intensity1(n,2) = CI1(2,1) - F1.i01;
%     mf_fitter.fit_data.intensity2(n,1) = F1.i02;
%     mf_fitter.fit_data.intensity2(n,2) = CI1(2,2) - F1.i02;
%     mf_fitter.fit_data.intensity3(n,1) = F1.i03;
%     mf_fitter.fit_data.intensity3(n,2) = CI1(2,3) - F1.i03;
%     mf_fitter.fit_data.center1(n,1) = F1.xc1;
%     mf_fitter.fit_data.center1(n,2) = CI1(2,4) - F1.xc1;
%     mf_fitter.fit_data.center2(n,1) = XC2;
%     %mf_fitter.fit_data.center2(img_num,2) = %((fitobject.a-ci(1)) + (ci(2)-fitobject.a))/2; % + (((fitobject.b-ci(3))+(ci(4)-fitobject.b))/2)^2 );
%     mf_fitter.fit_data.center3(n,1) = F1.xc3; 
%     mf_fitter.fit_data.center3(n,2) = CI1(2,5) - F1.xc3;
%     
% end
% 
%     if(view2) 
%     else close(fig_h) 
%     end
%     fig_h = [];
%     
%     %  FITTING CYCLE 2
% Prompts user to re-examine bad peaks
% Select the last peak you were confident in, attempt to refit once, then set what value to fix it to.
% 
% name = ['Fitting Cycle: 2 of ', num2str(Nc), '...']
% notification = waitbar(0, name, 'Name', 'Please wait','units', 'Normalized', 'Position', wait_position);
% 
% mf_fitter.refit.continue = 0;
% refit_GUI_window('initialize');
% while(mf_fitter.refit.continue ~= 1)
%    pause(2)
% end
% 
% close(notification);
% close(mf_fitter.handles.table);
% close(mf_fitter.handles.cplot);
% 
% mf_fitter_table;
% 
% mf_fitter_callbacks('center_store');
% mf_fitter_callbacks('center_plot')
% 
%  FITTING CYCLE 3
% FIT ALL w/ 3 peak with grasp fitter, fixed fwhm & x2

