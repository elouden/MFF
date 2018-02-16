function [] = mf_export(opt)
%Export fit data to txt

% MFF - Allan - 8/8/14

global mf_fitter;


choice = questdlg('Would you like to export these parameters as a .txt file?', ...
	'Export Parameters?','Yes','No', 'No');

switch choice
    
    case 'Yes'
        %prompt user for export file and add to path
        [destination, pathname] = uiputfile ( '*.txt', 'Select an export file');
        addpath (pathname);

        %Name export file
        fileName= [pathname destination];

        %Creating Table parameters
        Numors = mf_fitter.fit_data.names;
        if(isempty(mf_fitter.fit_data.cycles))
            Cycle = zeros(mf_fitter.depth,1); 
        else
            Cycle = mf_fitter.fit_data.cycles';
        end
        Background = mf_fitter.fit_data.background;
        FWHM = mf_fitter.fit_data.fwhm;
        Break = zeros(mf_fitter.depth,1);
        I1 = mf_fitter.fit_data.intensity1;
        X1 = mf_fitter.fit_data.center1;
        I2 = mf_fitter.fit_data.intensity2;
        X2 = mf_fitter.fit_data.center2;
        I3 = mf_fitter.fit_data.intensity3;
        X3 = mf_fitter.fit_data.center3;
        Chi2 = mf_fitter.fit_data.chi2;
        if(isempty(mf_fitter.fit_data.fms)== 0)
            FMS = mf_fitter.fit_data.fms';
            FMS2 = mf_fitter.error_prop.sig_frac;
            FGS = mf_fitter.fit_data.fgs';
            ExportTable = table(Numors, Cycle, Chi2, Break, Background, FWHM, Break, I1, X1, Break, I2, X2, Break, I3, X3, FMS, FMS2, FGS);
        else
            %Creating table
            ExportTable = table(Numors, Cycle, Chi2, Break, Background, FWHM, Break, I1, X1, Break, I2, X2, Break, I3, X3);
        end
        
        %exporting table
        writetable(ExportTable, fileName);
end

end

