function[cm_han, ps_han, int_han, int_tot_han, ad_han] = mf_fitter_plotMaker()
% makes plots
% handles returned:
%   cm_han      -   Colormap
%   ps_han      -   Peak Separation
%   int_han     -   Intensity of Peaks (individually)
%   int_tot_han -   Intensity of Peaks (summed)
%   ad_han      -   Angle Decay

        %% Initialize
        global mf_fitter
        
        phi = mf_fitter.SmoothedData.phi;
        Int = mf_fitter.SmoothedData.Int;
        cyc = mf_fitter.fit_data.cycles;
        
        
        %% Make Colormap
        try
            cm_han = CM(phi, Int, cyc);
        catch
            cm_han = 0;
        end
            
        
        %% Peak Separation
        try
            mf_fitter_callbacks('center_separation');
            y  = abs(mf_fitter.fit_data.center3(:,1) - mf_fitter.fit_data.center1(:,1));
            y_err = sqrt((mf_fitter.fit_data.center3(:,2)).^2+(mf_fitter.fit_data.center1(:,2)).^2);
            ps_han = figure;
            plotData(cyc, y, y_err, 'Applied || AC Cycles', '\Delta \phi (degrees)', ['Peak Separation - ' mf_fitter.folder])
            set(gca,'xScale','log')
        catch
            ps_han = 0;
        end
        
        
        %% Intensity
        try
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
        catch
            int_han = 0;    
        end
        
        try
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
        catch
            int_tot_han = 0;
        end

        
        %% Angle Decay
        try
            set(mf_fitter.handles.plot_dropdown,'Value',6) 
            mf_GUI_callbacks( 'plot_dropdown' )
            ad_han = gcf;
        catch
            ad_han = 0;
        end
