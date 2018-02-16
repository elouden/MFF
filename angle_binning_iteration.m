function [ ] = angle_binning_iteration( abifun )
    
    global mf_fitter;
    global grasp_handles;
    global grasp_env;

    
    switch abifun
        
        case 'open'
            %% Draw Window
            fig_position = [250, 250, 150, 250];
            mf_fitter.handles.ab_window = figure(....
                'units','pixels',....
                'Position',fig_position,....
                'Name','Angle Binning Iteration' ,....
                'NumberTitle', 'off',....
                'Tag','AB_window',...
                'color',grasp_env.background_color,....
                'menubar','none',....
                'resize','off',....
                'closerequestfcn','sector_callbacks2(''close'');closereq');

            handle = mf_fitter.handles.ab_window;

            % initialize values
            mf_fitter.abi.initial = 1.0;
            mf_fitter.abi.final = 2.2;
            mf_fitter.abi.stepsize = 0.2;
            current_ab = 0;


            % Initial Angle Binning
            uicontrol(handle,'units','normalized','Position',[0.15 0.85 0.7 0.1],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','center','Style','text','String','Initial Angle Bin:','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
            mf_fitter.handles.initial_ab = uicontrol(handle,'units','normalized','Position',[0.35 0.8 0.3 0.08],'tooltip','Sets the starting value for the angle binning.','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','edit','String',num2str(mf_fitter.abi.initial),'HorizontalAlignment','left','Tag','angular_bin_Start','Visible','on');

            % Final Angle Binning
            uicontrol(handle,'units','normalized','Position',[0.15 0.65 0.7 0.1],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','center','Style','text','String','Final Angle Bin:','BackgroundColor',grasp_env.background_color,'ForegroundColor',[1 1 1]);
            mf_fitter.handles.final_ab = uicontrol(handle,'units','normalized','Position',[0.35 0.6 0.3 0.08],'tooltip','Sets the final value for the angle binning.','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','edit','String',num2str(mf_fitter.abi.final),'HorizontalAlignment','left','Tag','angular_bin_End','Visible','on');

            % Angle Binning Step Size
            uicontrol(handle,'units','normalized','Position',[0.15 0.45 0.7 0.1],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','center','Style','text','String','Step Size:','BackgroundColor', grasp_env.background_color, 'ForegroundColor', [1 1 1]);
            mf_fitter.handles.step_size = uicontrol(handle,'units','normalized','Position',[0.35 0.4 0.3 0.08],'tooltip','Sets the step size for the angle binning iteration.','FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'Style','edit','String',num2str(mf_fitter.abi.stepsize),'HorizontalAlignment','left','Tag','angular_bin_StepSize','Visible','on');

            % Clear Button
             uicontrol(handle,'units','normalized','Position',[0.1 0.125 0.35 0.2],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','center','Style','pushbutton','String','Clear AB Values','BackgroundColor',[1 1 1],'Callback','angle_binning_iteration(''clear'')');
            
            % Go Button
            uicontrol(handle,'units','normalized','Position',[0.5 0.125 0.35 0.2],'FontName',grasp_env.font,'FontSize',grasp_env.fontsize,'HorizontalAlignment','center','Style','pushbutton','String','GO!','BackgroundColor',[1 1 1],'Callback','angle_binning_iteration(''run'')');

        case 'clear'
            mf_fitter.abi.ab = [];
            mf_fitter.abi.chi2 = [];
        
        case 'run'
            current_ab = str2num(get(mf_fitter.handles.initial_ab,'String'));
            if( isfield(mf_fitter.abi, 'ab') )
                i = length(mf_fitter.abi.ab) + 1;
                if(i == 0) i=1; end
            else
                i = 1;
            end
            
            while( current_ab <= str2num(get(mf_fitter.handles.final_ab,'String')) )
                set(grasp_handles.window_modules.radial_average.azimuth_bin,'String',num2str(current_ab));
                radial_average_callbacks2('azimuth_bin');
                mf_GUI_callbacks('apply_sector');
                mf_fitter.folder = ['AngleBinning_' num2str(current_ab) '/'];
                dir = [mf_fitter.extension mf_fitter.folder];
                mkdir(dir);
                mf_fitter.data.intensity = [];
                mf_fitter.data.intensity_err = [];
                mf_fitter.data.angle = [];

                mf_fitter_main();
               
                mf_fitter_callbacks('fraction');
                    mf_fitter.file = ['FMS_' num2str(current_ab)];
                    angle_binning_iteration('save');
                    close(gcf);
                mf_fitter_callbacks('colormapLin');
                    mf_fitter.file = ['ColormapLin_' num2str(current_ab)];
                    angle_binning_iteration('save');
                    close(gcf);
                mf_fitter_callbacks('colormapLog');
                    mf_fitter.file = ['ColormapLog_' num2str(current_ab)];
                    angle_binning_iteration('save');
                    close(gcf);
                mf_fitter_callbacks('peak_decay');
                    mf_fitter.file = ['PeakDecay_' num2str(current_ab)];
                    angle_binning_iteration('save');
                    close(gcf);
               % mf_fitter_callbacks('center_plot');
              %      mf_fitter.file = ['CenterPlot_' num2str(current_ab)];
               %     angle_binning_iteration('save');  
                %    close(gcf);
                    
                mf_fitter.file = ['Data_' num2str(current_ab)];
                angle_binning_iteration('export');
                
                mf_fitter.abi.ab(i) = current_ab;
                mf_fitter.abi.chi2(i) = mean(mf_fitter.fit_data.chi2)
                i = i+1;
                current_ab = current_ab + str2num(get(mf_fitter.handles.step_size,'String'));
                
            end
         
            ab = mf_fitter.abi.ab;
            chi2 = mf_fitter.abi.chi2;
            
            ab_div = chi2 ./ ab;
            sqrt_div = chi2 ./ sqrt(ab);
            h = figure
            hold on
            scatter(ab, chi2, 'ro');
            scatter(ab, ab_div);
            scatter(ab, sqrt_div);
            title('\fontsize{24}\color{black} Angle Binning Optimization')
            xlabel('\fontsize{16}\color{black} Angle Binning')
            ylabel('\fontsize{16}\color{black} Error')
            set(gca,'Color',[1,1,1])
            set(gca,'XColor',[0,0,0])
            set(gca,'YColor',[0,0,0])
            legend(gca, {'Mean Error', 'Divided by AB', ' Divided by Sqrt of AB'}, 'Color', [1,1,1], 'Fontcolor', [0,0,0])
            hold off
   
            beep
            beep
            
            
        case 'export'    
            %Name export file
           % current_ab = str2num(get(mf_fitter.handles.initial_ab,'String'));
           current_ab = 1.29
            mf_fitter.folder = ['AngleBinning_' num2str(current_ab) '/'];
            mf_fitter.file = ['Data_' num2str(current_ab)];
            fileName= [mf_fitter.extension mf_fitter.folder mf_fitter.file '.txt'];

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
            
        case 'save'
            current_ab = str2num(get(grasp_handles.window_modules.radial_average.azimuth_bin,'String'));
            mf_fitter.folder = ['AngleBinning_' num2str(current_ab) '/'];
            %mf_fitter.file = ['Data_' num2str(current_ab)];
            dir = [mf_fitter.extension mf_fitter.folder];
            mkdir(dir);
            figfile = [mf_fitter.extension mf_fitter.folder mf_fitter.file '.fig'];
            jpgfile = [mf_fitter.extension mf_fitter.folder mf_fitter.file '.jpg'];
            saveas(mf_fitter.handles.fig,figfile);
            saveas(mf_fitter.handles.fig,jpgfile);
            
    end
    
end

