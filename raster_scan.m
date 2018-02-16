function [ ] = raster_scan()

global mf_fitter;
global status_flags;
global grasp_handles;


mf_fitter_callbacks('initialize');
  % case 'set_params'
        % case run before every fit
        % called by all of the various set_ cases
        % these are parameters strictly determined by num_peaks
        num_peaks = 3;
        
        status_flags.fitter.number1d = num_peaks;
        status_flags.fitter.function_info_1d.no_parameters = 4;
        status_flags.fitter.function_info_1d.variable_names = repmat({'y0', 'i0', 'xc', 'fwhm'},1,num_peaks);
        status_flags.fitter.function_info_1d.long_names = repmat({'Background', 'Integrated Intensity', 'Center', 'FWHM'},1,num_peaks);
        
        % clear out old data/options
        status_flags_names = {'group','fix','values','err_values'}
        for i=1:length(status_flags_names)
            variable = char(status_flags_names(i));
           status_flags.fitter.function_info_1d.(variable) = []; 
        end
        
     
        status_flags.fitter.function_info_1d.group = repmat([1,0,0,1],1,num_peaks);
        status_flags.fitter.function_info_1d.err_values = repmat([0,0,0,0],1,num_peaks);
        status_flags.fitter.function_info_1d.values = repmat([0],1,4*num_peaks);
        status_flags.fitter.function_info_1d.fix = [0,0,1,1,1,0,1,1,1,0,1,1];
        
        status_flags.fitter.function_info_1d.values(1) = 0.04;
        status_flags.fitter.function_info_1d.values(2) = 5;
        status_flags.fitter.function_info_1d.values(3) = 239.2;
        status_flags.fitter.function_info_1d.values(4) = 6.37204; 
        status_flags.fitter.function_info_1d.values(6) = 5;
        status_flags.fitter.function_info_1d.values(7) = 246.38;
        status_flags.fitter.function_info_1d.values(10) = 5;
        status_flags.fitter.function_info_1d.values(11) = 253.93;
        
       % case 'fit'
        % load desired depth file and update main grasp GUI
        % the +1 is necessary as file 1 represents a sum
        
        
     for img_num = 1:mf_fitter.depth
        status_flags.selector.fd = img_num+1;
        grasp_update;

        % plot I vs xi and create the current figure, store handle
        radial_average_callbacks('averaging_control','azimuthal');
        mf_fitter.handles.plot_handle = gcf;
        grasp_update;

          grasp_plot_fit_callbacks_2('auto_guess','off');
          
        % perform the fit
        set(grasp_handles.window_modules.curve_fit1d.curve_number,'value',1)
        try
            grasp_plot_fit_callbacks_2('fit_it');
            pause(2)
            % store the temporary fit parameters to a permanent structure
             mf_fitter_callbacks('data_storage',num_peaks,img_num)
             close(gcf);
        catch
            mf_fitter.fit_data.intensity1 = 0;
            mf_fitter.fit_data.intensity2 = 0;
            mf_fitter.fit_data.intensity3 = 0;
        end
     end
     
     %% Image
     % 8 x 9 pixels

% create pixel matrix to fill
% p(height, width, rgb)
p = zeros(9,8,3) 
m = zeros(9,8,3) 
g = zeros(9,8,3) 
test = zeros(9,8,3)
IMS = zeros(9,8)
IGS1 = zeros(9,8) 
IGS2 = zeros(9,8)
 
count = 3;

   for y = 1:9
      for x = 1:8
         for c = 1:3
            
            if( (y==1) || (y==2 && x<5) || (y==3 && (x<4 || x>7)) || (y==4 && (x<3 || x>7)) || (y==5 && (x<3 || x>6)) || (y==6 && (x<3 || x>7)) || (y==7 && (x<3 || x>7)) || (y==8 && (x<3 || x>6)) || (y==9 && (x<4 || x>4)) )
               p(y,x,c) = 0;
               m(y,x,c) = 0;
               g(y,x,c) = 0;
               test(y,x,c) = 0;
               IMS(y,x,c) = 0;
               IGS1(y,x,c) = 0;
               IGS2(y,x,c) = 0;
            else
                    count = count + 1;
                    index = count/3;
                    element_num = floor(index)
                   
                   I1 = mf_fitter.fit_data.intensity1(element_num);
                   I2 = mf_fitter.fit_data.intensity2(element_num);
                   I3 = mf_fitter.fit_data.intensity3(element_num);


                    if(I1 < 0)
                        I1 = 0;
                    end
                    if(I2 < 0)
                        I2 = 0;
                    end
                    if(I3 < 0)
                        I3 = 0;
                    end

                  %  if( (I1 + I2 + I3) < 1)
                  %      I1 = 0;
                  %      I2 = 0;
                  %      I3 = 0;
                  %  end

                   fp1 = I1 / (I1 + I2 + I3);
                   fp2 = I2 / (I1 + I2 + I3);
                   fp3 = I3 / (I1 + I2 + I3);

                   fp = [fp1, fp3, fp2];
                   ms = [fp2, 1, 1]
                   gs = [1, 1, fp1]
                   test = [1, fp3, 1]
                   
               p(y, x, c) = fp(c);
               m(y, x, c) = ms(c);
               g(y, x, c) = gs(c);
               test(y,x,c) = test(c);
               IMS(y,x,c) = I2;
               IGS1(y,x,c) = I1;
               IGS2(y,x,c) = I3;
            end
            
         end
      end
   end
    
   IMS
   IGS1
   IGS2
   figure(12)
    image(p)
    figure(13)
    image(m)
    figure(14)
    image(g)
    figure(15)
    image(test)
    %figure(11)
    %image(IMS / max(IMS))


end