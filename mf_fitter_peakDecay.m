function [ ] = mf_fitter_peakDecay( xlow, xhigh)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

global mf_fitter
 disp('Peak Decay Plot')
       % mf_fitter.handles.fig = figure('PaperSize',[8.3 11.7],...
        %    'Color',[0.80 0.80 0.80]);

        h = figure('PaperSize',[8.3 11.7],...
            'Color',[0.80 0.80 0.80]);
        
        hold on
%         mf_fitter_callbacks('avg',0,0,'center1');
%         mf_fitter_callbacks('avg',0,0,'center2');
%         mf_fitter_callbacks('avg',0,0,'center3');
%         
%         x1 = mf_fitter.averages.center1;
%         x2 = mf_fitter.averages.center2;
%         x3 = mf_fitter.averages.center3;
%         xlow = x2-20;
%         xhigh = x2+20;

        
        
        scale = 1.25*max(max(mf_fitter.SmoothedData.Int));
        
         subplot2 = subplot(1,1,1,'Parent',h,...
            'XLim', [xlow, xhigh],...
            'XDir','reverse',...
            'YLim', [-(scale*mf_fitter.depth), 10],...
            'Fontname','Timesnewroman',...
            'FontSize',12);
        
        set(gca,'Color',[1, 1, 1], 'XColor',[0, 0, 0], 'YColor',[0, 0, 0]);
        
        phi = mf_fitter.SmoothedData.phi
        for i = 1:mf_fitter.depth
           intensity2 = mf_fitter.SmoothedData.Int(i,:) - scale*i;
           errorbar(phi,intensity2,mf_fitter.SmoothedData.Int_err(i,:),...
                              'bo','LineWidth',1,...
                              'MarkerEdgeColor','k',...
                              'MarkerFaceColor','k',...
                              'MarkerSize',5);
                          
           y0 = mf_fitter.fit_data.background(i)-scale*i;
           fwhm = mf_fitter.fit_data.fwhm(i);
           i01 = mf_fitter.fit_data.intensity1(i);
           i02 = mf_fitter.fit_data.intensity2(i)
           i03 = mf_fitter.fit_data.intensity3(i);
           x1 = mf_fitter.fit_data.center1(i);
           x2 = mf_fitter.fit_data.center2(i);
           x3 = mf_fitter.fit_data.center3(i);
   
           fplot(@(x) y0 + (i01/(fwhm * sqrt(pi/2) / sqrt(log(4))))*exp(-2*((x-x1)^2/(fwhm^2/log(4)))) + (i02/(fwhm * sqrt(pi/2) / sqrt(log(4))))*exp(-2*((x-x2)^2/(fwhm^2/log(4)))) + (i03/(fwhm * sqrt(pi/2) / sqrt(log(4))))*exp(-2*((x-x3)^2/(fwhm^2/log(4))) ),...
                        [xlow xhigh], 2, 'r'); 
        end
        titlename = [mf_fitter.folder ' Peak Decay ']% get(grasp_handles.window_modules.radial_average.azimuth_bin,'String')];
        title(titlename,'FontSize',16,'Fontname','Arial','Color','black');
        xlabel('\phi - \phi_0 (degrees)','FontSize',12,'Fontname','Arial','Color','black');
        ylabel('Relative Intensity (arb. units)','FontSize',12,'Fontname','Arial','Color','black');
       
        set(gca,'View',[-90,90]);
        set(gca,'YLim',[-scale*(mf_fitter.depth+1),scale]);

        hold off
        
        Position = [1.5 1.5 (5.7/2.5)*mf_fitter.depth 1.5*4.45]
        set(gca,'Units','centimeters','Position',Position)
end

