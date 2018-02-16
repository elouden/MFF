function [ output_args ] = fit_anlaysis( field )
%This function analyzes the quality of the fits by generating plots
%exploring the chi^2
%field should be a string 

global D11

bin = {'anglebin1_0','anglebin1_2','anglebin1_4','anglebin1_6','anglebin1_8','anglebin2_0'};
name = {'Angle Bin 1.0','Angle Bin 1.2', 'Angle Bin 1.4', 'Angle Bin 1.6', 'Angle Bin 1.8', 'Angle Bin 2.0'};

  axes1 = axes('Parent',figure,...
            'ZColor',[0 0 0],'YColor',[0 0 0],'XColor',[0 0 0],...
            'FontName','Times New Roman','FontColor',[0,0,0],...
            'Color',[0 0 0]); 
         hold(axes1,'all');
  
         
         
for i=1:6
    bin2 = bin(i);
    bin2 = char(bin2);
    Numors = D11.(field).(bin2).names;
    error = D11.(field).(bin2).chi2;
    legend_name = name(i);
    legend_name = char(legend_name);
    scatter(Numors,error,'DisplayName',legend_name)
end

title('\fontsize{18}\color{black} Deterimination of Angle Binning for 125 G: Distribution')
xlabel('\fontsize{14}\color{black} Numor')
ylabel('\fontsize{14}\color{black} Chi^2')
legend(gca,'show')

hold off;


 axes2 = axes('Parent',figure,...
            'ZColor',[0 0 0],'YColor',[0 0 0],'XColor',[0 0 0],...
            'Color',[1 1 1]); 
         hold(axes2,'all');

for i=1:6
    bin2 = bin(i);
    bin2 = char(bin2);
    avg_chi2(i) = mean(D11.(field).(bin2).chi2);
end
bin = 1.0:0.2:2.0;
scatter(bin,avg_chi2,'bo')

title('\fontsize{18}\color{black} Deterimination of Angle Binning for 125 G:  Avg Chi^2')
xlabel('\fontsize{14}\color{black} Angle Binning')
ylabel('\fontsize{14}\color{black} Average Chi^2')


end

