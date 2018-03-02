function [] = plot_template( scale )
% 10/31/2017 - E R Louden

% This function sets all of the preferred options for formatting plots for publications
% It should be called when the current figure is essentially complete and still the "active" figure
% Note:  it turns "hold" off, so if you have a plot command immediately following plot_template, 
% it will overwrite your figure!

% The scale changes how large the figure is:
%   for publication: scale = 1
%   for visualization during analysis: scale = 2 (it is easier to see what you are doing!!)

% The current preferences are:
%   Aspect Ratio/Size - 5.7 cm x 4.45 cm
%   Font Name - Arial
%   Font Size - 8 pt
%   The axes should be a closed box


% If scale unspecified, use 2
if(nargin<1)
    scale = 2;
end


% Figure Size
fig_size = [5.7 4.45];
fig_position = [1.25 1.25];
position = scale*[fig_position fig_size];
hold on


% Figure Properties
set(gca,...
    'Color',[1 1 1],...
    'FontName','Arial','FontSize',scale*8,...
    'Xcolor',[0 0 0],'Ycolor',[0 0 0],...   
    'Units','centimeters','Position', position,...
    'XMinorTick', 'on', 'YMinorTick', 'on');

set(gca,'Box','on')


% Legend
h = findobj('tag','legend');
if(isempty(h))
else
    set(h,'Color',[1 1 1],'FontName','Arial','FontSize',scale*8,'EdgeColor', [1 1 1], 'Textcolor', [0 0 0])
end


% Format Axes
h_title = get(gca,'title');
set(h_title,'FontName','Arial','FontSize',scale*10,'Color',[0 0 0]);

h_xlabel = get(gca,'XLabel');
set(h_xlabel,'FontName','Arial','FontSize',scale*8,'Color',[0 0 0]);

h_ylabel = get(gca,'YLabel');
set(h_ylabel,'FontName','Arial','FontSize',scale*8,'Color',[0 0 0]);


hold off

% Remove internal variables
clear h_title h_xlabel h_ylabel fig_size fig_position position f
end