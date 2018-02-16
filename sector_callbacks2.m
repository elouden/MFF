function smask = sector_callbacks2(to_do,options)

%Same as sector_callbacks, but with messy parts commented out

% v 7.4
% 7/22/14 MFF Allan

if nargin<2; options =[]; end
if nargin<1; to_do =''; end


global grasp_handles
global status_flags
global displayimage
global inst_params
global grasp_env

smask = []; %dummy value for when not used


switch to_do
    
    case 'inner_radius'
        temp = str2num(get(grasp_handles.window_modules.sector.inner_radius,'string'));
        if not(isempty(temp)); status_flags.analysis_modules.sectors.inner_radius = temp; end
        
    case 'outer_radius'
        temp = str2num(get(grasp_handles.window_modules.sector.outer_radius,'string'));
        if not(isempty(temp)); status_flags.analysis_modules.sectors.outer_radius = temp; end
        
    case 'theta'
        temp = str2num(get(grasp_handles.window_modules.sector.theta,'string'));
        if not(isempty(temp));
            temp = check_angle(temp);
            status_flags.analysis_modules.sectors.theta = temp; 
        end
        
    case 'delta_theta'
        temp = str2num(get(grasp_handles.window_modules.sector.delta_theta,'string'));
        if not(isempty(temp));
            status_flags.analysis_modules.sectors.delta_theta = temp; end
        
    case 'angle_plus'
        status_flags.analysis_modules.sectors.theta = status_flags.analysis_modules.sectors.theta + options;
        status_flags.analysis_modules.sectors.theta = check_angle(status_flags.analysis_modules.sectors.theta);
        
    case 'mirror_sectors'
        status_flags.analysis_modules.sectors.mirror_sectors = get(grasp_handles.window_modules.sector.mirror_sectors,'value');
        
    case 'color'
        color_string = get(grasp_handles.window_modules.sector.sector_color,'string');
        position = get(grasp_handles.window_modules.sector.sector_color,'value');
        color = color_string{position};
        status_flags.analysis_modules.sectors.sector_color = color;
        
    case 'anisotropy'
        temp = str2num(get(grasp_handles.window_modules.sector.anisotropy,'string'));
        if not(isempty(temp)); status_flags.analysis_modules.sectors.anisotropy = temp; end
        
    case 'anisotropy_angle'
        temp = str2num(get(grasp_handles.window_modules.sector.anisotropy_angle,'string'));
        if not(isempty(temp)); status_flags.analysis_modules.sectors.anisotropy_angle = temp; end
        
        
    case 'close'
        %Delete old sectors
        temp = find(ishandle(grasp_handles.window_modules.sector.sketch_handles));
        delete(grasp_handles.window_modules.sector.sketch_handles(temp));
        grasp_handles.window_modules.sector.sketch_handles = [];
        grasp_handles.window_modules.sector.window = [];
        return
        
    case 'radial_average'
        status_flags.analysis_modules.radial_average.sector_mask_chk = 1;
        radial_average_window;
        return
        
    case 'build_sector_mask'
                
        for det = 1:inst_params.detectors
            
            %Find current detector distance for particular detector pannel
            det_current = displayimage.(['params' num2str(det)])(inst_params.vectors.det); %Default unless otherwise
            if strcmp(status_flags.q.det,'detcalc');
                if isfield(inst_params.vectors,'detcalc')
                    if not(isempty(inst_params.vectors.detcalc));
                        det_current = displayimage.(['params' num2str(det)])(inst_params.vectors.detcalc);
                    end
                end
            end
            if isfield(inst_params.vectors,'det_pannel');
                det_current = displayimage.(['params' num2str(det)])(inst_params.vectors.det_pannel);
            end
            
            %Keep a memory of the main detector(1) Det
            if det ==1; det1det = det_current; end
                    
                    
            %Prepare Sector Coordinates
            if isempty(options) %i.e. just coming from the sectors window
                inner_radius = status_flags.analysis_modules.sectors.inner_radius * inst_params.detector1.pixel_size(1)/1000;
                outer_radius = status_flags.analysis_modules.sectors.outer_radius * inst_params.detector1.pixel_size(1)/1000;
                theta_set = status_flags.analysis_modules.sectors.theta;
                delta_theta = status_flags.analysis_modules.sectors.delta_theta;
                mirrors = status_flags.analysis_modules.sectors.mirror_sectors;
            else %i.e specific coordinates coming from the Sector BOX window
                inner_radius = options(1) * inst_params.detector1.pixel_size(1)/1000;
                outer_radius = options(2) * inst_params.detector1.pixel_size(1)/1000;
                theta_set = options(3);
                delta_theta = options(4);
                mirrors = options(5);
            end
            %Check for r1 and r2 the wrong way round
            if inner_radius > outer_radius;
                temp = inner_radius; inner_radius = outer_radius; outer_radius = temp;
            end
            
            %For conical radius scaling of multiple detectors
            eff_outer_radius = outer_radius * det_current / det1det;
            eff_inner_radius = inner_radius * det_current / det1det;
            
            %Empty sector masks for all detectors and partial masks, smask1, smask2 etc.
            smask.(['det' num2str(det)]) = zeros(inst_params.(['detector' num2str(det)]).pixels(2),inst_params.(['detector' num2str(det)]).pixels(1));
            smask1 = zeros(inst_params.(['detector' num2str(det)]).pixels(2),inst_params.(['detector' num2str(det)]).pixels(1));
            smask2 = zeros(inst_params.(['detector' num2str(det)]).pixels(2),inst_params.(['detector' num2str(det)]).pixels(1));
            smask3 = zeros(inst_params.(['detector' num2str(det)]).pixels(2),inst_params.(['detector' num2str(det)]).pixels(1));
            smask4 = zeros(inst_params.(['detector' num2str(det)]).pixels(2),inst_params.(['detector' num2str(det)]).pixels(1));
            
            %Define dectors by inner and outer radius.
            %Note:  Sectors are currently defined in pixel radii.
            %qmatrix(14),(15),(16) define x, y, mod radii in distance (m).
            smask1(displayimage.(['qmatrix' num2str(det)])(:,:,16)>= eff_inner_radius & displayimage.(['qmatrix' num2str(det)])(:,:,16) <= eff_outer_radius) = 1;
            
            %Now define sectors by angle (for multiple mirrors)
            angle_array = displayimage.(['qmatrix' num2str(det)])(:,:,6); %pixel angle around beam centre
            
            for n = 1:mirrors
                theta = theta_set + (n-1)*(360/mirrors);
                theta = check_angle(theta);
                %Angle
                smask2(angle_array >= (theta - (delta_theta/2)) & angle_array <= (theta + (delta_theta/2))) = 1;
                %Correct for sector straddling 0 or 360
                if (theta - delta_theta/2) < 0; smask3(angle_array > ((theta - delta_theta /2) + 360)) = 1; end
                if (theta + delta_theta/2) > 360; smask4(angle_array < ((theta + delta_theta /2) - 360)) = 1; end
                
                %Final logical mask
                smask.(['det' num2str(det)]) = or(smask.(['det' num2str(det)]),(smask1.*(smask2+smask3+smask4)));
                smask.(['det' num2str(det)]) = +smask.(['det' num2str(det)]); %Converts from logical to integer
            end
            
            %Include Current Mask conditions in the box
            smask.(['det' num2str(det)]) = smask.(['det' num2str(det)]).*displayimage.(['mask' num2str(det)]);
            
            
            
            
            % figure
            % pcolor(smask.(['det' num2str(det)]))
        end
        return
end

%make a check for inner and outer radii the wrong way round
if status_flags.analysis_modules.sectors.outer_radius < status_flags.analysis_modules.sectors.inner_radius;
    temp = status_flags.analysis_modules.sectors.inner_radius;
    status_flags.analysis_modules.sectors.inner_radius = status_flags.analysis_modules.sectors.outer_radius;
    status_flags.analysis_modules.sectors.outer_radius = temp;
end


%Update the window items
% set(grasp_handles.window_modules.sector.window.inner_radius,'string',num2str(status_flags.analysis_modules.sectors.inner_radius));
% 
% set(grasp_handles.window_modules.sector.window.outer_radius,'string',num2str(status_flags.analysis_modules.sectors.outer_radius));
% 
% set(grasp_handles.window_modules.sector.window.theta,'string',num2str(status_flags.analysis_modules.sectors.theta));
% 
% set(grasp_handles.window_modules.sector.window.delta_theta,'string',num2str(status_flags.analysis_modules.sectors.delta_theta));
% 
% set(grasp_handles.window_modules.sector.anisotropy,'string',num2str(status_flags.analysis_modules.sectors.anisotropy));
% 
% set(grasp_handles.window_modules.sector.anisotropy_angle,'string',num2str(status_flags.analysis_modules.sectors.anisotropy_angle));

%Delete old sectors
temp = find(ishandle(grasp_handles.window_modules.sector.sketch_handles));
delete(grasp_handles.window_modules.sector.sketch_handles(temp));
grasp_handles.window_modules.sector.sketch_handles = [];

%Draw new sectors
for det = 1:inst_params.detectors
    cm = displayimage.cm.(['det' num2str(det)]);
    params = displayimage.(['params' num2str(det)]);
    
    %Find current detector distance for particular detector pannel
    det_current = params(inst_params.vectors.det); %Default unless otherwise
    if strcmp(status_flags.q.det,'detcalc');
        if isfield(inst_params.vectors,'detcalc')
            if not(isempty(inst_params.vectors.detcalc));
                det_current = params(inst_params.vectors.detcalc);
            end
        end
    end
    if isfield(inst_params.vectors,'det_pannel');
        det_current = params(inst_params.vectors.det_pannel);
    end
    
    %Keep a memory of the main detector(1) Det
    if det ==1; det1det = det_current; end
    
    %For conical radius scaling of multiple detectors
    warning off
    eff_outer_radius = status_flags.analysis_modules.sectors.outer_radius * det_current / det1det;
    eff_inner_radius = status_flags.analysis_modules.sectors.inner_radius * det_current / det1det;
    warning on 

    
    %Pixel distances from Beam Centre
    if isfield(inst_params.vectors,'ox')
        cx_eff = cm.cm_pixels(1) - ((params(inst_params.vectors.ox) - cm.cm_translation(1))/inst_params.detector1.pixel_size(1));
    else
        cx_eff = cm.cm_pixels(1);
    end
    
    if isfield(inst_params.vectors,'oy')
        cy_eff = cm.cm_pixels(2) - ((params(inst_params.vectors.oy) - cm.cm_translation(2))/inst_params.detector1.pixel_size(1));
    else
        cy_eff = cm.cm_pixels(2);
    end

if strcmp(grasp_env.inst,'ILL_d33') && strcmp(grasp_env.inst_option,'D33')
    if det ==1; %Rear
        cx_eff = cm.cm_pixels(1);
        cy_eff = cm.cm_pixels(2);
    elseif det == 2; % Right
        cx_eff = cm.cm_pixels(1) - ((params(inst_params.vectors.oxr) - cm.cm_translation(1)))/inst_params.(['detector' num2str(det)]).pixel_size(2); %horizontal distance from beam centre to pixel (m)
        cy_eff = cm.cm_pixels(2);
    elseif det == 3; % Left
        cx_eff = cm.cm_pixels(1) + ((params(inst_params.vectors.oxl) - cm.cm_translation(1)))/inst_params.(['detector' num2str(det)]).pixel_size(2); %horizontal distance from beam centre to pixel (m)
        cy_eff = cm.cm_pixels(2);
    elseif det == 4; %Bottom
        cx_eff = cm.cm_pixels(1);
        cy_eff = cm.cm_pixels(2)  + ((params(inst_params.vectors.oyb) - cm.cm_translation(2)))/inst_params.(['detector' num2str(det)]).pixel_size(2); %vertical distance from beam centre to pixel (m)
    elseif det == 5; %Top
        cx_eff = cm.cm_pixels(1);
        cy_eff = cm.cm_pixels(2) - ((params(inst_params.vectors.oyt) - cm.cm_translation(2)))/inst_params.(['detector' num2str(det)]).pixel_size(2); %vertical distance from beam centre to pixel (m)
    end
end

        
    
    
    
    mirrors = status_flags.analysis_modules.sectors.mirror_sectors;
    theta_set = status_flags.analysis_modules.sectors.theta;
    for n = 1:mirrors
        theta = theta_set + (n-1)*(360/mirrors);
        theta = check_angle(theta);
        
        sector_handles = circle(det,eff_outer_radius,....
           eff_inner_radius,....
            cx_eff,cy_eff,....
            theta-status_flags.analysis_modules.sectors.delta_theta/2,....
            theta+status_flags.analysis_modules.sectors.delta_theta/2,....
            status_flags.analysis_modules.sectors.sector_color,....
            status_flags.analysis_modules.sectors.anisotropy,....
            status_flags.analysis_modules.sectors.anisotropy_angle);
        grasp_handles.window_modules.sector.sketch_handles = [grasp_handles.window_modules.sector.sketch_handles; sector_handles];
    end
end


function theta = check_angle(theta)

while theta < 0
    theta = theta + 360;
end
while theta >= 360
    theta = theta - 360;
end

