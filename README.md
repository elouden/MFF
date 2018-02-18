# MFF
Multi-file Fitter - automated data analysis within the Graphical Reduction and Analysis SANS Program (GRASP).
This software suite requires matlab and GRASP to run, the matlab version of GRASP is freely avaiable from the Institut Laue-Langevin:
https://www.ill.eu/instruments-support/instruments-groups/groups/lss/grasp/download/

To add Multi-File Fitter to GRASP, start by adding the MFF folder (all of the code stored in this repositor) to your MATLAB path. This can be done by right clicking on the folder in the Current Folder window and selecting Add to Path → Selected Folders and Subfolders.
Next, we need to add a button in GRASP so that we can run MFF. One simple way of doing this is to edit the modify_main_menu_items.m file, which can be found in the main_interface folder of GRASP. After opening this file, scroll down to the User Module section (approximately line 512) and add the following line of code:
uimenu(grasp_handles.menu.user_modules.root,'separator','on','label','M ulti-File Fitter','callback','mf_GUI_window','enable','on');
Save modify_main_menu_items.m after adding this line of code. This should be all the MATLAB coding required to use MFF. (Note: If you cannot/ do not want to modify GRASP code, MFF can be accessed and ran by returning the command mf_GUI_window in the Command Window while running GRASP).
