# MFF
Multi-file Fitter - automated data analysis within the Graphical Reduction and Analysis SANS Program (GRASP).
This software suite requires matlab and GRASP to run, the matlab version of GRASP is freely avaiable from the Institut Laue-Langevin:
https://www.ill.eu/instruments-support/instruments-groups/groups/lss/grasp/download/

To add Multi-File Fitter to GRASP, start by adding the MFF folder (all of the code stored in this repositor) to your MATLAB path. This can be done by right clicking on the folder in the Current Folder window and selecting Add to Path → Selected Folders and Subfolders.
Next, we need to add a button in GRASP so that we can run MFF. One simple way of doing this is to edit the modify_main_menu_items.m file, which can be found in the main_interface folder of GRASP. After opening this file, scroll down to the User Module section (approximately line 512) and add the following line of code:
uimenu(grasp_handles.menu.user_modules.root,'separator','on','label','M ulti-File Fitter','callback','mf_GUI_window','enable','on');
Save modify_main_menu_items.m after adding this line of code. This should be all the MATLAB coding required to use MFF. (Note: If you cannot/ do not want to modify GRASP code, MFF can be accessed and ran by returning the command mf_GUI_window in the Command Window while running GRASP).

MFF Data Structure Organization:
name  - mf_fitter

main fields
algorithm_options - stores all specifics related to the selected fitting algorithm, i.e. which fitting cycles should be viewed
fit_data - contains gaussian (or whatever fit type) fits to the data as a function of file number ("Numor"), which should directly  correspond with the control parameter.
handles - stores the pointers for all windows and input fields, such as the GUI window or the control parameter text entry box.
numors - file numors that the software should be run over.
save_options - includes the folder name and full path for storing all outputs from the software suite.  It also stores figure extension selections (i.e. should figures be stored as .eps, .fig, .jpg, and/or .pdf.
user_inputs - contains information input through the GUI by the user. Includes the subfields control_parameter, int_cutoff, dphi_cutoff, reference_files, and smoothing.  The control parameter is what is changing in each numor and will be the x-axis for final plotting, this could be number of applied AC cycles or temperature for example. The int_cutoff and dphi_cutoff specify what minimum peak intensity and peak separation is necessary for being included in code to determine guess values for future cycles. Reference_files are peaks that are clearly resolvable, in terms of e.g. the intensity and separation; fit parameters from these are used as initial guess values. Smoothing includes options for the fwhm and step size for binning / gaussian smoothing.

MFF Fitting Algorithms
mf_fitter_mpf_2K() - ideal for fitting the supercooled vortex lattice data where equilibrium state peaks nucleate at their final orientaiton and grow.
