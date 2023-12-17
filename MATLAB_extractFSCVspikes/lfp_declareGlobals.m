%lfp_declareGlobals
% A script that declares global all the variables listed in the cell row
% vector lfp_declareGlobals_global_names stored in the file
% 'lfp_globals.mat', and then assigns that list value to the global
% 'lfp_GlobalNames'.  That list is the canonical list of all globals used
% by lfp_lib.  Also sets the value of the global lfp_GlobalNames to a cell
% row vector containing the names of all the global variables.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

load('lfp_globals.mat');
lfp_declareGlobals_globals_string = '';
for lfp_declareGlobals_name = lfp_declareGlobals_global_names
    lfp_declareGlobals_globals_string = ...
        [lfp_declareGlobals_globals_string ' ' char(lfp_declareGlobals_name)];
end
lfp_declareGlobals_cmd = ['global ' lfp_declareGlobals_globals_string];
eval(lfp_declareGlobals_cmd);
lfp_GlobalNames = lfp_declareGlobals_global_names;
clear lfp_declareGlobals_globals_string lfp_declareGlobals_name ...
lfp_declareGlobals_cmd lfp_declareGlobals_global_names