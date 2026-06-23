function p = pert_path()
% pert_path  Return the absolute path to the msg_pert repository root.
%
%   p = pert_path()
%
% Used internally by pt_add_functions to locate the toolbox root regardless
% of the current working directory.

p = fileparts(mfilename('fullpath'));
end
