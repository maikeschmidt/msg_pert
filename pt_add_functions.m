% pt_add_functions - Initialise msg_pert toolbox and check dependencies
%
% Adds the msg_pert root to the MATLAB path and verifies that the required
% companion toolboxes (msg_coreg, msg_fwd, HBF, FieldTrip via SPM) are
% accessible.
%
% USAGE:
%   pt_add_functions
%
% DEPENDENCIES CHECKED:
%   msg_coreg    — sibling directory; provides cr_add_functions and HBF
%   msg_fwd      — sibling directory; provides BEM/FEM pipeline and functions/
%   SPM          — must be on path (provides FieldTrip wrappers)
%   HBF          — Helsinki BEM Framework (via msg_coreg/hbf_lc_p)
%
% REPOSITORY:
%   https://github.com/maikeschmidt/msg_pert
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk

root = pert_path();
addpath(root);
fprintf('msg_pert root added to path: %s\n', root);

% Check SPM
if isempty(which('spm'))
    warning('pt_add_functions: SPM not found on path. Add SPM before running scripts.');
else
    fprintf('SPM found: %s\n', fileparts(which('spm')));
end

% Check and initialise msg_coreg
coreg_path = fullfile(fileparts(root), 'msg_coreg');
if isfolder(coreg_path)
    coreg_init = fullfile(coreg_path, 'cr_add_functions.m');
    if exist(coreg_init, 'file')
        run(coreg_init);
        fprintf('msg_coreg initialised: %s\n', coreg_path);
    else
        addpath(coreg_path);
        fprintf('msg_coreg added to path: %s\n', coreg_path);
    end
else
    warning('pt_add_functions: msg_coreg not found at %s\nClone from https://github.com/maikeschmidt/msg_coreg', coreg_path);
end

% Check msg_fwd and add its functions/ subfolder
fwd_path = fullfile(fileparts(root), 'msg_fwd');
if isfolder(fwd_path)
    addpath(fwd_path);
    fwd_fns = fullfile(fwd_path, 'functions');
    if isfolder(fwd_fns)
        addpath(fwd_fns);
    end
    fprintf('msg_fwd added to path: %s\n', fwd_path);
else
    warning('pt_add_functions: msg_fwd not found at %s\nClone from https://github.com/maikeschmidt/msg_fwd', fwd_path);
end

% Check HBF
if isempty(which('hbf_CheckTriangleOrientation'))
    hbf_path = fullfile(fileparts(root), 'msg_coreg', 'hbf_lc_p');
    if isfolder(hbf_path)
        addpath(genpath(hbf_path));
        fprintf('HBF library added to path: %s\n', hbf_path);
    else
        warning('pt_add_functions: HBF library not found. Clone into msg_coreg/hbf_lc_p from https://github.com/MattiStenroos/hbf_lc_p');
    end
else
    fprintf('HBF library already on path.\n');
end

fprintf('pt_add_functions complete.\n');
