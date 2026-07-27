function out = pt_modality(action, value)
% pt_modality - get/set the active modality ('msg' or 'esg') for the pipeline
%
% The msg_pert analysis scripts each call config_pert at their top (after a
% clearvars), so there is no ordinary workspace variable an outer loop could use
% to tell them which modality to run. This function holds the choice in a
% persistent variable instead — clearvars does NOT clear persistents, so
% run_perturbation_analysis can loop over MSG and ESG, setting the modality once
% per iteration, and every sub-script picks it up via config_pert.
%
% USAGE:
%   pt_modality('set', 'esg')   % choose modality (done by the master loop)
%   m = pt_modality('get')      % read it ('' if never set)
%   pt_modality('clear')        % reset (standalone runs then default to msg)
%
% If nothing has been set, 'get' returns '' and config_pert falls back to 'msg'.
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk

persistent M

if nargin == 0 || strcmp(action, 'get')
    if isempty(M); out = ''; else; out = M; end
    return
end

switch action
    case 'set'
        if ~ismember(value, {'msg', 'esg'})
            error('pt_modality: modality must be ''msg'' or ''esg'', got ''%s''.', value);
        end
        M   = value;
        out = M;
    case 'clear'
        M   = '';
        out = '';
    otherwise
        error('pt_modality: unknown action ''%s'' (use get/set/clear).', action);
end
end
