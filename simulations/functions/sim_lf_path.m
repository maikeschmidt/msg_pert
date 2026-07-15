function p = sim_lf_path(model, geo, array)
% sim_lf_path - Build the leadfield .mat path for a model + geometry + array
%
% Centralises the msg_fwd file-naming conventions so every sim script resolves
% leadfield files the same way, for any geometry stem.
%
% USAGE:
%   p = sim_lf_path(model, geo, array)
%
% INPUT:
%   model  one entry of sim_models — uses .method ('bem'/'bslaw'), .root
%   geo    one entry of sim_geometries — uses .short, .kind ('standard'/'cond'),
%          and for 'cond': .cond_bundle, .cond_shift, .root (override)
%   array  'front' or 'back'
%
% NAMING (matching msg_fwd / run_conductivity_perturbation output):
%   bem  standard : <root>/geometries_<short>/leadfield_<short>_bem_<array>.mat
%   bem  cond     : <croot>/geometries_<short>/leadfield_<short>_bem_cond_bundleB_shiftS_<array>.mat
%   bslaw         : <root>/leadfield_geometries_<short>_bslaw_<array>.mat   (flat folder)
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk

short = geo.short;

switch model.method
    case 'bslaw'
        % Biot-Savart has no cond variant and lives in a flat folder
        p = fullfile(model.root, ...
            sprintf('leadfield_geometries_%s_bslaw_%s.mat', short, array));

    case 'bem'
        if isfield(geo, 'kind') && strcmp(geo.kind, 'cond')
            % Conductivity leadfields live under the MODEL's per-modality cond
            % root (MSG and ESG each have their own). A geo.root override, if
            % set, still wins.
            if isfield(geo, 'root') && ~isempty(geo.root)
                root = geo.root;
            elseif isfield(model, 'cond_root') && ~isempty(model.cond_root)
                root = model.cond_root;
            else
                error(['sim_lf_path: model "%s" has no cond_root for the ' ...
                       'conductivity variant "%s".'], model.id, geo.name);
            end
            fname = sprintf('leadfield_%s_bem_cond_bundle%d_shift%d_%s.mat', ...
                short, geo.cond_bundle, geo.cond_shift, array);
        else
            root  = model.root;
            fname = sprintf('leadfield_%s_bem_%s.mat', short, array);
        end
        p = fullfile(root, ['geometries_' short], fname);

    otherwise
        error('sim_lf_path: unknown model method "%s".', model.method);
end
end
