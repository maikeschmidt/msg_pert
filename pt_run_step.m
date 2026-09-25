function ok = pt_run_step(script_name)
% pt_run_step - Run one pipeline script in an isolated workspace
%
% Every msg_pert analysis script starts with clearvars. run() executes a
% script in the CALLER's workspace, so calling one directly from the master
% loop would wipe the loop's own variables. Running it through this function
% puts the clearvars in this function's workspace instead, and turns an error
% into a warning so one failed step does not stop the pipeline.
%
% USAGE:
%   ok = pt_run_step('pt_compute_metrics.m')
%
% OUTPUT:
%   ok  - true if the script finished without error
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk

try
    run_isolated(script_name);
    ok = true;
catch me
    ok = false;
    fprintf(2, 'WARNING: %s failed:\n  %s\n', script_name, me.message);
    if ~isempty(me.stack)
        fprintf(2, '  at %s, line %d\n', me.stack(1).name, me.stack(1).line);
    end
end
end


function run_isolated(script_file)
% The script's clearvars lands here, leaving pt_run_step's variables intact.
    run(script_file);
end
