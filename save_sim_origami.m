% Save the simulation results
%
% If you use this code in your work, please cite the following paper:
%    R. Li and G. P. T. Choi,
%    "Explosive rigidity percolation in origami."
%    Preprint, arXiv:2410.13945, 2024.
% 
% Copyright (c) 2024, Rongxuan Li and Gary P. T. Choi
% 
% https://github.com/garyptchoi/origami-explosive-percolation

function save_sim_origami(filename, k, rule, dof_all, n_sim)
 
    parts = strsplit(filename, '/');
    folder = parts{1};     
    name = parts{2};      

    save_folder = fullfile('simulation_results', folder, name);

    if ~exist(save_folder, 'dir')
        mkdir(save_folder);
    end

    save_path = fullfile(save_folder, ...
        [name, '_k_', num2str(k), '_rule_', num2str(rule), '_n_', num2str(n_sim), '.mat']);

    save(save_path, 'filename', 'dof_all', 'k', 'rule', 'n_sim');
end
