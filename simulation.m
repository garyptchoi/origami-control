% Main Simulation Code for DoF calculation

%% ===================== Load Origami Data =====================
clear;

filenames = {
    % 'auxetic_triangle/auxetic_triangle_2_25%_data';
    % 'auxetic_triangle/auxetic_triangle_2_50%_data';
    % 'auxetic_triangle/auxetic_triangle_2_75%_data';
    % 'auxetic_triangle/auxetic_triangle_3_25%_data';
    % 'auxetic_triangle/auxetic_triangle_3_50%_data';
    % 'auxetic_triangle/auxetic_triangle_3_75%_data';
    % 'auxetic_triangle/auxetic_triangle_4_25%_data';
    % 'auxetic_triangle/auxetic_triangle_4_50%_data';
    % 'auxetic_triangle/auxetic_triangle_4_75%_data';
    % 
    % 'honey_comb/honey_comb_2_25%_data';
    % 'honey_comb/honey_comb_2_50%_data';
    % 'honey_comb/honey_comb_2_75%_data';
    % 'honey_comb/honey_comb_3_25%_data';
    % 'honey_comb/honey_comb_4_25%_data';

    % 'lang_honey/lang_honey_2_25%_data';
    % 'lang_honey/lang_honey_2_50%_data';
    % 'lang_honey/lang_honey_2_75%_data';
    % 'lang_honey/lang_honey_3_25%_data';
    % 'lang_honey/lang_honey_3_50%_data';
    % 'lang_honey/lang_honey_3_75%_data';
    % 'lang_honey/lang_honey_4_25%_data';
    % 'lang_honey/lang_honey_4_50%_data';
    % 'lang_honey/lang_honey_4_75%_data';

    % 'langoval/langoval_2_25%_data';
    % 'langoval/langoval_2_50%_data';
    % 'langoval/langoval_2_75%_data';
    % 'langoval/langoval_3_25%_data';
    % 'langoval/langoval_3_50%_data';
    % 'langoval/langoval_3_75%_data';
    % 'langoval/langoval_4_25%_data';
    % 'langoval/langoval_4_50%_data';
    % 'langoval/langoval_4_75%_data';


    % 'triPerfTess/triPerfTess_2_25%_data';
    % 'triPerfTess/triPerfTess_2_50%_data';
    % 'triPerfTess/triPerfTess_2_75%_data';
    % 'triPerfTess/triPerfTess_3_25%_data';
    % 'triPerfTess/triPerfTess_3_50%_data';
    % 'triPerfTess/triPerfTess_3_75%_data';
    % 'triPerfTess/triPerfTess_4_25%_data';
    % 'triPerfTess/triPerfTess_4_50%_data';
    % 'triPerfTess/triPerfTess_4_75%_data';
    % 
    % 'HexTri/HexTri_2_25%_data';
    % 'HexTri/HexTri_2_50%_data';
    % 'HexTri/HexTri_2_75%_data';
    % 'HexTri/HexTri_3_25%_data';
    % 'HexTri/HexTri_3_50%_data';
    % 'HexTri/HexTri_3_75%_data';
    % 'HexTri/HexTri_4_25%_data';
    % 'HexTri/HexTri_4_50%_data';
    % 'HexTri/HexTri_4_75%_data';

    'miura-ori/miura-ori_2_25%_data';
    'miura-ori/miura-ori_2_50%_data';
    'miura-ori/miura-ori_2_75%_data';
    'miura-ori/miura-ori_3_25%_data';
    'miura-ori/miura-ori_3_50%_data';
    'miura-ori/miura-ori_3_75%_data';
    'miura-ori/miura-ori_4_25%_data';
    'miura-ori/miura-ori_4_50%_data';
    'miura-ori/miura-ori_4_75%_data';


};

for i = 1:length(filenames)
    filename = filenames{i};
    fprintf('\nProcessing: %s\n', filename);

    load(['pre_process/origami_data/', filename, '.mat']);  % loads V, E, F, L

    edge_num = size(E, 1); 
    vertex_num = size(V, 1);
    facet_num = size(F, 1);

%% ===================== Parameters =====================
k_all = [1,2,4,8,16,32];          % Number of constraints to sample per step
rule_all = [1,2];                 % 1 = most efficient, 2 = least efficient
n_sim = 100;                      % Number of simulations

out_folder = ['simulation_results/', filename];
if ~exist(out_folder, 'dir')
    mkdir(out_folder);
end

% ===================== Facet Adjacency =====================
% Build facet adjacency graph via label-1 (planar diagonal) edges
facet_groups = num2cell(1:facet_num);  

% Build edge to facet map
edge_to_facet = containers.Map; 

for i = 1:facet_num
    f = F(i, :);
    edges_in_facet = sortrows([f([1 2]); f([2 3]); f([3 1])], 2);  
    for j = 1:3
        e = sort(edges_in_facet(j, :));
        key = sprintf('%d_%d', e(1), e(2));
        if isKey(edge_to_facet, key)
            edge_to_facet(key) = [edge_to_facet(key), i];
        else
            edge_to_facet(key) = i;
        end
    end
end

% Merge adjacent facets if they share a label-1 edge
for i = 1:edge_num
    if L(i) == 1  
        e = sort(E(i, :));
        key = sprintf('%d_%d', e(1), e(2));
        if isKey(edge_to_facet, key)
            adjacent_facets = edge_to_facet(key);
            if numel(adjacent_facets) == 2

                f1 = adjacent_facets(1);
                f2 = adjacent_facets(2);
                g1 = find(cellfun(@(g) ismember(f1, g), facet_groups));
                g2 = find(cellfun(@(g) ismember(f2, g), facet_groups));

                if g1 ~= g2
                    facet_groups{g1} = unique([facet_groups{g1}, facet_groups{g2}]);
                    facet_groups(g2) = [];
                end
            end
        end
    end
end



%% ===================== Visualization of Merged Facet Groups =====================
% figure;
% hold on;
% axis equal;
% xlabel('X'); ylabel('Y'); zlabel('Z');
% title('Merged Facet Groups', 'Interpreter', 'latex', 'FontSize', 14);
% view(3);
% 
% % Generate random colors for each group
% rng(1);  % For reproducibility
% group_colors = rand(numel(facet_groups), 3);
% 
% for i = 1:numel(facet_groups)
%     group_facet_indices = facet_groups{i};
%     patch('Faces', F(group_facet_indices, :), ...
%           'Vertices', V, ...
%           'FaceColor', group_colors(i, :), ...
%           'EdgeColor', [0.3, 0.3, 0.3], ...
%           'FaceAlpha', 0.8, ...
%           'LineWidth', 0.5);
% end

%% ===================== DoF =====================

for k = k_all
    for rule = rule_all

        max_steps = length(facet_groups);
        dof_all = zeros(n_sim, max_steps + 1);

        parfor sim = 1:n_sim
            fprintf('Sim %d/%d | k=%d, rule=%d\n', sim, n_sim, k, rule);
            tic;

            % ===== Initial Rigidity Matrix =====
            [rows, cols, vals] = deal([]);
            row_idx = 1;
            for i = 1:edge_num
                v1 = E(i, 1); v2 = E(i, 2);
                col = [3*v1-2, 3*v1-1, 3*v1, 3*v2-2, 3*v2-1, 3*v2];
                diff = V(v1, :) - V(v2, :);
                val = [2 * diff, -2 * diff];
                rows = [rows, repmat(row_idx, 1, 6)];
                cols = [cols, col];
                vals = [vals, val];
                row_idx = row_idx + 1;
            end

            H = sparse(rows, cols, vals, row_idx - 1, 3 * vertex_num);
            DoF_initial = 3 * vertex_num - calc_rank(H) - 6;

            used = false(length(facet_groups), 1);
            dof_list = zeros(1, max_steps + 1);
            dof_list(1) = DoF_initial;
            step = 1;

            while any(~used)
                available_indices = find(~used);
                num_choices = min(k, length(available_indices));
                random_indices = available_indices(randperm(length(available_indices), num_choices));
                DoF_process = zeros(num_choices, 1);

                for p = 1:num_choices
                    facet_id = random_indices(p);
                    %===========================================================
                    R_new = Compute_row_cross(V, E, L, F, facet_groups{facet_id}); 
                    %===========================================================
                    A = [H; R_new];
                    DoF = 3 * vertex_num - calc_rank(A) - 6;
                    DoF = min(DoF, DoF_initial);
                    DoF = max(DoF, 1);
                    DoF_process(p) = DoF;
                end

                % Selection rule
                if rule == 1
                    [~, best_idx] = min(DoF_process);
                else
                    [~, best_idx] = max(DoF_process);
                end

                selected_facet_idx = random_indices(best_idx);
                
                %===========================================================
                R_new = Compute_row_cross(V, E, L, F, facet_groups{selected_facet_idx});
                %===========================================================
                H = [H; R_new];

                dof = 3 * vertex_num - calc_rank(H) - 6;
                dof = min(dof, DoF_initial);
                dof = max(dof, 1);

                dof_list(step + 1) = dof;
                step = step + 1;
                used(selected_facet_idx) = true;
            end

            dof_all(sim, :) = dof_list;
            toc;

        end
        
        save_sim_origami(filename, k, rule, dof_all, n_sim) 

    end
end



end



%% ===================== Plot DoF Evolution =====================

% 
% base_color = [0.1, 0.45, 0.6];
% light_color = base_color + (1 - base_color) * 0.5;
% 
% figure;
% hold on;
% grid on;
% box on;
% 
% title(sprintf('DoF: %s\n$k = %d$, rule = %d', ...
%     strrep(filename, '_', '\_'), k, rule), ...
%     'Interpreter', 'latex', 'FontSize', 14);
% xlabel('Step', 'Interpreter', 'latex');
% ylabel('Degrees of Freedom', 'Interpreter', 'latex');
% 
% h_sim = plot(0:max_steps, dof_all(1, :), 'Color', light_color, 'LineWidth', 1.2);
% 
% for sim = 2:n_sim
%     plot(0:max_steps, dof_all(sim, :), 'Color', light_color, 'LineWidth', 1.2);
% end
% 
% % Plot average curve and capture its handle
% avg_curve = mean(dof_all, 1);
% h_avg = plot(0:max_steps, avg_curve, 'Color', base_color, 'LineWidth', 2.0);
% 
% % Create legend with correct handles and colors
% legend([h_sim, h_avg], {'Simulations', 'Average'}, ...
%        'Location', 'best', 'Interpreter', 'latex');
% 


