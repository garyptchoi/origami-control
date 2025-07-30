clear; close all;

file_list = {'miura-ori';'recweave'; 'auxetic_triangle'; 'HexTri'; 'honey_comb'; ...
             'huffman_waterbomb'; 'lang_honey'; 'langoval'; 'leafy'; 'triPerfTess'};

for i = 1
    file = file_list{i};  
    
    fprintf('Processing file: %s\n', file);

if strcmp(file, 'recweave')
    plot_name = 'H. Rect Weave';
elseif strcmp(file, 'auxetic_triangle')
    plot_name = 'Auxetic Tri';
elseif strcmp(file, 'HexTri')
    plot_name = 'Hex/Tri';
elseif strcmp(file, 'honey_comb')
    plot_name = 'Kiri. Honeycomb';
elseif strcmp(file, 'huffman_waterbomb')
    plot_name = 'H. Waterbombs';
elseif strcmp(file, 'lang_honey')
    plot_name = 'Lang Honeycomb';
elseif strcmp(file, 'langoval')
    plot_name = 'Lang Oval';
elseif strcmp(file, 'leafy')
    plot_name = 'Leafy';
elseif strcmp(file, 'triPerfTess')
    plot_name = 'Perforated Tri';
elseif strcmp(file, 'miura-ori')
    plot_name = 'Miura-Ori';
else
    plot_name = 'Unknown Pattern';
end

if strcmp(file, 'recweave') || strcmp(file, 'huffman_waterbomb')
    L_all = [22,33,44];

else
    L_all = [2,3,4];
end


k_all = [1,2,4,8,16,32];     % Number of sampled constraints per step
rule_all = [1,2];                % Selection rules
n_sim = 100;                      % Number of simulations per setup
cmap = color_scheme(length(k_all));  % Color map for k values


%% Plot P vs rho for each L and rule
for rule = rule_all
    for L = L_all
        figure; hold on;

        for idx = 1:length(k_all)
            k = k_all(idx);

            filename = sprintf('simulation_results/%s/%s_%d_25%%_data/%s_%d_25%%_data_k_%d_rule_%d_n_100.mat', ...
                               file, file, L, file, L, k, rule);

            if exist(filename, 'file')              
                load(filename);  
                [n_sim, max_steps] = size(dof_all);
                r = linspace(0,1,max_steps);  
                dof_all = cummin(dof_all, 2);

                min_dof = max(dof_all(:, end));
                P = sum(dof_all <= min_dof, 1) / n_sim;

                plot(r, P, 'Color', cmap(idx,:), 'LineWidth', 2, ...
                     'DisplayName', ['k = ', num2str(k)]);
            else
                warning('File not found: %s', filename);
            end
        end

        title(sprintf('%s, Total Facets = %d, Rule %d',plot_name, max_steps-1, rule));
        xlabel('$\rho$', 'Interpreter', 'latex');
        ylabel('$P$', 'Interpreter', 'latex');
        set(gca, 'FontSize', 16, 'LineWidth', 1.5);
        axis([0 1 0 1]); grid off;
       

        save_folder = 'plots_P_k';
        if ~exist(save_folder, 'dir')
            mkdir(save_folder);
        end
        save_filename = sprintf('%s/P_%s_L_%d_rule_%d.svg', save_folder, file, L, rule);
        saveas(gcf, save_filename);

        close all;
    end   
end


end

% %% Plot k vs critical r for each L and rule
% r_star_all = zeros(length(L_all), length(k_all));
% r_0_all    = zeros(length(L_all), length(k_all));
% r_1_all    = zeros(length(L_all), length(k_all));
% r_01_all   = zeros(length(L_all), length(k_all));
% r_09_all   = zeros(length(L_all), length(k_all));
% r_025_all  = zeros(length(L_all), length(k_all));
% r_075_all  = zeros(length(L_all), length(k_all));
% 
% for rule = [2, 1]  % Plot rule 2 first
%     for L = L_all
%         id_L = find(L_all == L);
% 
%         for idx = 1:length(k_all)
%             k = k_all(idx);
%             id_k = find(k_all == k);
% 
%             filename = sprintf('simulation_results_office/%s/%s_%d_25%%_data/%s_%d_25%%_data_k_%d_rule_%d_n_100.mat', ...
%                                file, file, L, file, L, k, rule);
% 
%             if exist(filename, 'file')              
%                 load(filename);  % loads dof_all
%                 [n_sim, max_steps] = size(dof_all);
%                 r = linspace(0,1,max_steps);
% 
%                 min_dof = min(dof_all(:));
%                 P = sum(dof_all == min_dof, 1) / n_sim;
% 
%                 % Record critical rho thresholds
%                 idx_50 = find(P >= 0.5, 1);
%                 if ~isempty(idx_50), r_star_all(id_L, id_k) = r(idx_50); end
% 
%                 idx_0 = find(P == 0, 1, 'last');
%                 if ~isempty(idx_0), r_0_all(id_L, id_k) = r(idx_0); end
% 
%                 idx_1 = find(P == 1, 1);
%                 if ~isempty(idx_1), r_1_all(id_L, id_k) = r(idx_1); end
% 
%                 idx_01 = find(P <= 0.1, 1, 'last');
%                 if ~isempty(idx_01), r_01_all(id_L, id_k) = r(idx_01); end
% 
%                 idx_09 = find(P >= 0.9, 1);
%                 if ~isempty(idx_09), r_09_all(id_L, id_k) = r(idx_09); end
% 
%                 idx_025 = find(P <= 0.25, 1, 'last');
%                 if ~isempty(idx_025), r_025_all(id_L, id_k) = r(idx_025); end
% 
%                 idx_075 = find(P >= 0.75, 1);
%                 if ~isempty(idx_075), r_075_all(id_L, id_k) = r(idx_075); end
%             else
%                 warning('File not found: %s', filename);
%             end
%         end
% 
%         % Plot k vs rho*
%         figure;
%         hold on;
%         if rule == 1
%             plot([0, 34], [(4*L - 4)/L^2, (4*L - 4)/L^2], 'k--', 'LineWidth', 1);
%         end
% 
%         plot(k_all, r_star_all(id_L,:), '-o', ...
%              'Color', [201 0 22]/255, 'MarkerFaceColor', [201 0 22]/255, ...
%              'LineWidth', 2);
% 
%         xlabel('k');
%         ylabel('$\rho^*$', 'Interpreter', 'latex');
%         title(['L = ', num2str(L), ', Rule ', num2str(rule)]);
%         set(gca, 'FontSize', 18, 'LineWidth', 2);
%         xlim([0 34]); ylim([0 1]);
%     end
% end
% 
% 
% 



