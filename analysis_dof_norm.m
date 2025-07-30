
clear; close all;

file_list = {'miura-ori';'recweave'; 'auxetic_triangle'; 'HexTri'; 'honey_comb'; ...
           'huffman_waterbomb'; 'lang_honey'; 'langoval'; 'triPerfTess'};

k_list = [1, 2, 4, 8, 16, 32];
rule_list = [1, 2];

n = length(k_list);
color_list = color_scheme(n);
style_list = {'-', '-', '-'};

for i = 1:numel(file_list)
    file = file_list{i};

    % Set plot name
    if strcmp(file, 'miura-ori')
        plot_name = 'Miura-Ori';
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
    elseif strcmp(file, 'recweave')
        plot_name = 'H. Rect Weave';
    elseif strcmp(file, 'triPerfTess')
        plot_name = 'Perforated Tri';
    else
        plot_name = 'Unknown Pattern';
    end

    % Set L_all
    if strcmp(file, 'recweave') || strcmp(file, 'huffman_waterbomb')
        L_all = [22,33,44];
    else
        L_all = [2,3,4];
    end

    for rule = rule_list
        for L = L_all
            figure; hold on;
            max_dof = 0;

            for idx = 1:length(k_list)
                k = k_list(idx);

                filename = sprintf('simulation_results/%s/%s_%d_25%%_data/%s_%d_25%%_data_k_%d_rule_%d_n_100.mat', ...
                                   file, file, L, file, L, k, rule);

                if exist(filename, 'file')
                    load(filename);  % loads dof_all
                    [n_sim, max_steps] = size(dof_all);
                    rho = (0:max_steps-1)/(max_steps-1);

                    DoF_initial = max(dof_all(:, 1));
                    dof_all_norm = (dof_all - 1) / (DoF_initial - 1);

                    % Plot each individual curve
                    light_color = color_list(idx,:) + (1 - color_list(idx,:)) * 0.1;
                    for sim = 1:n_sim
                        plot(rho, dof_all_norm(sim,:), ...
                             'Color', [light_color, 0.1], 'LineWidth', 1, ...
                             'HandleVisibility', 'off');
                    end

                    % Plot average
                    avg_curve = mean(dof_all_norm, 1);
                    plot(rho, avg_curve, ...
                         'Color', color_list(idx,:), ...
                         'LineStyle', style_list{mod(idx-1, length(style_list)) + 1}, ...
                         'LineWidth', 2.2, ...
                         'DisplayName', sprintf('k = %d', k));

                    max_dof = max(max_dof, max(avg_curve(:)));
                else
                    warning('File not found: %s', filename);
                end
            end

            % Labels and title
            xlabel('$\rho$', 'Interpreter', 'latex');
            ylabel('$\widetilde{d}$', 'Interpreter', 'latex');

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%

            xlim([0, 1]);        
            xticks([0 0.25 0.5 0.75 1]);


            ylim([0, 1]);
            yticks([0 0.25 0.5 0.75 1]);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            title(sprintf('%s', plot_name));
            set(gca, 'FontSize', 20, 'LineWidth', 1.5);

            grid off;

            % Set figure size and save
            set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 10.5, 8]);
            save_folder = 'plots_k_norm';
            if ~exist(save_folder, 'dir')
                mkdir(save_folder);
            end
           save_filename_svg = sprintf('%s/%s_L_%d_rule_%d.svg', save_folder, file, L, rule);
           print(gcf, save_filename_svg, '-dsvg');
           close all;
        end
    end
end



