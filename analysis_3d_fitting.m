clear; close all;

%% === 3D plot ===

cmap = [...
    0.8500, 0.3250, 0.0980;
    0.4660, 0.6740, 0.1880;
    0.0000, 0.4470, 0.7410;
    0.4940, 0.1840, 0.5560;
    0.9290, 0.6940, 0.1250;
    0.3010, 0.7450, 0.9330;
    0.0000, 0.5000, 0.0000;
    0.6350, 0.0780, 0.1840;
    1.0000, 0.4118, 0.1608;
];

file_list = {'miura-ori';
              'recweave'; 
              'lang_honey';
             'auxetic_triangle'; 
             'HexTri'; ...
             'honey_comb'; 
             'huffman_waterbomb'; 
             'langoval'; 
             'triPerfTess'
             };

load_folder = 'P_summary_data';
rule_seq = [1, 2];

figure('Visible','on', 'Position', [0, 0, 300, 250]);
hold on;

%[-90,0]
%[-40, 20];

view([-40, 20]);
rotate3d on;
grid on;

plot_data = containers.Map;

for rule = rule_seq
    if rule == 1
        k_seq = [32, 16, 8, 4, 2, 1];  
    else
        k_seq = [2, 4, 8, 16, 32];  
    end

    for k_idx = 1:length(k_seq)
        k = k_seq(k_idx);
        slice_y = -log2(k);  
        if rule == 2
            slice_y = log2(k);  
        end

        for i = 1:numel(file_list)
            file = file_list{i};
            base_color = cmap(i, :);

            switch file
                case 'recweave'
                    plot_name = 'H. Rect Weave'; L_all = [22, 33, 44]; marker_shape = 's';
                case 'auxetic_triangle'
                    plot_name = 'Auxetic Tri'; L_all = [2, 3, 4]; marker_shape = '^';
                case 'sim'
                    plot_name = 'Miura-Ori'; L_all = [2, 3, 4]; marker_shape = 's';
                case 'HexTri'
                    plot_name = 'Hex/Tri'; L_all = [2, 3, 4]; marker_shape = 'h';
                case 'huffman_waterbomb'
                    plot_name = 'H. Waterbombs'; L_all = [22, 33, 44]; marker_shape = 'v';
                case 'lang_honey'
                    plot_name = 'L. Honeycomb'; L_all = [2, 3, 4]; marker_shape = 'o';
                case 'langoval'
                    plot_name = 'Lang Oval'; L_all = [2, 3, 4]; marker_shape = 'o';
                case 'triPerfTess'
                    plot_name = 'Perforated Tri'; L_all = [2, 3, 4]; marker_shape = 'v';
                case 'honey_comb'
                    plot_name = 'Kiri. Honeycomb'; L_all = [2, 3, 4]; marker_shape = 's';
                otherwise
                    plot_name = 'Unknown'; L_all = [2, 3, 4]; marker_shape = 'o';
            end

            for L = L_all
               
                if k == 1
                    load_filename1 = sprintf('%s/P_%s_L%d_k%d_rule1.mat', load_folder, file, L, k);
                    load_filename2 = sprintf('%s/P_%s_L%d_k%d_rule2.mat', load_folder, file, L, k);

                    if ~exist(load_filename1, 'file') || ~exist(load_filename2, 'file')
                        continue;
                    end

                    load(load_filename1);  
                    idx_crit_1 = find(P >= 0.5, 1);
                    if isempty(idx_crit_1), continue; end
                    critical_r_1 = r(idx_crit_1);
                    tri_ratio_avg = tri_ratio;

                    load(load_filename2);  
                    idx_crit_2 = find(P >= 0.5, 1);
                    if isempty(idx_crit_2), continue; end
                    critical_r_2 = r(idx_crit_2);

                    critical_r = (critical_r_1 + critical_r_2) / 2;
                    tri_ratio = tri_ratio_avg;  
                else
                    load_filename = sprintf('%s/P_%s_L%d_k%d_rule%d.mat', ...
                                            load_folder, file, L, k, rule);
                    if ~exist(load_filename, 'file')
                        continue;
                    end

                    load(load_filename);  

                    idx_crit = find(P >= 0.5, 1);
                    if isempty(idx_crit), continue; end
                    critical_r = r(idx_crit);
                end

                if ismember(L, [2, 22, 10])
                    brightness = 0.5;
                elseif ismember(L, [4, 44, 20])
                    brightness = 1;
                else
                    brightness = 0.75;
                end
                color_final = base_color * brightness + (1 - brightness);

  
                scatter3(tri_ratio, slice_y * ones(size(tri_ratio)), critical_r, ...
                         80, marker_shape, ...
                         'MarkerFaceColor', color_final, ...
                         'MarkerEdgeColor', 'w');
                
                key = sprintf('%s_L%d', plot_name, L);
                if ~isKey(plot_data, key)
                    plot_data(key) = [];
                end
                pts = plot_data(key);
                new_pts = [tri_ratio(:), ...
                           slice_y * ones(numel(tri_ratio), 1), ...
                           critical_r(:), ...
                           repmat(color_final, numel(tri_ratio), 1)];
                plot_data(key) = [pts; new_pts];
            end
        end
    end
end


legend_keys = keys(plot_data);
for i = 1:length(legend_keys)
    key = legend_keys{i};
    pts = plot_data(key);
    if size(pts, 1) > 1
        pts = sortrows(pts, 2);  
        xyz = pts(:, 1:3);
        rgb = pts(1, 4:6);
        plot3(xyz(:, 1), xyz(:, 2), xyz(:, 3), '-', ...
              'Color', rgb, 'LineWidth', 2);
    end
end


set(gca, 'FontSize', 16, 'LineWidth', 1.5);
xlim([0 1]); ylim([-5 5]); zlim([0 1.12]);
axis vis3d;
grid off;


%% === Fit rho^* = a * tanh(b * y) + c * x + d ===

fprintf('\nPattern-by-pattern fitting results:\n');

pattern_keys = keys(plot_data);
base_map = containers.Map();  

% Group resolutions by base pattern
for i = 1:length(pattern_keys)
    full_key = pattern_keys{i};
    parts = split(full_key, '_');
    base_key = parts{1};
    pts = plot_data(full_key);
    if isKey(base_map, base_key)
        base_map(base_key) = [base_map(base_key); pts];
    else
        base_map(base_key) = pts;
    end
end

param_all = [];
label_all = {};
marker_all = {};
color_all = [];
rmse_all = [];

% Fit model for each base pattern
base_keys = keys(base_map);
for i = 1:length(base_keys)
    key = base_keys{i};
    pts = base_map(key);
    X_local = pts(:,1);
    Y_local = pts(:,2);
    Z_local = pts(:,3);
    color_final = mean(pts(:,4:6), 1);
    marker_shape = 'o';

    if contains(key, 'Miura') || contains(key, 'Rect') || contains(key, 'Kiri')
        marker_shape = 's';
    elseif contains(key, 'Auxetic')
        marker_shape = '^';
    elseif contains(key, 'Hex')
        marker_shape = 'h';
    elseif contains(key, 'Waterbomb') || contains(key, 'Perforated')
        marker_shape = 'v';
    elseif contains(key, 'Oval') || contains(key, 'L. Honey')
        marker_shape = 'o';
    end

    xy_local = [X_local, Y_local];
    model_fun = @(params, input) ...
        params(1) * tanh(params(2) * input(:,2) + params(3)) + ...
        params(4) * input(:,1) + params(5);
    init_guess = [0.5, 0.5, 0, 0, 0];
    opts = optimoptions('lsqcurvefit', 'Display', 'off');

    try
        params_fit = lsqcurvefit(model_fun, init_guess, xy_local, Z_local, [], [], opts);
        Z_fit_local = model_fun(params_fit, xy_local);
        rmse_local = sqrt(mean((Z_fit_local - Z_local).^2));

        fprintf('%s:\n', key);
        fprintf('  a = %.6f\n  b = %.6f\n  c = %.6f\n  d = %.6f\n  f = %.6f\n', ...
                params_fit(1), params_fit(2), params_fit(3), params_fit(4), params_fit(5));
        fprintf('  RMSE = %.6f\n\n', rmse_local);

        param_all = [param_all; params_fit];
        rmse_all = [rmse_all; rmse_local];
        label_all{end+1} = key;
        marker_all{end+1} = marker_shape;
        color_all = [color_all; color_final];

    catch
        fprintf('%s: Fitting failed.\n', key);
    end
end

% === 3D Plot ===
figure('Position', [0, 0, 300, 250]);
hold on; 
view([-40, 20]);
set(gca, 'FontSize', 16, 'LineWidth', 1.5);
xlim([0 1]); ylim([-5 5]); zlim([0 1.12]);
axis vis3d; pbaspect([1 1 1]); rotate3d on; grid off;

for i = 1:length(label_all)
    key = label_all{i};
    pts = base_map(key);
    X_local = pts(:,1);
    Y_local = pts(:,2);
    color_local = mean(pts(:,4:6), 1);
    marker_shape = marker_all{i};
    params = param_all(i, :);

    model_fun_local = @(input) ...
        params(1) * tanh(params(2) * input(:,2) + params(3)) + ...
        params(4) * input(:,1) + params(5);
    xy_local = [X_local, Y_local];
    Z_fit_local = model_fun_local(xy_local);

end

% smooth curves within each resolution
legend_keys = keys(plot_data);
for j = 1:length(legend_keys)
    subkey = legend_keys{j};
    pts = plot_data(subkey);
    if size(pts, 1) <= 1, continue; end

    pts = sortrows(pts, 2);  % by Y
    X_local = pts(:,1);
    Y_local = pts(:,2);
    color_local = pts(1, 4:6);

    Y_dense = linspace(min(Y_local), max(Y_local), 100)';
    X_dense = interp1(Y_local, X_local, Y_dense, 'linear');
    xy_dense = [X_dense, Y_dense];

    parts = split(subkey, '_');
    base_key = parts{1};
    idx = find(strcmp(label_all, base_key));
    if isempty(idx), continue; end

    params = param_all(idx, :);
    model_fun_local = @(input) ...
        params(1) * tanh(params(2) * input(:,2) + params(3)) + ...
        params(4) * input(:,1) + params(5);
    Z_dense = model_fun_local(xy_dense);

    plot3(X_dense, Y_dense, Z_dense, '-', ...
        'Color', color_local, ...
        'LineWidth', 3.5);
end

% === Print again summary ===
fprintf('\n=== Fitted Parameters and RMSE per Pattern ===\n');
for i = 1:length(label_all)
    key = label_all{i};
    params = param_all(i, :);
    rmse = rmse_all(i);
    fprintf('%s:\n', key);
    fprintf('  a = %.6f\n  b = %.6f\n  c = %.6f\n  d = %.6f\n  f = %.6f\n', ...
            params(1), params(2), params(3), params(4), params(5));
    fprintf('  RMSE = %.6f\n\n', rmse);
end


%% Legends information

% Base color for each pattern
cmap = [
    0.8500, 0.3250, 0.0980;
    0.4660, 0.6740, 0.1880;
    0.0000, 0.5000, 0.0000;
    0.6350, 0.0780, 0.1840;
    0.9290, 0.6940, 0.1250;
    0.0000, 0.4470, 0.7410;
    0.3010, 0.7450, 0.9330;
    1.0000, 0.4118, 0.1608;
    0.4940, 0.1840, 0.5560;
];

file_list = {'Miura-Ori', 'H. Rect Weave','H. Waterbombs', ...
             'Lang Oval', 'Hex/Tri', 'L. Honeycomb', 'Kiri. Honeycomb', ...
              'Perforated Tri','Auxetic Tri'};

marker_list = {'s','s','v','o','h','o','s','v','^'};
brightness_levels = [0.5, 0.75, 1];

figure('Position', [0, 0, 100, 200]); hold on;

for i = 1:numel(file_list)
    base_color = cmap(i,:);
    marker = marker_list{i};

    for j = 1:3  
        brightness = brightness_levels(j);
        color_final = base_color * brightness + (1 - brightness);
        
        if j == 3
            label = file_list{i};
        else
            label = '';  
        end

        plot(NaN, NaN, marker, ...
            'MarkerFaceColor', color_final, ...
            'MarkerEdgeColor', 'w', ...
            'MarkerSize', 10, ...
            'DisplayName', label);
    end
end

lgd = legend('Location', 'northoutside', ...
             'Orientation', 'horizontal', ...
             'NumColumns', 3);
set(lgd, 'FontSize', 11);

axis off;
title('Legend of Origami Structures', ...
      'FontWeight', 'bold', 'FontSize', 16)




