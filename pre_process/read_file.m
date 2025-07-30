% Reading file: Origami visualization and information (vertices, edges, facets) saving code.

% Visualize and Save the Origami/Kirigami Vertices, Edges
% and Facets information from https://origamisimulator.org/

% Download the '.obj' file of the pattern from https://origamisimulator.org/ and change to '.m' .
% Input the file '.m' file at line 32.

% Output: matrix V (all vertices (x_i,y_i,z_i)), [size=n*3] n: # of vertices)
%         matrix F (all facets   (v1, v2, v3),   [size=n*3] n: # of facets)
%         matrix E (all edges    (v1,v2),        [size=n*2] n: # of edges) 
%         matrix L (edges' label (v1,v2),        [size=n*1] n: # of edges).
%
%         matrix L Labeling rule:
%         Label 0 (marginal edge)       : light gray
%         Label 1 (planar diagonal edge): red
%         Label 2 (mountain edge)       : green
%         Label 3 (valley edge)         : blue

%%

clear; close all;

vertices = [];
facets = [];
edges = [];
labels = [];

% List of input files
file_list = {
    'origami_source/miura-ori/miura-ori_2_25%.m';
    'origami_source/miura-ori/miura-ori_2_50%.m';
    'origami_source/miura-ori/miura-ori_2_75%.m';
    'origami_source/miura-ori/miura-ori_3_25%.m';
    'origami_source/miura-ori/miura-ori_3_50%.m';
    'origami_source/miura-ori/miura-ori_3_75%.m';
    'origami_source/miura-ori/miura-ori_4_25%.m';
    'origami_source/miura-ori/miura-ori_4_50%.m';
    'origami_source/miura-ori/miura-ori_4_75%.m';
};

for file_idx = 1:length(file_list)
    filename = file_list{file_idx};

    vertices = [];
    facets = [];
    edges = [];
    labels = [];

    fileID = fopen(filename, 'r');
    [~, name, ~] = fileparts(filename);
    name_disp = strrep(name, '_', '\_');

    % Read the file line by line
    tline = fgetl(fileID);
    while ischar(tline)
        if startsWith(tline, 'v')
            tline_clean = strtrim(erase(tline, 'v'));
            vertex = str2double(strsplit(tline_clean));
            if numel(vertex) == 3 && all(~isnan(vertex))
                vertices = [vertices; vertex];
            end
        elseif startsWith(tline, 'f')
            tline_clean = strtrim(erase(tline, 'f'));
            parts = strsplit(tline_clean);
            indices = [];
            for i = 1:length(parts)
                idx_parts = strsplit(parts{i}, '/');
                indices(end+1) = str2double(idx_parts{1});
            end
            if numel(indices) == 3 && all(~isnan(indices))
                facets = [facets; indices];
            end
        elseif startsWith(strtrim(tline), '#e')
            tline_clean = strtrim(erase(tline, '#e'));
            parts = str2double(strsplit(tline_clean));
            if numel(parts) >= 3 && all(~isnan(parts(1:3)))
                edges = [edges; parts(1:2)];
                labels = [labels; parts(3)];
            end
        end
        tline = fgetl(fileID);
    end
    fclose(fileID);

    % Save to variables
    V = vertices;
    F = facets;
    E = edges;
    L = labels;

    %% Plotting
    figure;
    hold on;

    patch('Faces', F, 'Vertices', V, ...
          'FaceColor', [1,1,1], 'EdgeColor', [1,1,1], ...
          'LineWidth', 0.5, 'FaceAlpha', 1.0);

    labelColors = [
        0.8 0.8 0.8;
        1.0 0.0 0.0;
        0.0 0.6 0.0;
        0.0 0.4 1.0;
    ];

    for i = 1:size(E,1)
        v1 = V(E(i,1), :);
        v2 = V(E(i,2), :);
        label = L(i);
        if label >= 0 && label <= 3
            color = labelColors(label + 1, :);
            plot3([v1(1), v2(1)], [v1(2), v2(2)], [v1(3), v2(3)], '-', ...
                  'Color', color, 'LineWidth', 1.5);
        end
    end

    minVals = min(V, [], 1);
    maxVals = max(V, [], 1);
    center = (minVals + maxVals) / 2;
    range = (maxVals - minVals) * 1.3 / 2;

    xlim([center(1) - range(1), center(1) + range(1)]);
    ylim([center(2) - range(2), center(2) + range(2)]);
    zlim([center(3) - range(3), center(3) + range(3)]);

    axis equal;
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(name_disp, 'Interpreter', 'latex', 'FontSize', 16); 
    view(3);

    h1 = plot3(nan, nan, nan, '-', 'Color', labelColors(1,:), 'LineWidth', 1.5);
    h2 = plot3(nan, nan, nan, '-', 'Color', labelColors(2,:), 'LineWidth', 1.5);
    h3 = plot3(nan, nan, nan, '-', 'Color', labelColors(3,:), 'LineWidth', 1.5);
    h4 = plot3(nan, nan, nan, '-', 'Color', labelColors(4,:), 'LineWidth', 1.5);

    legend([h1, h2, h3, h4], ...
           {'marginal edge', 'planar diagonal edge', 'mountain edge', 'valley edge'}, ...
           'Location', 'bestoutside', 'Interpreter', 'latex','FontSize', 10);

    % Save results
    name_save = strrep(name, '\_', '_');
    save(['origami_data/', name_save, '_data.mat'], 'V', 'E', 'F', 'L');
end



%%  Color Name	RGB Triplet

% Red	[1.0, 0.0, 0.0]
% Green	[0.0, 1.0, 0.0]
% Blue	[0.0, 0.0, 1.0]
% Yellow	[1.0, 1.0, 0.0]
% Cyan (Aqua)	[0.0, 1.0, 1.0]
% Magenta	[1.0, 0.0, 1.0]
% Orange	[1.0, 0.5, 0.0]
% Purple	[0.5, 0.0, 0.5]
% Pink	[1.0, 0.4, 0.7]
% Teal	[0.0, 0.5, 0.5]
% Light Gray	[0.8, 0.8, 0.8]
% Dark Gray	[0.3, 0.3, 0.3]
% Black	[0.0, 0.0, 0.0]
% White	[1.0, 1.0, 1.0]
% Sky Blue	[0.4, 0.7, 1.0]
