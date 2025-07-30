function R = Compute_row_cross(V, E, L, F, facet_merged)
% Compute planarity constraints for a merged group of facets
% using edges labeled 1 (planar diagonal).
%
% Output R is a sparse matrix. Each row represents one planarity constraint.

% vp-- v2
% |  / |
% |/   |
% v1 --vq  

R =[];

    vertice_num = size(V, 1);
    R = sparse(0, 3 * vertice_num);  % initialize empty matrix

    if numel(facet_merged) == 1
        return;  % no constraint needed
    end

    % Step 1: collect all vertices in the merged group
    all_facet_vertices = unique(F(facet_merged, :));

    % Step 2: collect all label-1 edges from the merged facets
    label1_edges = [];
    for i = 1:size(E, 1)
        if L(i) == 1
            e = E(i, :);
            % if both ends are in the merged group
            if all(ismember(e, all_facet_vertices))
                label1_edges = [label1_edges; e];
            end
        end
    end

    % Step 3: for each label-1 edge, find the 2 facets from the merged group that contain it
    for i = 1:size(label1_edges, 1)
        edge = label1_edges(i, :);
        v_shared = sort(edge);

        % Find two facets in the merged group that contain the edge
        containing_facets = [];
        for f_idx = facet_merged
            f = F(f_idx, :);
            if sum(ismember(f, v_shared)) == 2
                containing_facets = [containing_facets; f];
            end
        end

        if size(containing_facets, 1) ~= 2
            continue;  % skip if not exactly two facets contain this edge
        end

        % Get the third vertex in each triangle (not in the edge)
        vp = setdiff(containing_facets(1, :), v_shared);
        vq = setdiff(containing_facets(2, :), v_shared);

        v1 = v_shared(1);
        v2 = v_shared(2);        

        if isempty(vp) || isempty(vq)
            continue;  % safety check
        end

        diff4 = V(vp,:) - V(v1,:);
        diff2 = V(vq,:) - V(v1,:);
        diff3 = V(v2,:) - V(v1,:);

        % Compute coefficients
        x21 = diff2(1);  x31 = diff3(1);  x41 = diff4(1);
        y21 = diff2(2);  y31 = diff3(2);  y41 = diff4(2);
        z21 = diff2(3);  z31 = diff3(3);  z41 = diff4(3);

        a1 = -(y21 * z41 - y41 * z21) + (y31 * z41 - y41 * z31) - (y31 * z21 - y21 * z31);
        a2 =  (x21 * z41 - x41 * z21) - (x31 * z41 - x41 * z31) + (x31 * z21 - x21 * z31);
        a3 = -(x21 * y41 - x41 * y21) + (x31 * y41 - x41 * y31) - (x31 * y21 - x21 * y31);

        b1 = -(y31 * z41 - y41 * z31);
        b2 =  x31 * z41 - x41 * z31;
        b3 = -(x31 * y41 - x41 * y31);

        c1 =  y21 * z41 - y41 * z21;
        c2 = -(x21 * z41 - x41 * z21);
        c3 =  x21 * y41 - x41 * y21;

        d1 =  y31 * z21 - y21 * z31;
        d2 = -(x31 * z21 - x21 * z31);
        d3 =  x31 * y21 - x21 * y31;

        % Indices of the 4 involved vertices: v1, vp, v2, vq
        col = [3*v1-2, 3*v1-1, 3*v1, ...
               3*vq-2, 3*vq-1, 3*vq, ...
               3*v2-2, 3*v2-1, 3*v2, ...
               3*vp-2, 3*vp-1, 3*vp];

        val = [a1 a2 a3 b1 b2 b3 c1 c2 c3 d1 d2 d3];

        row = [ones(1, 12)];

        R = [R; sparse(row, col, val, 1, 3 * vertice_num)];


    end
end

