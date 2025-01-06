

function quadtree = build_quadtree(layout, max_points)
    [rows, cols] = find(layout > 0); % 找到暴露點
    points = [rows, cols]; % 點的座標
    quadtree = subdivide(points, [1, size(layout, 1)], [1, size(layout, 2)], max_points);
end

function node = subdivide(points, row_range, col_range, max_points)
    if size(points, 1) <= max_points
        node.isLeaf = true;
        node.points = points;
        node.row_range = row_range;
        node.col_range = col_range;
        return;
    end
    % 劃分為四個子區域
    row_mid = mean(row_range);
    col_mid = mean(col_range);
    regions = {
        points(points(:, 1) <= row_mid & points(:, 2) <= col_mid, :),
        points(points(:, 1) <= row_mid & points(:, 2) > col_mid, :),
        points(points(:, 1) > row_mid & points(:, 2) <= col_mid, :),
        points(points(:, 1) > row_mid & points(:, 2) > col_mid, :)
    };
    node.isLeaf = false;
    node.children = cell(1, 4);
    for i = 1:4
        node.children{i} = subdivide(regions{i}, ...
            [row_range(1), row_mid] * (i <= 2) + [row_mid, row_range(2)] * (i > 2), ...
            [col_range(1), col_mid] * mod(i, 2) + [col_mid, col_range(2)] * (~mod(i, 2)), ...
            max_points);
    end
end

function multipole_weights = compute_s2m(node, psf, chebyshev_points)
    % S2M：計算葉節點的多極展開
    if node.isLeaf
        [row_indices, col_indices] = ind2sub(size(psf), node.points);
        % 切比雪夫展開
        multipole_weights = zeros(size(chebyshev_points));
        for i = 1:numel(chebyshev_points)
            multipole_weights(i) = sum(psf(row_indices, col_indices) .* chebyshev_interp(node.points, chebyshev_points(i)));
        end
    else
        % M2M：合併子節點的多極展開
        multipole_weights = zeros(size(chebyshev_points));
        for i = 1:4
            multipole_weights = multipole_weights + compute_s2m(node.children{i}, psf, chebyshev_points);
        end
    end
end

function local_expansion = compute_m2l(node, interaction_list, psf)
    local_expansion = zeros(size(psf));
    for other_node = interaction_list
        % 使用多極展開計算遠場影響
        local_expansion = local_expansion + multipole_to_local(node, other_node);
    end
end

function energy_distribution = gather_results(node, psf)
    % 最終收集結果
    if node.isLeaf
        energy_distribution = calculate_leaf_energy(node, psf);
    else
        energy_distribution = zeros(size(psf));
        for i = 1:4
            energy_distribution = energy_distribution + gather_results(node.children{i}, psf);
        end
    end
end


% Step 1: Load data
data = load('test_wo_ans.mat');
PSF = data.kernel; % PSF matrix
layout = data.mask; % layout matrix

% Step 2: Initialize parameters
max_points = 100; % Maximum points per quadtree node
chebyshev_points = linspace(-1, 1, 4); % Chebyshev nodes

% Step 3: Build quadtree and compute FMM
quadtree = build_quadtree(layout, max_points);
multipole_weights = compute_s2m(quadtree, PSF, chebyshev_points);
local_expansion = compute_m2l(quadtree, quadtree.children, PSF);
energy_distribution = gather_results(quadtree, local_expansion);

% Step 4: Visualization
figure;
subplot(2, 2, 1);
imagesc(PSF);
title('PSF');
colorbar;

subplot(2, 2, 2);
imagesc(layout);
title('Layout Mask');
colorbar;

subplot(2, 2, 3);
imagesc(energy_distribution);
title('Energy Distribution (FMM)');
colorbar;

subplot(2, 2, 4);
direct_convolution = conv2(layout, PSF, 'same');
error = abs(direct_convolution - energy_distribution);
imagesc(error);
title('Error Between FMM and Direct Convolution');
colorbar;
