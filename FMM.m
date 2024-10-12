%% Initialization
function quadtree = build_quadtree(particles, boundary, max_particles)
    % particles: 粒子的位置信息
    % boundary: 節點的邊界 [xmin, xmax, ymin, ymax]
    % max_particles: 每個節點中最多能容納的粒子數
    if numel(particles) > max_particles
        % 若超出最大粒子數，則將區域分成四個子區域
        quadtree = struct('children', {}, 'particles', particles, 'boundary', boundary);
        % 分割區域，並遞迴建立子樹
        quadtree.children{1} = build_quadtree(...); % 子樹1
        quadtree.children{2} = build_quadtree(...)  % 子樹2
        quadtree.children{3} = build_quadtree(...); % 子樹3
        quadtree.children{4} = build_quadtree(...); % 子樹4
    else
        % 若節點中的粒子數小於 max_particles，則不再分割
        quadtree = struct('children', {}, 'particles', particles, 'boundary', boundary);
    end
end

%% Upward pass Downward pass
function multipole_expansion = upward_pass(node)
    if isempty(node.children)
        % 若是葉節點，直接計算多極展開
        multipole_expansion = compute_multipole(node.particles);
    else
        % 若有子節點，遞迴計算子節點的多極展開，並合併到父節點
        for i = 1:4
            child_expansion = upward_pass(node.children{i});
            % 合併多極展開結果
            multipole_expansion = merge_multipole(multipole_expansion, child_expansion);
        end
    end
end



%%  Gather
parfor i = 1:N
    result(i) = compute_local_energy(particles(i), quadtree);
end

 