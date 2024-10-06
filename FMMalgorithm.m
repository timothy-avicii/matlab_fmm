%% Initialization
function quadtree = initialize_quadtree(particles, boundary, max_particles)
    % particles: 粒子的位置信息
    % boundary: 節點的邊界 [xmin, xmax, ymin, ymax]
    % max_particles: 每個節點中最多能容納的粒子數
    
    % 如果當前節點的粒子數量超過 max_particles，則分割區域
    if numel(particles) > max_particles
        quadtree = struct('children', {}, 'particles', [], 'boundary', boundary);
        % 分割成四個區域
        [region1, region2, region3, region4] = split_boundary(boundary);
        % 根據位置將粒子分配到各個區域
        particles1 = filter_particles(particles, region1);
        particles2 = filter_particles(particles, region2);
        particles3 = filter_particles(particles, region3);
        particles4 = filter_particles(particles, region4);
        
        % 遞迴地為每個子區域建立四叉樹
        quadtree.children{1} = initialize_quadtree(particles1, region1, max_particles);
        quadtree.children{2} = initialize_quadtree(particles2, region2, max_particles);
        quadtree.children{3} = initialize_quadtree(particles3, region3, max_particles);
        quadtree.children{4} = initialize_quadtree(particles4, region4, max_particles);
    else
        % 如果粒子數量在可接受範圍內，則不分割
        quadtree = struct('children', {}, 'particles', particles, 'boundary', boundary);
    end
end

function [region1, region2, region3, region4] = split_boundary(boundary)
    % 根據當前邊界將其分成四個子區域
    mid_x = (boundary(1) + boundary(2)) / 2;
    mid_y = (boundary(3) + boundary(4)) / 2;
    
    region1 = [boundary(1), mid_x, boundary(3), mid_y];
    region2 = [mid_x, boundary(2), boundary(3), mid_y];
    region3 = [boundary(1), mid_x, mid_y, boundary(4)];
    region4 = [mid_x, boundary(2), mid_y, boundary(4)];
end

function filtered_particles = filter_particles(particles, region)
    % 過濾出位於 region 內的粒子
    filtered_particles = particles(particles(:,1) >= region(1) & particles(:,1) <= region(2) & ...
                                   particles(:,2) >= region(3) & particles(:,2) <= region(4), :);
end

%% Upward pass
function multipole_expansion = upward_pass(node)
    % 如果是葉節點，計算多極展開
    if isempty(node.children)
        multipole_expansion = compute_multipole(node.particles);
    else
        % 否則，遞迴地計算每個子節點的多極展開，並合併結果
        multipole_expansion = 0;
        for i = 1:4
            child_expansion = upward_pass(node.children{i});
            multipole_expansion = merge_multipole(multipole_expansion, child_expansion);
        end
    end
end

function multipole = compute_multipole(particles)
    % 根據粒子的位置和劑量值計算多極展開
    % 此處可自訂 Chebyshev 插值或其他方法
    multipole = sum(particles(:, 3));  % 假設第三列是劑量值
end

function result = merge_multipole(parent_multipole, child_multipole)
    % 合併父節點和子節點的多極展開
    result = parent_multipole + child_multipole;
end


%% Downward pass
function local_expansion = downward_pass(node, parent_local_expansion)
    % 計算當前節點的局部展開，從父節點向下傳遞
    if isempty(node.children)
        local_expansion = parent_local_expansion + compute_local(node.particles);
    else
        % 遞迴地對每個子節點進行局部展開的傳遞
        for i = 1:4
            node.children{i}.local_expansion = downward_pass(node.children{i}, parent_local_expansion);
        end
    end
end

function local = compute_local(particles)
    % 計算粒子的局部展開
    local = sum(particles(:, 3));  % 假設第三列是劑量值
end


%% Gather
function energy_distribution = gather(node)
    if isempty(node.children)
        % 如果是葉節點，返回能量分佈
        energy_distribution = compute_energy(node.particles);
    else
        % 否則，遞迴地收集每個子節點的能量分佈
        energy_distribution = 0;
        for i = 1:4
            energy_distribution = energy_distribution + gather(node.children{i});
        end
    end
end

function energy = compute_energy(particles)
    % 計算每個粒子的最終能量分佈
    % 假設能量與多極展開和局部展開的結果相關
    energy = sum(particles(:, 3));  % 假設第三列是劑量值
end

