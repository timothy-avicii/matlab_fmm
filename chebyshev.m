classdef chebyshev
    properties
        N          % Number of Chebyshev nodes
        D          % Number of dimensions
        nodes      % Cached Chebyshev nodes
    end

    methods
        % Constructor
        function obj = chebyshev(N, D)
            % Initialize Chebyshev object
            obj.N = N;
            obj.D = D;
            obj.nodes = obj.generate_nodes(); % Generate nodes during initialization
        end

        % Generate Chebyshev nodes
        function nodes = generate_nodes(obj)
            ms = 0:obj.N-1;
            onedim = cos((ms + 0.5) * pi / obj.N); % 1D Chebyshev nodes
            nodes = obj.cartesian_product(onedim, obj.D); % Cartesian product for D dimensions
        end

        % Similarity function
        % 相似矩陣similar matrix與a與b的大小有關，row是a的row，column是b的row，D維度是a, b的column
        % 若 [a1 a2] = size(a) [b1 b2] = size(b),則 a2 = b2, [s1 s2] = size(s) s1 = a1, s2 = b1
        function S = similarity(obj, a, b)
            % Ensure points are within [-1, 1]
            assert(all(a(:) >= -1 & a(:) <= 1), 'Points in a must be in [-1, 1]');
            assert(all(b(:) >= -1 & b(:) <= 1), 'Points in b must be in [-1, 1]');

            % Dimensions
            [n_a, D] = size(a); % a: (n_a, D)
            [n_b, ~] = size(b); % b: (n_b, D)

            % Precompute k values
            ks = reshape(1:obj.N-1, [1, 1, obj.N-1]); % (1, 1, N-1)

            % Initialize similarity matrix
            S = ones(n_a, n_b);

            % Iterate over each dimension
            for d = 1:D
                % Compute Chebyshev angles
                theta_a = acos(a(:, d)); % (n_a, 1)
                theta_b = acos(b(:, d)); % (n_b, 1)

                % Expand dimensions for broadcasting
                theta_a = reshape(theta_a, [n_a, 1, 1]); % (n_a, 1, 1)
                theta_b = reshape(theta_b, [1, n_b, 1]); % (1, n_b, 1)

                % Compute tensor product terms
                terms = cos(ks .* theta_a) .* cos(ks .* theta_b); % (n_a, n_b, N-1)

                % Sum over Chebyshev polynomial orders
                terms_sum = sum(terms, 3); % Sum over 3rd dimension (N-1)

                % Compute similarity for current dimension
                S_d = (1 / obj.N) + (2 / obj.N) * terms_sum; % (n_a, n_b)

                % Multiply with similarity matrix
                S = S .* S_d; % Element-wise multiplication
            end

            % Squeeze result to remove singleton dimensions
            S = squeeze(S);
        end

        %down_coeffs會是一個四維大小的矩陣 大小為[父節點大小 子節點大小 N D]
        function down_coeffs = downwards_coeffs(obj)

            % Generate child node positions
            shifts = [-0.5, 0.5]; %四個子節點的位置偏移量
            offsets = obj.cartesian_product(shifts, obj.D); % Offset vectors (2^D, D)

            % Compute child nodes by applying offsets to each parent node
            n_parents = size(obj.nodes, 1); % Number of parent nodes
            n_children = size(offsets, 1); % Number of children per parent

            % Preallocate for 4D tensor output
            %down_coeffs = zeros(n_parents, n_children, obj.N, obj.D);
       
            % For each parent node, compute similarity with children
            for i = 1:n_parents
                % Child nodes for the current parent
                child_nodes = obj.nodes(i, :) / 2 + offsets;

                % Compute similarity matrix for this parent
                for d = 1:obj.D
                    down_coeffs(i, :, :, d) = obj.similarity(child_nodes, obj.nodes);
                end
            end
        end

        %up_coeffs會是一個四維大小的矩陣 大小為[子節點大小 父節點大小 N D]
        function up_coeffs = upwards_coeffs(obj)
            % Generate child node positions
            % 子節點與父節點的偏移量
            shifts = [-0.5, 0.5];
            offsets = obj.cartesian_product(shifts, obj.D); % Offset vectors (2^D, D)

            % Compute child nodes by applying offsets to each parent node
            % 用來預分配空間
            n_parents = size(obj.nodes, 1); % Number of parent nodes
            n_children = size(offsets, 1); % Number of children per parent

            % Preallocate for 4D tensor output
            % 這個上行係數矩陣的大小應該要和理論上一樣
            % 但在python的那個程式裡執行出來的大小也是如此，矩陣的entry順序也是正確的
            % 下行係數矩陣的大小也和理論上的大小有出入
            %up_coeffs = zeros(n_children, n_parents, obj.N, obj.D);

            % For each child node set, compute similarity with parents
            for i = 1:n_children
                % Parent nodes for the current child set
                child_nodes = obj.nodes / 2 + offsets(i, :);

                % Compute similarity matrix for this child
                for d = 1:obj.D
                    up_coeffs(i, :, :, d) = obj.similarity(obj.nodes, child_nodes);
                end
            end
        end

        % Interpolation weights for arbitrary points
        function values = interpolate(obj, points, weights)
            % points: Interpolation points (n_points, D)
            % weights: Interpolation weights (n_nodes, 1)

            % Compute similarity between points and nodes
            % 用相似矩陣去做插值
            S = obj.similarity(points, obj.nodes); % (n_points, n_nodes)

            % Interpolate values
            values = S .* weights; % (n_points, 1)
        end

        % Compute interpolation weights from known values
        function weights = anterpolate(obj, points, values)
            % points: Known points (n_points, D)
            % values: Known values (n_points, 1)

            % Compute similarity between nodes and points
            S = obj.similarity(obj.nodes, points); % (n_nodes, n_points)
            % 用相似矩陣去做反插值

            % Compute weights
            weights = S .* values; % (n_nodes, 1)
        end


    end

    methods (Static)
        % Cartesian product helper function
        function result = cartesian_product(vec, D)
            grids = cell(1, D);
            [grids{:}] = ndgrid(vec);
            result = reshape(cat(D+1, grids{:}), [], D);
        end
    end
end
