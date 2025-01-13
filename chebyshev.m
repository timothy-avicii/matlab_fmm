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
        function S = similarity(obj, a, b)
            % Ensure points are within [-1, 1]
            assert(all(a(:) >= -1 & a(:) <= 1), 'Points in a must be in [-1, 1]');
            assert(all(b(:) >= -1 & b(:) <= 1), 'Points in b must be in [-1, 1]');

            % Compute Chebyshev angles
            theta_a = acos(a); % Angles for a
            theta_b = acos(b); % Angles for b

            % Precompute k values
            ks = (1:obj.N-1); % Chebyshev polynomial orders (column vector)

            % Expand dimensions for broadcasting
            theta_a = reshape(theta_a, [size(a, 1), 1]); % Ensure theta_a is column
            theta_b = reshape(theta_b, [size(b, 1), 1]);% Ensure theta_b is row
            
            % Tensor product terms
            terms = cos(theta_a * ks) .* cos(theta_b * ks) % Element-wise product
           
            % Compute similarity matrix
            S = prod((1/obj.N) + (2/obj.N) * sum(terms, 1), 2);
        end


        % Interpolation weights for downwards propagation
        function down_coeffs = downwards_coeffs(obj)
            shifts = [-0.5, 0.5];
            child_nodes = obj.cartesian_product(shifts, obj.D) + obj.nodes / 2; % Child nodes
            down_coeffs = obj.similarity(child_nodes, obj.nodes);
        end

        % Interpolation weights for upwards propagation
        function up_coeffs = upwards_coeffs(obj)
            shifts = [-0.5, 0.5];
            child_nodes = obj.cartesian_product(shifts, obj.D) + obj.nodes / 2; % Child nodes
            up_coeffs = obj.similarity(obj.nodes, child_nodes);
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
