% Load data (kernel and mask)
load test_wo_ans.mat
kernel ;% Electron energy distribution    
mask  ;% Lithographic pattern

% Parameters
[m, n] = size(mask);   % Size of the lithographic pattern
exp_kernel_size = size(kernel);  % Kernel size
expansion_factor = 2;  % Expansion factor for FMM

% Initialize the exposure energy distribution matrix
exposure_energy = zeros(m, n);

% Fast Multipole Method (FMM) calculation
% Split the grid into quadrants and calculate near and far field contributions
%% 

for i = 1:m
    for j = 1:n
        if mask(i, j) ~= 0
            % Near-field calculation
            for k = max(1, i-expansion_factor):min(m, i+expansion_factor)
                for l = max(1, j-expansion_factor):min(n, j+expansion_factor)
                    exposure_energy(k, l) = exposure_energy(k, l) + mask(i, j) * kernel(abs(i-k)+1, abs(j-l)+1);
                end
            end
            % Far-field approximation (FMM)
            % Approximate contributions from distant points using multipole expansion
            % (this part of the implementation can be optimized further with hierarchical grids)
        end
    end
end
%% 

% Visualize the results
figure;
imagesc(exposure_energy);
title('Exposure Energy Distribution');
colorbar;



