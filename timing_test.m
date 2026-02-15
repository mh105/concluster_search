%% Timing test for conCluster_search()
close all
clear all
clc

boxSize = [25, 40, 35];
lc_voxels = rand(boxSize);
[x, y, z] = ndgrid(1:boxSize(1), 1:boxSize(2), 1:boxSize(3));

% Linearize the voxels and coordinates
input_voxels = lc_voxels(:);
input_coords = [x(:), y(:), z(:)];

cluster_sizes = 5:20;
connectivities = [6, 18, 26];
average_method = @mean;
plot_on = false;

% Preallocate timing results
T = nan(length(cluster_sizes), length(connectivities));

% --- Warm-up (important) ---
conCluster_search( ...
    input_coords, input_voxels, ...
    cluster_sizes(1), connectivities(1), ...
    average_method, plot_on);

% --- Timing loop ---
for c = 1:length(connectivities)
    conn = connectivities(c);

    for k = 1:length(cluster_sizes)
        cs = cluster_sizes(k);
        fprintf('Timing conCluster_search | connectivity = %d | cluster_size = %d\n', conn, cs)
        
        t = tic;
        conCluster_search(input_coords, input_voxels, cs, conn, average_method, plot_on);
        T(k, c) = toc(t);
    end
end

% --- Plot ---
figure; hold on
plot(cluster_sizes, T(:,1), '-o', 'LineWidth', 2)
plot(cluster_sizes, T(:,2), '-s', 'LineWidth', 2)
plot(cluster_sizes, T(:,3), '-^', 'LineWidth', 2)

xlabel('Cluster size')
ylabel('Execution time (seconds)')
title('conCluster\_search runtime scaling')

legend( ...
    '6-connected', ...
    '18-connected', ...
    '26-connected', ...
    'Location', 'northwest')

set(gca, 'FontSize', 20)
grid on
box on
