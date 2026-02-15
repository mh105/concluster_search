close all
clear all
clc

%% Simulate a box of voxels
boxSize = [25, 40, 35];
lc_voxels = rand(boxSize);
[x, y, z] = ndgrid(1:boxSize(1), 1:boxSize(2), 1:boxSize(3));

% Linearize the voxels and coordinates
input_voxels = lc_voxels(:);
input_coords = [x(:), y(:), z(:)];

%% Contiguous cluster search algorithm - example use
cluster_size = 5;
connectivity = 6; % 6, 18, 26
average_method = @mean; % @mean, @median, @sum, @max, @min, etc.
plot_on = false; % '3d', 'patch', or false to turn off plotting

t = tic;
[maxClusterValue, maxClusterIndices, maxClusterCoords] = conCluster_search(...
    input_coords, input_voxels, cluster_size, connectivity, average_method, plot_on);
toc(t)
