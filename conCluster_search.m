function [maxClusterValue, maxClusterIndices, maxClusterCoords] = conCluster_search(input_coords, input_voxels, cluster_size, connectivity, average_method, plot_mode)
%%CONCLUSTER_SEARCH Find the strongest contiguous cluster of a specified size among 3D voxels
% input_coords: [N x 3] integer coords
% input_voxels: [N x 1] values aligned with input_coords
% cluster_size: scalar (e.g., 5)
% connectivity: 6, 18, or 26
% average_method: function handle, e.g. @mean, @sum
% plot_mode: visualize the cluster ('3d', 'patch', or false)
%
% Outputs:
%   maxClusterValue   Scalar. The score of the selected cluster computed as:
%                     average_method(values_of_selected_voxels).
%                     Example: if average_method=@mean, this is the mean value
%                     across the returned cluster voxels.
%
%   maxClusterIndices [cluster_size x 1]. Indices into the original input arrays
%                     (input_coords / input_voxels) identifying which input voxels
%                     belong to the selected cluster.
%
%   maxClusterCoords  [cluster_size x 3]. The coordinates (in the same original
%                     coordinate system as input_coords) of the selected cluster
%                     voxels, ordered consistently with maxClusterIndices.

if ~exist('cluster_size', 'var') || isempty(cluster_size)
    cluster_size = 5;
end
if ~exist('connectivity', 'var') || isempty(connectivity)
    connectivity = 6;
else
    assert(ismember(connectivity, [6, 18, 26]), 'Incorrect input: 3D connectivity must be 6, 18, or 26.')
end
if ~exist('average_method', 'var') || isempty(average_method)
    average_method = @mean;
end
if ~exist('plot_mode', 'var') || isempty(plot_mode)
    plot_mode = false;
elseif ~isequal(plot_mode, false)
    assert(ismember(plot_mode, {'3d', 'patch'}), "Incorrect input: specified plotting mode must be '3d' or 'patch'.")
end

%% ---- 1) Build a local bounding box & local coords ----
mins = min(input_coords, [], 1);
maxs = max(input_coords, [], 1);
boxSizeLocal = maxs - mins + 1;
N_voxels = numel(input_voxels);
assert(N_voxels >= cluster_size, 'Incorrect input: cluster_size exceeds number of input voxels.')

coordsL = input_coords - mins + 1; % local coords starting from 1

% Map each input voxel -> local linear index
linAll = sub2ind(boxSizeLocal, coordsL(:,1), coordsL(:,2), coordsL(:,3));

% Build a value volume (NaN where no voxel exists)
valVol = nan(boxSizeLocal);
valVol(linAll) = input_voxels;

% Create a linear indexing offset for connected voxels depending on connectivity
connOffsets = makeConnOffsets(size(valVol), connectivity);

% Mapping from local linear index -> input index (lets us return input indices)
lin2input = zeros(prod(boxSizeLocal), 1, 'uint32');
lin2input(linAll) = uint32(1:N_voxels);

%% ---- 2) Threshold search ----
sorted_voxels = sort(input_voxels, 'descend');
thresh_idx = cluster_size; % initialize index
found = false;

while ~found && thresh_idx <= N_voxels
    thresh = sorted_voxels(thresh_idx);

    % suprathreshold local linear indices
    above_lin = linAll(input_voxels >= thresh);

    mask = false(boxSizeLocal);
    mask(above_lin) = true;

    CC = bwconncomp(mask, connectivity);
    c_sizes = cellfun(@numel, CC.PixelIdxList);
    keep = find(c_sizes >= cluster_size);

    if ~isempty(keep)
        found = true;
    else
        thresh_idx = thresh_idx + 1;
    end
end

assert(found, ['Cannot find any cluster of size ', num2str(cluster_size), '. Please double check inputs!' ])

%% ---- 3) Process each candidate cluster ----
N_clusters = numel(keep);
candidate_cluster_val = nan(1, N_clusters);
candidate_cluster_lin = cell(1, N_clusters);

for ii = 1:N_clusters
    cl_lin = CC.PixelIdxList{keep(ii)};
    cl_vals = valVol(cl_lin);

    % pick the top cluster_size voxels that preserve connectivity
    if numel(cl_lin) > cluster_size
        idx_local = shrinkCluster(cl_lin, cl_vals, connOffsets, cluster_size, average_method);
        cl_lin  = cl_lin(idx_local);
        cl_vals = cl_vals(idx_local);
    end

    candidate_cluster_val(ii) = average_method(cl_vals);
    candidate_cluster_lin{ii} = cl_lin;
end

%% ---- 4) Max cluster ----
[maxClusterValue, bestIdx] = max(candidate_cluster_val);
best_lin = candidate_cluster_lin{bestIdx};
maxClusterIndices = lin2input(best_lin); % Return indices into the original input list
maxClusterCoords = input_coords(maxClusterIndices,:);

%% ---- (Optional) Visualize the cluster ----
if ~plot_mode
    return
else
    switch plot_mode
        case '3d'
            clusterMask = false(size(valVol));
            clusterMask(best_lin) = true;
            valVol_plot = valVol;
            valVol_plot(~isfinite(valVol_plot)) = min(input_voxels);

            % --- Swap X<->Y so volume display matches patch coords (X,Y,Z) ---
            valVol_plot  = permute(valVol_plot,  [2 1 3]); % (Y,X,Z) -> (X,Y,Z) for display
            clusterMask  = permute(clusterMask,  [2 1 3]);

            viewer = viewer3d;
            viewer.BackgroundColor = [1 1 1]; % white

            % --- Background volume ---
            hVol = volshow(valVol_plot, ...
                'Parent', viewer, ...
                'Colormap', gray(256));

            % Make background semi-transparent
            hVol.Alphamap = linspace(0, 0.05, 256);

            % --- Cluster overlay ---
            hClust = volshow(clusterMask, ...
                'Parent', viewer, ...
                'Colormap', [1 0 0]);

            % Fully opaque cluster
            hClust.Alphamap = [0; 1];

        case 'patch'
            figure
            ax = gca;
            hold on

            patchCubes(input_coords, input_voxels, ...
                'FaceAlpha', 0.02, 'EdgeColor', 'none');

            % grayscale mapping for background
            colormap(ax, gray);
            cb = colorbar(ax);
            cb.Label.String = 'Voxel value';

            % Set color limits
            clim(ax, prctile(input_voxels, [2 98]));

            % red cubes for identified cluster
            patchCubes(maxClusterCoords, [], ...
                'FaceColor', [1 0 0], ... % red
                'FaceAlpha', 1.0, 'EdgeColor', 'none');

            % viewing / lighting
            axis(ax, 'equal');
            view(ax, 3);
            camlight(ax, 'headlight');
            lighting(ax, 'gouraud');
            grid(ax, 'on');

            % labels and formatting
            xlabel('X')
            ylabel('Y')
            zlabel('Z')
            set(ax, 'FontSize', 20)
    end
end
end


%% HELPER FUNCTIONS
function connOffsets = makeConnOffsets(volSize, connectivity)
nx = volSize(1);
ny = volSize(2);

[i,j,k] = ndgrid(-1:1, -1:1, -1:1);
d = [i(:) j(:) k(:)];
d(all(d==0,2),:) = [];

switch connectivity
    case 6
        d = d(sum(abs(d),2)==1,:);
    case 18
        d = d(sum(abs(d),2)<=2,:);
    case 26
        % keep all
    otherwise
        error('connectivity must be 6, 18, or 26');
end

% MATLAB linear index: +1 in i is +1; +1 in j is +nx; +1 in k is +nx*ny
connOffsets = d(:,1) + d(:,2)*nx + d(:,3)*nx*ny;
end


function idx_best = shrinkCluster(cl_lin, cl_vals, connOffsets, cluster_size, average_method)
%%SHRINKCLUSTER Shrink a subcluster to the specified size while maintaining connectivity
% cl_lin: Mx1 global linear indices (one connected component)
% cl_vals: Mx1 values aligned with cl_lin
% connOffsets: Px1 linear neighbor offsets for the chosen connectivity
% cluster_size: specified cluster size

M = numel(cl_lin);
cluster_size = min(cluster_size, M);
if cluster_size <= 0
    idx_best = [];
    return;
end

% Map global linear idx -> local (1..M) for O(1) membership
mp = containers.Map('KeyType','double','ValueType','int32');
for i = 1:M
    mp(cl_lin(i)) = int32(i);
end

% pick top seeds by value
[~, seeds] = sort(cl_vals, 'descend');

bestScore = -Inf;
idx_best = [];

for s = seeds(:)'
    chosen = false(M,1);
    chosen(s) = true;
    set = s;

    % grow connected set until size K
    while numel(set) < cluster_size
        % frontier = all neighbors of current set not already chosen
        frontier = [];

        for u = set(:)'
            base = cl_lin(u);
            cand = base + connOffsets(:);

            % collect candidates that are inside this component
            for t = 1:numel(cand)
                if isKey(mp, cand(t))
                    j = mp(cand(t));
                    if ~chosen(j)
                        frontier(end+1,1) = j; %#ok<AGROW>
                    end
                end
            end
        end

        if isempty(frontier)
            break; % can't expand from this seed to size K
        end

        frontier = unique(frontier);
        [~, mx] = max(cl_vals(frontier));
        j = frontier(mx);

        chosen(j) = true;
        set(end+1,1) = j; %#ok<AGROW>
    end

    if numel(set) == cluster_size
        sc = average_method(cl_vals(set));
        if sc > bestScore
            bestScore = sc;
            idx_best = set;
        end
    end
end
end


function p = patchCubes(coords, vals, varargin)
%PATCHCUBES Draw voxel cubes at integer coords using a single patch object.
% coords: N x 3 voxel indices [i j k]
% vals  : N x 1 values for colormap scaling; pass [] for solid color

coords = double(coords);
N = size(coords, 1);

% ----- define unit cube (8 vertices, 6 faces) -----
V0 = [ ...
    0 0 0;
    1 0 0;
    1 1 0;
    0 1 0;
    0 0 1;
    1 0 1;
    1 1 1;
    0 1 1];

F0 = [ ...
    1 2 3 4;
    5 6 7 8;
    1 2 6 5;
    2 3 7 6;
    3 4 8 7;
    4 1 5 8];

% center cube on voxel index
V0 = V0 - 0.5;

% ----- build vertices/faces for all cubes -----
V = zeros(8*N, 3);
F = zeros(6*N, 4);

for i = 1:N
    vIdx = (i-1)*8 + (1:8);
    fIdx = (i-1)*6 + (1:6);

    V(vIdx,:) = V0 + coords(i,:);
    F(fIdx,:) = F0 + (i-1)*8;
end

% ----- create patch -----
if ~isempty(vals)
    vals = double(vals(:));
    if numel(vals) ~= N
        error('vals must be N x 1 to match coords.');
    end

    % For FaceColor='flat', patch expects one CData per face.
    % Each voxel has 6 faces, so repeat each value 6 times.
    C = repelem(vals, 6);

    p = patch('Vertices', V, ...
        'Faces', F, ...
        'FaceColor', 'flat', ...
        'FaceVertexCData', C, ...
        varargin{:});
else
    p = patch('Vertices', V, ...
        'Faces', F, ...
        varargin{:});
end
end
