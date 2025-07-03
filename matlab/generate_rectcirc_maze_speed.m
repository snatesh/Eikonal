% Generate a speed function for a maze within [-1,1]^2 using disks to fill maze walls
% Each wall cell is replaced by a circular obstacle with speed eps_obst

function [f_handle, obstacles, goal, x0] = generate_rectcirc_maze_speed(M, eps_obst)
% M         : number of cells along one side (odd number preferred)
% eps_obst  : speed value inside maze wall disks

% Step 1: Generate maze
maze = generate_maze_recursive_backtracker(M);

% Step 2: Map maze to geometry
xmin = -1; xmax = 1;
ymin = -1; ymax = 1;
dx = (xmax - xmin) / M;
dy = (ymax - ymin) / M;

% Step 3: Extract wall disk centers and radii
centers = [];
radius = min(dx, dy) * 0.45;  % wall disk radius (slightly smaller than cell)
for i = 1:M
    for j = 1:M
        if maze(i,j) == 0  % wall cell
            xc = xmin + (j - 0.5)*dx;
            yc = ymax - (i - 0.5)*dy;
            centers = [centers [xc; yc]];
        end
    end
end

r0s = radius * ones(1, size(centers, 2));
obstacles = [centers; r0s];  % 3 x N: [x; y; r]

% Step 4: Define speed function
f_handle = @(x,y) maze_speed_disks(x, y, centers, r0s, eps_obst);

% Step 5: Entry and Exit (top-left and bottom-right path cells)
[entry_row, entry_col] = deal(1,2);  % second cell on top row
[exit_row, exit_col] = deal(M,M-1); % second-last cell on bottom row

goal = [xmin + (entry_col-0.5)*dx; ymax - (entry_row-0.5)*dy];
x0   = [xmin + (exit_col -0.5)*dx; ymin + (M-exit_row+0.5)*dy];

% Ensure goal and x0 are not inside any disk
for k = 1:size(centers,2)
    if norm(goal - centers(:,k)) < r0s(k)
        error('Goal point is inside an obstacle disk.');
    end
    if norm(x0 - centers(:,k)) < r0s(k)
        error('Start point is inside an obstacle disk.');
    end
end

end

function f = maze_speed_disks(x, y, centers, r0s, eps_obst)
% Piecewise speed function: eps_obst inside disks, 1 outside
x = x(:); y = y(:);
XY = [x y]';
N = length(x);
f = ones(N, 1);

for i = 1:size(centers, 2)
    c = centers(:,i);
    r = r0s(i);
    inside = sum((XY - c).^2, 1) < r^2;
    f(inside') = eps_obst;
end
end

function maze = generate_maze_recursive_backtracker(M)
% Recursive backtracker maze generation (binary matrix)
maze = zeros(M, M);
visited = false(M, M);

% Directions: [drow, dcol]
dirs = [0 2; 2 0; 0 -2; -2 0];

% Ensure M is odd
if mod(M,2) == 0
    error('M must be odd');
end

% Initialize: mark all walls, then start at (1,1)
maze(1:2:end,1:2:end) = 1;  % mark potential path cells
stack = [1, 1];
visited(1,1) = true;

while ~isempty(stack)
    [r,c] = deal(stack(end,1), stack(end,2));
    neighbors = [];
    for d = 1:4
        nr = r + dirs(d,1);
        nc = c + dirs(d,2);
        if nr >= 1 && nr <= M && nc >= 1 && nc <= M && ~visited(nr,nc)
            neighbors = [neighbors; nr nc dirs(d,:)];
        end
    end
    if ~isempty(neighbors)
        idx = randi(size(neighbors,1));
        nr = neighbors(idx,1); nc = neighbors(idx,2);
        dr = neighbors(idx,3); dc = neighbors(idx,4);
        maze(r + dr/2, c + dc/2) = 1;  % open wall between
        visited(nr,nc) = true;
        stack = [stack; nr, nc];
    else
        stack(end,:) = [];  % backtrack
    end
end
end

