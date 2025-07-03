% Generate a speed function for a maze within [-1,1]^2
% using recursive backtracker and map it to thin walls

function [f_handle, obstacles, goal, x0] = generate_maze_speed(M, eps_obst)
% M         : number of cells along one side (odd number preferred)
% eps_obst  : speed value inside maze walls

% Step 1: Generate maze
maze = generate_maze_recursive_backtracker(M);

% Step 2: Map maze to geometry
xmin = -1; xmax = 1;
ymin = -1; ymax = 1;
dx = (xmax - xmin) / M;
dy = (ymax - ymin) / M;

obstacles = {};  % list of [xlo xhi ylo yhi] wall rectangles
for i = 1:M
    for j = 1:M
        if maze(i,j) == 0  % wall
            xlo = xmin + (j-1)*dx;
            xhi = xlo + dx;
            ylo = ymin + (M-i)*dy;
            yhi = ylo + dy;
            obstacles{end+1} = [xlo xhi ylo yhi];
        end
    end
end

% Step 3: Define speed function handle
f_handle = @(x,y) speed_fn_maze(x, y, obstacles, eps_obst);

% Step 4: Entry and Exit (top-left and bottom-right path cells)
[entry_row, entry_col] = deal(1,2);  % second cell on top row
[exit_row, exit_col] = deal(M,M-1); % second-last cell on bottom row

% Find cell centers and ensure they are in path (maze=1)
goal = [xmin + (entry_col-0.5)*dx; ymax - (entry_row-0.5)*dy];
x0   = [xmin + (exit_col -0.5)*dx; ymin + (M-exit_row+0.5)*dy];

if maze(entry_row, entry_col) == 0
    error('Goal is inside a wall cell. Adjust maze parameters.');
end
if maze(exit_row, exit_col) == 0
    error('Start point is inside a wall cell. Adjust maze parameters.');
end

end

function f = speed_fn_maze(x, y, obstacles, eps_obst)
% Evaluate speed function at arrays x, y
f = ones(size(x));
for k = 1:length(obstacles)
    rect = obstacles{k};
    in_rect = x >= rect(1) & x <= rect(2) & ...
              y >= rect(3) & y <= rect(4);
    f(in_rect) = eps_obst;
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

