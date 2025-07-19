function [f_handle, obstacles, goal, x0] = generate_maze_squares(n_obstacles, min_width, max_width, eps_obst)
% Generates a speed function with square obstacles on [-1,1]^2.
%
%   n_obstacles : number of square obstacles
%   min_width   : minimum half-width of square obstacles
%   max_width   : maximum half-width of square obstacles
%   eps_obst    : speed inside obstacles (e.g., 1e-3)

% Domain box
L = 1;

% Random centers and widths for obstacles
rng(1);  % For reproducibility
cx = 2*L*rand(1, n_obstacles) - L;
cy = 2*L*rand(1, n_obstacles) - L;
half_widths = min_width + (max_width - min_width) * rand(1, n_obstacles);

% Store obstacle bounding boxes
obstacles = cell(1, n_obstacles);
for i = 1:n_obstacles
  hw = half_widths(i);
  obstacles{i} = [cx(i)-hw, cx(i)+hw;
                  cy(i)-hw, cy(i)+hw];  % [xmin,xmax; ymin,ymax]
end

% Speed function handle
f_handle = @(x,y) evaluate_speed_function(x, y, obstacles, eps_obst);

% Start and goal
goal = [0.9; -0.9];  % bottom-right
x0   = [-0.9; 0.9];  % top-left
end

function f = evaluate_speed_function(x, y, obstacles, eps_obst)
% Evaluate piecewise constant speed function over (x,y) arrays

f = ones(size(x));  % Default speed = 1
for k = 1:length(obstacles)
  box = obstacles{k}; % 2x2: [xmin,xmax; ymin,ymax]
  mask = (x >= box(1,1)) & (x <= box(1,2)) & ...
         (y >= box(2,1)) & (y <= box(2,2));
  f(mask) = eps_obst;
end
end
