function goal_pts = generate_goal_region(goal, r_goal, ngoal)
% Generate ngoal + 1 points on a circle centered at goal
% Includes the original goal point at the center

if ngoal == 0
  goal_pts = goal(:);
  return;
end
theta = linspace(0, 2*pi, ngoal + 1);
theta(end) = [];  % avoid duplicating first point

dx = r_goal * cos(theta);
dy = r_goal * sin(theta);

circle_pts = goal(:) + [dx; dy];  % 2 x ngoal
goal_pts = [goal(:), circle_pts]; % 2 x (ngoal + 1)



end
