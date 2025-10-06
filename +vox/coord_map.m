%create coordinate map

function map = coord_map(voxel_size, chamfer)
    persistent cache_v cache_c cache_map
    if nargin==0
        % Back-compat: pull defaults once
        lat = get_lattice(1);
        voxel_size = lat.voxel_size;
        chamfer    = lat.chamfer;
    end
    
    if ~isempty(cache_map) && cache_v==voxel_size && cache_c==chamfer
        map = cache_map; return;
    end

    lat = get_lattice(1); %gets template parameters
    s = lat.voxel_size/2; %m
    c = lat.chamfer; %m
    % Corners (8 x 3), relative to voxel center
    P = [ -s -s -s;   % 1
           s -s -s;   % 2
           s  s -s;   % 3
          -s  s -s;   % 4
          -s -s  s;   % 5
           s -s  s;   % 6
           s  s  s;   % 7
          -s  s  s ]; % 8

    %face centers
    F = [ 0,-s, 0;  %1
          s, 0, 0;  %2
          0, s, 0;  %3
         -s, 0, 0;  %4
          0, 0, s;  %5
          0, 0,-s ];%6

    % 12 cube edges as pairs of corner indices (12 x 2)
    edge_idx = [ 1 2; 2 3; 3 4; 4 1;   % top square (z = -s)
                 1 5; 2 6; 3 7; 4 8;   % verticals
                 5 6; 6 7; 7 8; 8 5 ]; % bottom square (z = +s)

    % Sample parameters
    tc = c/(2*s);                     % chamfer fraction along edge
    t_vals = [tc, 0.5, 1-tc];         % 1x3

    % Precompute points for each edge and slot: edge_points(e,slot,:) = [x y z]
    edge_points = zeros(12, 3, 3);    % 12 edges x 3 slots x 3 coords 
    for e = 1:12
        p0 = P(edge_idx(e,1), :);     % 1x3
        p1 = P(edge_idx(e,2), :);     % 1x3
        for k = 1:3
            t = t_vals(k);
            edge_points(e, k, :) = (1 - t)*p0 + t*p1;
        end
    end
   
    map.P = P;                        % 8x3 corners (relative)
    map.edge_idx = edge_idx;          % 12x2 indices into P
    map.t_vals = t_vals;              % [tc, 0.5, 1-tc]
    map.edge_points = edge_points;    % 12x3x3 sampled coords (relative)
    map.edge_names = ["E1","E2","E3","E4","E5","E6","E7","E8","E9","E10","E11","E12"];
    map.face_points = F;
    map.chamfer = c;
    map.edge_length = s*2;

    cache_v = voxel_size; cache_c = chamfer; cache_map = map;
    
    % figure; hold on; axis equal
    % grid on
    % xlabel('X'); ylabel('Y'); zlabel('Z');
    % 
    % % Plot cube corners
    % plot3(map.P(:,1), map.P(:,2), map.P(:,3), 'ko', 'MarkerFaceColor','k');
    % 
    % % Plot face centers
    % plot3(map.face_points(:,1), map.face_points(:,2), map.face_points(:,3), ...
    %       'ro', 'MarkerFaceColor','r');
    % 
    % % Plot edge points (flatten out 3rd dim)
    % edge_pts = reshape(permute(map.edge_points,[1 3 2]), [], 3);  % (12*3) × 3
    % plot3(edge_pts(:,1), edge_pts(:,2), edge_pts(:,3), 'bo', 'MarkerFaceColor','b');
end