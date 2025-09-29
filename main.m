%% Stiffness Matrix
clear;clc;
A =         [1,1,1
            ];
A(:,:,2) =  [1,1,1
            ];
% A = 1;

% vox_arr = init_arr(A);

% % get_nodes()
% conn = get_conn_map(A);
% conn.x, conn.z
% get_voxel_config(get_nodes(), 2);
lat = get_lattice(A);
cfg = config_array(lat);
stiffness_matrix = compute_matrix(cfg);
n = define_nodes(1);
%% code for displacement BC
% use sensor readings to find global node displacements by assuming top
% face only moves in z rx ry and bottom face is fixed.
% estimate_displacements(stiffness_matrix, [n.f1.tm, n.f2.tm, n.f3.tm], [n.f1.bm, n.f2.bm, n.f3.bm], [250e-6, 260e-6, 240e-6])

sensor_beams = [n.f1.c, n.f1.bm
                n.f2.c, n.f2.bm
                n.f3.c, n.f3.bm];

sensor_readings = -[250e-6, 250e-6, 250e-6];
u = estimate_displacements(stiffness_matrix,sensor_beams(:,1), sensor_beams(:,2), sensor_readings);
coords = beamCoords(cfg);
plotBeamDeform(cfg, coords, u, 1000, sensor_beams);
%% template code for force BC

% n2 = define_nodes(2);

%Define applied Force
F_value = 1000; %Define load condition: 1 Newton in y1 direction
F_input = zeros(length(stiffness_matrix)/6,6); %index using F(node, dofs)
top_nodes = [get_face_nodes(n1,"f5",1),];
% top_nodes = top_nodes(1:end-1);
for node = top_nodes
    F_input(node,1:3) = [0,0,-F_value];
end

%Define constraints
constrained_nodes = [get_face_nodes(n1, "f6"),];
Displacements = find_displacements(stiffness_matrix, F_input, constrained_nodes);

% % -- after you have beam_config and Displacements -----------------
coords = beamCoords(cfg);
plotBeamDeform(cfg, coords, Displacements, 20);  % scale=50×
% plotBeamDeform(cfg, coords, zeros(size(coords,1),6), 1); % no deformation
% showBeamAxes(cfg, coords, 0.12); 


%find displacements from force bc
function x_full = find_displacements(K_actual, F, constrained_nodes)
    % finds numerical displacement values
    % Inputs:
    % k_assembly - 3D Stiffness matrix
    % loads - Force vector with applied loads
    % constraints - fixed degrees of freedom
    % Outputs: 
    % 6 by n matrix where each column is a node
    % index using (node, dof)
    F = F(:);
    C = [];
    for node = constrained_nodes
        C = [C, 6*node-5:6*node]; % x, y, theta DOFs
    end

    % Define free DOFs (all DOFs minus constrained DOFs)
    all_dofs = 1:size(K_actual, 1);
    free_dofs = setdiff(all_dofs, C);

    % Reduce stiffness matrix and force vector
    K_reduced = K_actual(free_dofs, free_dofs); % Remove constrained rows/columns
    F_reduced = F(free_dofs); % Remove constrained rows
    x_reduced = K_reduced \ F_reduced; % Solve for displacements at free DOFs
    
    % Reconstruct the full displacement vector
    x_full = zeros(size(K_actual, 1), 1);
    x_full(free_dofs) = x_reduced; % Fill in free DOFs (constrained DOFs remain zero)
    x_full = reshape(x_full, 6, [])';
end

% Display beam_config
function nodeCoords = beamCoords(beam_config)
% Returns an array  [x y z]  for every node index that appears.
%
% beam_config cols: [node1 node2 anglex angley anglez L]

    ex = [1;0;0];                           % local x unit
    maxNode = max( beam_config(:,1:2) ,[],'all');
    nodeCoords = nan(maxNode,3);            % pre‑fill NaNs
    nodeCoords(1,:) = get_xyz(1);              % start at origin

    unresolved = true;
    while unresolved
        unresolved = false;
        for i = 1:size(beam_config,1)
            n1 = beam_config(i,1);      n2 = beam_config(i,2);
            ax = beam_config(i,3:5);      ay = beam_config(i,6:8);
            az = beam_config(i,9:11);      L  = beam_config(i,12);

            % rotation matrix (intrinsic X‑Y‑Z, matches your code)
            R  = [ax' ay' az'];
            vec = (R*ex*L).';           % 1×3

            if  ~isnan(nodeCoords(n1,1)) && isnan(nodeCoords(n2,1))
                nodeCoords(n2,:) = nodeCoords(n1,:) + vec;
                unresolved = true;
            elseif ~isnan(nodeCoords(n2,1)) && isnan(nodeCoords(n1,1))
                nodeCoords(n1,:) = nodeCoords(n2,:) - vec;
                unresolved = true;
            end
        end
    end
end

function plotBeamDeform(beam_config, nodeCoords, U, scale, sensor_id)
% Plots undeformed (gray) and deformed (blue) beams.
% U : N×6 displacement matrix, translational DOFs in columns 1..3

    if nargin < 4, scale = 1; end
    if nargin < 5, sensor_id = [];end
    U = U(:,1:3) * scale;          % use only translations

    figure; hold on; grid on; axis equal
    view(3);                       % 3‑D view
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('Undeformed (gray) and deformed (blue)');

    for i = 1:size(beam_config,1)
        n1 = beam_config(i,1); n2 = beam_config(i,2);

        P1 = nodeCoords(n1,:);   P2 = nodeCoords(n2,:);
        % undeformed
        plot3([P1(1) P2(1)], [P1(2) P2(2)], [P1(3) P2(3)], 'k-', ...
              'Color',[0.6 0.6 0.6],'LineWidth',1);

        % deformed
        P1d = P1 + U(n1,:);  P2d = P2 + U(n2,:);
        plot3([P1d(1) P2d(1)], [P1d(2) P2d(2)], [P1d(3) P2d(3)], ...
              'b-','LineWidth',1.5);
    end

    % node labels
    for n = 1:size(nodeCoords,1)
        p = nodeCoords(n,:);
        text(p(1),p(2),p(3),sprintf(' %d',n), 'FontSize',8);
    end

    %add sensors
    if ~isempty(sensor_id)
        sensor_coords = arrayfun(@get_xyz,sensor_id,'UniformOutput',false);
        for i = 1:size(sensor_coords,1)
            p1 = sensor_coords{i,1} + U(sensor_id(i,1),1:3);  % first endpoint (row vector [x y z])
            p2 = sensor_coords{i,2} + U(sensor_id(i,2),1:3);  % second endpoint
            plot3([p1(1) p2(1)], [p1(2) p2(2)], [p1(3) p2(3)], 'r-', 'LineWidth', 2);
        end
    end
end

function showBeamAxes(beam_config, nodeCoords, scale)
% Draws a short green arrow = local +Y
%        and a short red   arrow = local +Z
% at the mid‑point of every beam.
% scale : arrow length (default = 10 % of beam length)

    if nargin<3, scale = 0.10; end
    hold on
    for i = 1:size(beam_config,1)
        n1 = beam_config(i,1);   n2 = beam_config(i,2);
        ax = beam_config(i,3);   ay = beam_config(i,4);   az = beam_config(i,5);
        L  = beam_config(i,6);

        P1 = nodeCoords(n1,:);   P2 = nodeCoords(n2,:);
        P0 = (P1+P2)/2;                       % mid‑point

        % build the same rotation used in your stiffness
        Rx = [1 0 0; 0 cos(ax) -sin(ax); 0 sin(ax) cos(ax)];
        Ry = [cos(ay) 0 sin(ay); 0 1 0; -sin(ay) 0 cos(ay)];
        Rz = [cos(az) -sin(az) 0; sin(az) cos(az) 0; 0 0 1];
        R  = Rx*Ry*Rz;                        % extrinsic X‑Y‑Z

        ey = R*[0;1;0];                       % local +y
        ez = R*[0;0;1];                       % local +z

        aLen = scale*L;
        quiver3(P0(1),P0(2),P0(3), ey(1),ey(2),ey(3), aLen,...
                'Color',[0 0.7 0],'LineWidth',1,'MaxHeadSize',0.5);
        quiver3(P0(1),P0(2),P0(3), ez(1),ez(2),ez(3), aLen,...
                'Color',[0.8 0 0],'LineWidth',1,'MaxHeadSize',0.5);
    end
    legend({'undeformed','deformed','local +Y','local +Z'},'Location','best');
end

function nodes = get_face_nodes(n, face, exclude_face)
    %returns all nodes ids of a face
    if nargin<3
        exclude_face = 0;
    end

    nodes = [];
    for name = n.face_nodes
        nodes = [nodes, n.(face).(name)];
    end
    if exclude_face
        nodes = nodes(1:end-1);
    end
end
