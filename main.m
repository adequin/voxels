%remaining issues
% 1. Bending inertias of stub beams. currently they are implemented as two
%rectangular beams that have a cross section of +, but probably should be
%offset by setting individual inertias offset using parallel axis.
% 2. Length parameters are still a bit arbitrary and only approximate all
% the fillets and non-uniform geometry
% 3. Better UI for constraining dofs instead of doing 6*node-5 or smth
% 4. Sensor placement and interpretation to for bending
% 5. Optimizing sensor placement to approximate general state feedback...

%% Generate a bridge assembly and display
clear;clc;
A =         [1,1,1,1;
             1,0,0,1
            ];
A(:,:,2) =  [1,1,1,1
            1,0,0,1];
A = [1];
A(:,:,2) = 1; 
% A = [1,1]

lat = get_lattice(A);
cfg = config_array(lat);
stiffness_matrix = compute_matrix(cfg);
coords = beamCoords(lat,cfg);
size(coords)
plotBeams(cfg, coords);

%% code for compression/tension test
% use sensor readings to find global node displacements by assuming top
% face only moves in z rx ry and bottom face is fixed.
clear;clc;

A = [1;1;1];

lat = get_lattice(A);
cfg = config_array(lat);
stiffness_matrix = compute_matrix(cfg);
n = define_nodes(1);

sensor_beams = [n.f1.c, n.f1.bm
                n.f2.c, n.f2.bm
                n.f3.c, n.f3.bm];

sensor_readings = -[250e-6, 250e-6, 250e-6];
u = estimate_with_sensors(lat, stiffness_matrix, sensor_beams(:,1), sensor_beams(:,2), sensor_readings);
coords = beamCoords(lat,cfg);
plotBeams(cfg, coords, u, 50, sensor_beams);

%% code for bending test
clear;clc;
A = [1,1,1];
lat = get_lattice(A);
cfg = config_array(lat);
K = compute_matrix(cfg);

z_displ = -3e-2; %move by 3cm

%strategy: 
% middle vox gets pushed down from node tm, c, bm on top face -> imposed
% displacement in z
% first and last vox are supported from the same nodes on the bottom face->
% constrained to 0 displacement in xyz
bc_idx = []; bc_val = [];
n1 = define_nodes(1);
n2 = define_nodes(2);
n3 = define_nodes(3);
bot_nodes = [n1.f6.tm, n1.f6.c, n1.f6.bm, ...
             n3.f6.tm, n3.f6.c, n3.f6.bm];
% bot_nodes = [n1.f6.tr, n1.f6.br, ...
%              n3.f6.tl, n3.f6.bl]; %alternative test.

for node = bot_nodes
    bc_idx = [bc_idx, 6*node-5:6*node-3]; %x offset to z offset
    bc_val = [bc_val, [0,0,0]];
end

top_nodes = [n2.f5.tm, n2.f5.c, n2.f5.bm];
bc_idx = [bc_idx, top_nodes*6-3]; %z offset
bc_val = [bc_val, z_displ*[1,1,1]];

u = solve_with_dirichlet(K,bc_idx,bc_val);

coords = beamCoords(lat,cfg);
plotBeams(cfg, coords, u, 1);

%% template code for force BC

function nodeCoords = beamCoords(lat, beam_config)
    max_id = max(beam_config(:,1:2),[],'all');
    nodeCoords = (1:max_id)';
    out = arrayfun(@(x) get_xyz(x,lat), nodeCoords,'UniformOutput',false);
    nodeCoords = vertcat(out{:});
end

function plotBeams(beam_config, nodeCoords, U, scale, sensor_id)
% Plots undeformed (gray) and deformed (blue) beams.
% U : N×6 displacement matrix, translational DOFs in columns 1..3
    if nargin < 3, U = zeros(length(nodeCoords),6); end
    if nargin < 4, scale = 1; end
    if nargin < 5, sensor_id = [];end
    if width(U)~=6,U = reshape(U', 6,[])';end
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
