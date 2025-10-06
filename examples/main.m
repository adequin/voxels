%remaining issues
% 1. Bending inertias of stub beams. currently they are implemented as two
%rectangular beams that have a cross section of +, but probably should be
%offset by setting individual inertias offset using parallel axis.
% 2. Length parameters are still a bit arbitrary and only approximate all
% the fillets and non-uniform geometry
% 3. Better UI for constraining dofs instead of doing 6*node-5 or smth
% 4. Sensor placement and interpretation to for bending
% 5. Optimizing sensor placement to approximate general state feedback...
%% 
clear;clc;
runtime = [];
              
for dim = 1:12
    tic % start timer
    z_displ = -3e-2;
    
    A   = ones(dim,dim,dim);
    lat = get_lattice(A);
    
    % Collect nodes in cells (cheap), concatenate once, then deduplicate
    bot_cells = cell(dim,dim);
    top_cells = cell(dim,dim);
    
    for i = 1:dim
        for j = 1:dim
            n  = define_nodes(lat, dim, i, j);
            n2 = define_nodes(lat, 1,   i, j);
            bot_cells{i,j} = get_face_nodes(n,  'f6');
            top_cells{i,j} = get_face_nodes(n2, 'f5');
        end
    end
    
    bot_nodes = unique([bot_cells{:}], 'stable');
    top_nodes = unique([top_cells{:}], 'stable');
    
    % Build DOF indices vectorized (x,y,z DOFs for bottom, z-DOF for top)
    bot = bot_nodes(:);
    top = top_nodes(:);
    
    bot_xyz = [bot*6-5, bot*6-4, bot*6-3].';   % 3-by-N
    bot_xyz = bot_xyz(:).';                    % row vector
    
    top_z   = (top*6-3).';                     % row vector
    
    bc_idx = [bot_xyz, top_z];
    
    % Matching values: zeros for bottom xyz, z_displ for top z
    bc_val = [zeros(1, numel(bot)*3), repmat(z_displ, 1, numel(top))];
    
    % Build config once
    cfg    = config_array(lat);
    
    % If K is huge, ensure compute_matrix returns sparse
    K      = compute_matrix(cfg);
    
    % coords only if you need them
    % coords = beamCoords(lat,cfg);
    u      = solve_with_dirichlet(K, bc_idx, bc_val);
    % plotBeams(cfg, coords, u);
    elapsedTime = toc; % stop timer and return elapsed time in seconds
    runtime = [runtime, elapsedTime];
end
plot(runtime);
improvePlot;
%%
%% Runtime scaling benchmark
clear; clc;

dims = 1:4;                  % test sizes
nrep = 6;                     % repeat per size; take median to reduce noise

runtime_total   = zeros(size(dims));
runtime_config  = zeros(size(dims));
runtime_K       = zeros(size(dims));
runtime_solve   = zeros(size(dims));
ndof            = zeros(size(dims));
nelem           = zeros(size(dims));

for ii = 1:numel(dims)
    dim = dims(ii);

    t_total = zeros(nrep,1);
    t_cfg   = zeros(nrep,1);
    t_K     = zeros(nrep,1);
    t_sol   = zeros(nrep,1);

    for r = 1:nrep
        z_displ = -3e-2;
        A   = ones(dim,dim,dim);
        lat = get_lattice(A);

        % --- build BCs (same as your script) ---
        bot_cells = cell(dim,dim);
        top_cells = cell(dim,dim);
        for i = 1:dim
            for j = 1:dim
                n  = define_nodes(lat, dim, i, j);
                n2 = define_nodes(lat, 1,   i, j);
                bot_cells{i,j} = get_face_nodes(n,  'f6');
                top_cells{i,j} = get_face_nodes(n2, 'f5');
            end
        end
        bot_nodes = unique([bot_cells{:}], 'stable');
        top_nodes = unique([top_cells{:}], 'stable');

        bot = bot_nodes(:);
        top = top_nodes(:);
        bot_xyz = [bot*6-5, bot*6-4, bot*6-3].';
        bot_xyz = bot_xyz(:).';
        top_z   = (top*6-3).';
        bc_idx  = [bot_xyz, top_z];
        bc_val  = [zeros(1, numel(bot)*3), repmat(z_displ, 1, numel(top))];

        t0 = tic;

        % --- config_array ---
        t1 = tic;
        cfg = config_array(lat);
        t_cfg(r) = toc(t1);

        % --- K assembly ---
        t2 = tic;
        K  = compute_matrix(cfg);    % should return sparse
        t_K(r) = toc(t2);

        % --- solve ---
        t3 = tic;
        u  = solve_with_dirichlet(K, bc_idx, bc_val);
        t_sol(r) = toc(t3);

        t_total(r) = toc(t0);

        if r==1
            nelem(ii) = size(cfg,1);
            ndof(ii)  = size(K,1);
        end
    end

    % robust aggregate
    runtime_total(ii)  = median(t_total);
    runtime_config(ii) = median(t_cfg);
    runtime_K(ii)      = median(t_K);
    runtime_solve(ii)  = median(t_sol);
end

% --- Plots ---
figure; 
subplot(1,2,1);
plot(dims, runtime_total, '-o'); hold on;
plot(dims, runtime_config, '-s');
plot(dims, runtime_K, '-^');
plot(dims, runtime_solve, '-d');
xlabel('dim'); ylabel('time (s)'); grid on;
legend('total','config','assemble K','solve','Location','northwest');
title('Runtime vs dim');

subplot(1,2,2);
loglog(ndof, runtime_total, '-o'); hold on;
loglog(ndof, runtime_config, '-s');
loglog(ndof, runtime_K, '-^');
loglog(ndof, runtime_solve, '-d');
xlabel('DOFs'); ylabel('time (s)'); grid on;
legend('total','config','assemble K','solve','Location','northwest');
title('Runtime vs DOFs (log-log)');

%% Generate a bridge assembly and display
clear;clc;

% A = [1];
% A(:,:,2) = 1; 
% A = [1,1]
dim = 12;
z_displ = -3e-2;
A = ones(dim,dim,dim);
lat = get_lattice(A);
bot_nodes = [];
top_nodes = [];
for i = 1:dim
    for j = 1:dim
        n = define_nodes(lat,dim,i,j);
        n2= define_nodes(lat,1,i,j);
        bot_nodes = [bot_nodes, get_face_nodes(n,'f6')]; 
        top_nodes = [top_nodes, get_face_nodes(n2,'f5')];
    end
end

bc_idx = [];
bc_val = [];
for node = bot_nodes
    bc_idx = [bc_idx, 6*node-5:6*node-3]; %x offset to z offset
    bc_val = [bc_val, [0,0,0]];
end

bc_idx = [bc_idx, top_nodes*6-3]; %z offset
bc_val = [bc_val, z_displ*ones(1,length(top_nodes))];

%% 
clear;clc;
A =         [1,1,1,1;
             1,0,0,1
            ];
A(:,:,2) =  [1,1,1,1
            1,0,0,1];
lat = get_lattice(A);
cfg = config_array(lat);
K = compute_matrix(cfg);
coords = beamCoords(lat,cfg);
% u = solve_with_dirichlet(K,bc_idx,bc_val);
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
    map = coord_map(lat.voxel_size, lat.chamfer);
    out = arrayfun(@(x) get_xyz(x,lat,map), nodeCoords,'UniformOutput',false);
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
    if nargin<3, exclude_face = false; end
    names = n.face_nodes;
    if exclude_face, names = names(1:end-1); end

    nodes = zeros(1, numel(names));
    for j = 1:numel(names)
        nodes(j) = n.(face).(names{j});
    end
end
