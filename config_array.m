%% what the hell am I doing
% define voxel assembly here:
% x   o   x
% x   x   x

%what are all the things to keep track of?
% attributes
% - node names
% - connectors locations
% - voxel edges
% 
% functions
% - create voxel array ids
% - find connector map using array
% - create a voxel using id
% - create connectors using map



%%test code
% something like this
% clear;clc;
A =         [1,0;
            1,1 ];
A(:,:,2) =  [0,0;
            1,1];

% vox_arr = get_lattice(A);

% % define_nodes(A,1,1,1)
% conn = get_conn_map(A);
% conn.x, conn.z
% get_voxel_config(define_nodes(A,1,1,1), 2);
% arr = get_array_config(A);
% size(arr)

%new test code
lat = get_lattice(A);
n = define_nodes(lat.id(1,1,1));
c = get_voxel_config(n);
size(c)
conn = get_conn_map(lat);
conn_config = get_connector_edges(lat)
array_config = config_arr(A)

function array_config = config_arr(A) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   <----------------------
    lat = get_lattice(A);
    array_config = get_connector_edges(lat);
    for i = 1:nnz(A)
        nodes = define_nodes(i);
        array_config = [array_config; get_voxel_config(nodes)]; %we can just use i because spacial relations are determined by connectors
    end

end

function conn_map = get_conn_map(lat)
    A = lat.A;
    % z connectors (rows)
    [Rows,Cols,Slices] = size(A);
    mask = false(size(A));

    conn_map.x = mask;
    conn_map.y = mask;
    conn_map.z = mask;

    if Rows > 1
        A1 = A(1:end-1,:,:); B1 = A(2:end,:,:); % get two copies shifted by one
        M1 = (A1>0) & (B1>0);
        conn_map.z(1:Rows-1,:,:) = M1;
    end

    % x connectors (cols)
    if Cols > 1
        A2 = A(:,1:end-1,:); B2 = A(:,2:end,:);
        M2 = (A2>0) & (B2>0);
        conn_map.x(:,1:Cols-1,:) = M2;
    end

    % y connectors (layers)
    if Slices > 1
        A3 = A(:,:,1:end-1); B3 = A(:,:,2:end);
        M3 = (A3>0) & (B3>0);
        conn_map.y(:,:,1:Slices-1) = M3;
    end
    assert(isequal(size(conn_map.x),lat.dims)) 
    assert(isequal(size(conn_map.y),lat.dims))
    assert(isequal(size(conn_map.z),lat.dims))
end

function voxel_config = get_voxel_config(n)
    %nodes: voxel node struct
    
    %DEFINE PARAMETERS
    b = 6e-3;% width in m
    h = 1.6e-3; % height in m
    b_stub = 8.556e-3;
    e = 25e9; %youngs modulus (25-30 GPA)in Pa
    g = 10e9; %shear modulus estimate in Pa


    % CONSTRUCT VOX CONFIG: [node1, node2, length, young's mod, shear mod, width, height, theta, face_normal, face_right]
    voxel_config  = [];
    for face = n.face_names
        f = n.(face);
        stub_edges = [
                     %small edges
                     f.bl, f.bm;
                     f.bm, f.br;
                     f.rb, f.rm;
                     f.rm, f.rt;
                     f.tr, f.tm;
                     f.tm, f.tl;
                     f.lt, f.lm;
                     f.lm, f.lb];
        main_edges = [
                     %cross edges
                     f.lm, f.c;
                     f.c,  f.rm;
                     f.tm, f.c;
                     f.c,  f.bm;
                     
                     %diag edges
                     f.br, f.rb;
                     f.rt, f.tr;
                     f.tl, f.lt;
                     f.lb, f.bl;];
        edges = [stub_edges; main_edges];

        p = zeros(size(edges,1),14);
        for i = 1:size(edges,1)
            n1 = edges(i,1); n2 = edges(i,2);
            p1 = get_xyz(n1); p2 = get_xyz(n2);     
            vec = p2-p1;
            normal = f.normal;
            up = f.up;
            [ex, ey, ez] = local_basis(unit(vec), normal, up);
            p(i,1:9) = [ex, ey, ez];
            p(i,10) = norm(vec);
            p(i,11) = e;
            p(i,12) = g;
            p(i,13) = b;
            p(i,14) = h;
            
        end
        for j = 1:size(stub_edges,1) %edits specific to stub parameters
            p(j,13) = b_stub;
        end
        face_config = [edges,p];
        voxel_config = [voxel_config ; face_config];
    end
    
end

function conn_config = find_conn_nodes(r,c,s,lat,cfg)
    %assumes you checked that there was a neighbor in x (cols), y (slices),
    %or z (rows)
    % r, c, s are the voxel coordinates
    % CONSTRUCT CONN CONFIG: [node1, node2, length, young's mod, shear mod, width, height, theta, face_normal]

    %unpack cfg parameters
    shift = cfg.shift;
    normal = cfg.normal;
    up = cfg.up;
    f1 = cfg.face1;
    f2 = cfg.face2;
    P = cfg.pairs; % cell array N×2 of feature names
    
    r2 = r + shift(1); c2 = c + shift(2); s2 = s + shift(3);

    id_arr = lat.id;
    vox_id1 = id_arr(r,c,s);
    vox_id2 = id_arr(r2,c2,s2);

    if any(~[vox_id1,vox_id2]) %bound check
        error('Expected neighbor but got zero id at (%d,%d,%d) or (%d,%d,%d).', r,c,s,r2,c2,s2);
    end

    n1 = define_nodes(vox_id1);
    n2 = define_nodes(vox_id2);

    rot = zeros(size(P,1),9);  
    id_pairs = zeros(size(P,1), 2);
    flip_flag = 1;% cfg starts with bottom connector then rotates ccw, default is upright so we want to flip first
    for k = 1:size(P,1)
        id_pairs(k,1) = n1.(f1).(P{k,1});   % e.g., f1.bm
        id_pairs(k,2) = n2.(f2).(P{k,2});   % e.g., f2.lm

        %get basis vectors
        axis = normal; %can't use local voxel coords for this, assume they stick out the face
        [ex, ey, ez] = local_basis(axis, normal, up, flip_flag);
        rot(k,:) = [ex,ey,ez]; %compute basis vectors
        flip_flag = ~flip_flag; %rotates every other connector, assumes they are listed consecutively
    end
    conn_config = [id_pairs, rot];
end

function conn_config = get_connector_edges(lat)
    cfg.x = struct( ...
        'shift', [0 1 0], ...
        'normal', [1 0 0], ...
        'up', [0 0 1], ...
        'face1', "f2", ...
        'face2', "f4", ...
        'pairs', { {'bm','bm'; 'rm','lm'; 'tm','tm'; 'lm','rm'} } );

    cfg.y = struct( ...
        'shift', [0 0 1], ...
        'normal', [0 1 0], ...
        'up', [0 0 1], ...
        'face1', "f3", ...
        'face2', "f1", ...
        'pairs', { {'bm','bm'; 'rm','lm'; 'tm','tm'; 'lm','rm'} } );

    cfg.z = struct( ...
        'shift', [1 0 0], ...
        'normal', [0 0 1], ...
        'up', [0 1 0], ...
        'face1', "f5", ...
        'face2', "f6", ...
        'pairs', { {'bm','tm'; 'rm','lm'; 'tm','bm'; 'lm','rm'} } );

    conn_config = [];
    conn_map = get_conn_map(lat);
    id_arr = lat.id;
    if any(conn_map.x,'all')
        for id = id_arr(logical(conn_map.x))'
            idx = find(id_arr==id, 1);            % linear index
            [rx,cx,sx] = ind2sub(size(id_arr), idx);
            conn_config = [conn_config; find_conn_nodes(rx,cx,sx,lat,cfg.x)];
        end
    end
    if any(conn_map.y,'all')
        for id = id_arr(logical(conn_map.y))'
            idx = find(id_arr==id, 1);            % linear index
            [ry,cy,sy] = ind2sub(size(id_arr), idx);
            conn_config = [conn_config; find_conn_nodes(ry,cy,sy,lat,cfg.y)];
        end
    end
    if any(conn_map.z,'all')
        for id = id_arr(logical(conn_map.z))'
            idx = find(id_arr==id, 1);            % linear index
            [rz,cz,sz] = ind2sub(size(id_arr), idx);
            conn_config = [conn_config; find_conn_nodes(rz,cz,sz,lat,cfg.z)];
        end
    end
    
    %connector parameters:
    b = lat.conn_b;% width in m
    h = lat.conn_h; % height in m

    e = lat.conn_e; % youngs modulus (25-30 GPA)in Pa
    g = lat.conn_g; %shear modulus estimate in Pa
    l = lat.conn_l; % length of connector in m
    param = [l, e, g, b, h];

    %append parameters
    param_arr = repmat(param, size(conn_config, 1), 1, 1);
    conn_config = [conn_config, param_arr];
end

function coords = get_xyz(id)
%returns relative coordinates from voxel center. Assumes you never need
%global coords
    id_loc = mod(id-1,42)+1; %shift back to base node coords (need to account for mod(42,42)->0)
    map = coord_map();

    if id_loc<=36
        edge  = ceil(id_loc/3);
        slot  = id_loc - 3*(edge-1);
        coords = reshape(map.edge_points(edge,slot,:),1,[]);

    elseif id_loc <=42
        coords = map.face_points(id_loc-36,:);
    else
        error('error mapping id = %d',id_loc);
    end

    assert(isequal(size(coords), [1 3]))
end

function [ex, ey, ez] = local_basis(axis,normal,up,flip_flag)
    %returns three basis row vectors to shift general beam stiffness matrix into
    %beam defined by p1, p2 on face defined by normal and up vectors.
    %handles connectors parallel to face normal as well.
    if nargin<5
        flip_flag = 0;
    end

    tol = 1e-9;
    axis = unit(axis);

    if iscolumn(normal), normal = normal.'; end
    if iscolumn(up),     up     = up.';     end
    if iscolumn(axis),   axis    = axis.';    end

    ex = axis/norm(axis);
    test = abs(dot(normal, ex));
    if test < tol
        ey = -normal; %normal case
        ez = unit(cross(ex, ey));
    elseif abs(test - 1) < tol
        ez = up; %connector case
        ey = unit(cross(ez, ex)); %should be equal to right face vector
    else
        error("off axis beam. dot product with face normal = %d, " + ...
            "ex = %d, %d, %d, " + ...
            "face normal = %d, %d, %d, ", test, ex(1), ex(2), ex(3), ...
            normal(1), normal(2), normal(3));
    end

    if flip_flag %rotate by pi/2 axially
            [ey, ez] = deal(ez, -ey); 
    end
    % checks
    assert(~all([dot(ex,ey), dot(ey,ez), dot(ez,ex)]))
    assert(norm(ex)==norm(ey)==norm(ez))

end

function v_normalized = unit(v)
    v_normalized = v/norm(v);
end