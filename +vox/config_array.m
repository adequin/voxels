function array_config = config_array(lat)
    % takes lattice struct created by get_lattice(A) and generates the config of the array
    % for each beam: [node1, node2, ex, ey, ez, length, young's mod, shear mod, width, height]
    if nargin == 0
        array_config = test();
    else
        nb_beams = 96; %number of beams in a voxel
        A = lat.A;
        conn_config = get_connector_edges(lat);
        n = define_nodes(1); % get base nodes
        base_config = get_voxel_config(n); %calculates parameters
        p = base_config(:, 3:end); %extract parameters
        array_config = [conn_config; base_config; zeros(nb_beams*(nnz(A)-1), 16)];
        offset = size(conn_config,1) + nb_beams; % skip connectors and base voxel

        for i = 2:(nnz(A)) % skip base voxel
            nodes = define_nodes(i);
            id_pairs = get_beam_ids(nodes);
            new_vox_config = [id_pairs, p]; % copy parameters to new id_pairs
            array_config(offset+(1:nb_beams),:) = new_vox_config; % write voxel into array_config
            offset = offset + nb_beams; %update index
        end
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
    % assert(isequal(size(conn_map.x),lat.dims)) 
    % assert(isequal(size(conn_map.y),lat.dims))
    % assert(isequal(size(conn_map.z),lat.dims))
end

function id_pairs = get_beam_ids(n)
    nb_beams = 96; %number of beams in a voxel
    id_pairs = zeros(nb_beams,2); % all id pairs in a voxel
    offset = 0;
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
        face_edges = [stub_edges; main_edges];
        nb_face_edges = size(face_edges,1);
        id_pairs(offset + (1:nb_face_edges),:) = face_edges;
        offset = offset + nb_face_edges;
    end
end

function voxel_config = get_voxel_config(n)
    %nodes: voxel node struct
    %DEFINE PARAMETERS
    b = 6e-3;% width in m
    h = 1.6e-3; % height in m
    % b_stub1 = 8.556e-3;
    % b_stub1 = 15e-3; 'actual' width
    b_stub = 23.8e-3; %equivalent b for PAT (probably affects J as well)
    % test values to iterate
    % I = h*b_stub1^3/12
    % I1 = h*b_stub^3/12
    % I2 = h*b_stub1^3/12 + b_stub1*h*(b_stub1/2)^2
    
    e = 25e9; %youngs modulus (25-30 GPA)in Pa
    g = 10e9; %shear modulus estimate in Pa
    face_beams = 16;
    offset = 0;

    % CONSTRUCT VOX CONFIG: [node1, node2, ex, ey, ez, length, young's mod, shear mod, width, height]
    voxel_config  = zeros(face_beams*6,16);
    map = coord_map();

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
            
            p1 = get_xyz(n1,[],map); p2 = get_xyz(n2,[],map);     
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

        p(1:size(stub_edges,1),13) = b_stub; % edits specific to stub parameters
        face_config = [edges,p];
        voxel_config(offset+(1:face_beams), :) = face_config; 
        offset = offset+face_beams; %update index offset
    end
    
end

function conn_config = find_conn_nodes(r,c,s,lat,cfg)
    %assumes you checked that there was a neighbor in x (cols), y (slices),
    %or z (rows)
    % r, c, s are the voxel coordinates

    %unpack cfg parameters
    shift = cfg.shift;
    normal = cfg.normal;
    up = cfg.up;
    f1 = cfg.face1;
    f2 = cfg.face2;
    P = cfg.pairs; % cell array N×2 of feature names
    
    r2 = r + shift(1); c2 = c + shift(2); s2 = s + shift(3);

    if any(~[r,c,s,r2,c2,s2]) %bound check
        error('Expected neighbor but got zero id at (%d,%d,%d) or (%d,%d,%d).', r,c,s,r2,c2,s2);
    end

    n1 = define_nodes(lat, r,c,s);
    n2 = define_nodes(lat,r2,c2,s2);

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
    cfgx = struct( ...
        'shift', [0 1 0], ...
        'normal', [1 0 0], ...
        'up', [0 0 1], ...
        'face1', "f2", ...
        'face2', "f4", ...
        'pairs', { {'bm','bm'; 'rm','lm'; 'tm','tm'; 'lm','rm'} } );

    cfgy = struct( ...
        'shift', [0 0 1], ...
        'normal', [0 1 0], ...
        'up', [0 0 1], ...
        'face1', "f3", ...
        'face2', "f1", ...
        'pairs', { {'bm','bm'; 'rm','lm'; 'tm','tm'; 'lm','rm'} } );

    cfgz = struct( ...
        'shift', [1 0 0], ...
        'normal', [0 0 1], ...
        'up', [0 1 0], ...
        'face1', "f6", ...
        'face2', "f5", ...
        'pairs', { {'bm','tm'; 'rm','rm'; 'tm','bm'; 'lm','lm'} } );

    conn_map = get_conn_map(lat);
    pairs_per = 4;
    nX = nnz(conn_map.x); nY = nnz(conn_map.y); nZ = nnz(conn_map.z);
    total = (nX+nY+nZ)*pairs_per;
    conn_config = zeros(total, 2+9);
    idx = 0;

    id_arr = lat.id;
    dims = size(id_arr);

    [rx,cx,sx] = ind2sub(dims, find(conn_map.x));
    for t = 1:numel(rx) %same number of elements for rx, cx, sx
        idx = idx(end) + (1:pairs_per);  %update indicies
        conn_config(idx,:) = find_conn_nodes(rx(t),cx(t),sx(t),lat,cfgx);
    end

    [ry,cy,sy] = ind2sub(dims, find(conn_map.y));
    for t = 1:numel(ry)
        idx = idx(end) + (1:pairs_per); %update indicies
        conn_config(idx,:) = find_conn_nodes(ry(t),cy(t),sy(t),lat,cfgy);
    end

    [rz,cz,sz] = ind2sub(dims, find(conn_map.z));
    for t = 1:numel(rz)
        idx = idx(end) + (1:pairs_per);  %update indicies
        conn_config(idx,:) = find_conn_nodes(rz(t),cz(t),sz(t),lat,cfgz);
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

function c = test
% put test code
% something like this.
% clear;clc;
% A =         [1,0;
%             1,1 ];
% A(:,:,2) =  [0,0;
%             1,1];
% 
A = ones(2,2,2);
lat = get_lattice(A);
n = define_nodes(lat.id(1,1,1));
% c = get_voxel_config(n)
% c = get_beam_ids(n);
% size(c)
% conn = get_conn_map(lat);
c = get_connector_edges(lat);
% c = config_array(lat);
end