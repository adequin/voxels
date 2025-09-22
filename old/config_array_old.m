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
% A =         [1,0;
%             1,0 ];
% A(:,:,2) =  [0,0;
%             1,1];

% vox_arr = init_arr(A);

% % define_nodes(A,1,1,1)
% conn = get_conn_map(A);
% conn.x, conn.z
% get_voxel_config(define_nodes(A,1,1,1), 2);
% arr = get_array_config(A);
% size(arr)


function array_config = get_array_config(A)
    id_arr = init_arr(A);
    nodes = define_nodes(A,1,1,1);
    conn_map = get_conn_map(A);
    array_config = get_connector_edges(conn_map, id_arr, nodes);
    for i = 1:nnz(A)
        array_config = [array_config; get_voxel_config(nodes, i)];
    end

end

function conn_map = get_conn_map(A)
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
end

function offset = get_offset(nodes, vox_id)
    all_vals = struct2array(nodes);         % dump all field values into one array
    nodes_per_vox = max(all_vals);  % assumes values are given sequentially
    offset = (vox_id - 1)*nodes_per_vox;% create offset
end

function voxel_config = get_voxel_config(nodes, vox_id)
    %id: voxel number
    %vox_id: sequential number assigned to vox in array (indexed by grid coords)
    id_off = get_offset(nodes, vox_id);% create offset

    %OFFSET NODES
    %face 1
    b1 = id_off+nodes.b1; %bottom
    cl1 = id_off+nodes.cl1; %corner left
    m1 = id_off+nodes.m1; %middle
    cr1 = id_off+nodes.cr1; %corner right
    t1 = id_off+nodes.t1; %top

    %face 2 (right face)
    b2 = id_off+nodes.b2; %bottom
    cf2 = cr1; %corner front
    m2 = id_off+nodes.m2; %middle
    cb2 = id_off+nodes.cb2; %corner back
    t2 = id_off+nodes.t2; %top

    %face 3 (back face)
    b3 = id_off+nodes.b3; %bottom
    cl3 = id_off+nodes.cl3; %corner left (viewed from the front)
    m3 = id_off+nodes.m3; %middle
    cr3 = cb2; %corner right
    t3 = id_off+nodes.t3; %top

    %face 4 (left face)
    b4 = id_off+nodes.b4; %bottom
    cf4 = cl1; %corner front
    m4 = id_off+nodes.m4; %middle
    cb4 = cl3; %corner back
    t4 = id_off+nodes.t4; %top

    %face 5 (top face)
    f5 = t1; %front
    l5 = t4; %left
    m5 = id_off+nodes.m5; %middle
    r5 = t2; %right
    bk5 = t3; %back

    %face 6 (bottom face)
    f6 = b1; %front
    l6 = b4; %left
    m6 = id_off+nodes.m6; %middle
    r6 = b2; %right
    bk6 = b3; %back

    %DEFINE PARAMETERS
    b = 1.6e-3;% width in m
    h = 6e-3; % height in m
    % b = h
    e = 25e9; %youngs modulus (25-30 GPA)in Pa
    g = 10e9; %shear modulus estimate in Pa
    l = 47.686e-3; %m
    l_45 = l*sqrt(2);
    offset = pi/2;

    %CONSTRUCT VOX CONFIG: [node1, node2, x_rot, y_rot, z_rot, length, young's mod, width, height]
    voxel_config  = [
    %face1
                b1,cl1,offset,0,3*pi/4, l_45, e, g, b, h;
                b1,m1,offset,0,pi/2, l, e, g, b, h;
                b1,cr1,offset,0,pi/4,l_45, e, g, b, h;
                cl1,t1,offset,0,pi/4,l_45, e, g, b, h;
                m1,t1,offset,0,pi/2, l, e, g, b, h;
                cr1,t1,offset,0,3*pi/4, l_45, e, g, b, h;
                cl1,m1,offset,0,0, l, e, g, b, h;
                m1,cr1,offset,0,0, l, e, g, b, h;
    %face 2
                cf2,b2,offset,pi/2,-pi/4,l_45, e, g, b, h;
                cf2,m2,offset,pi/2,0, l, e, g, b, h;
                cf2,t2,offset,pi/2,pi/4,l_45, e, g, b, h;

                b2,cb2,offset,pi/2,pi/4,l_45, e, g, b, h;
                m2,cb2,offset,pi/2,0, l, e, g, b, h;
                t2,cb2,offset,pi/2,-pi/4,l_45, e, g, b, h;

                b2,m2,offset,pi/2,pi/2, l, e, g, b, h;
                m2,t2,offset,pi/2,pi/2, l, e, g, b, h;

   %face 3
                cl3,b3,offset,0,-pi/4,l_45, e, g, b, h;
                cl3,m3,offset,0,0,l, e, g, b, h;
                cl3,t3,offset,0,pi/4,l_45, e, g, b, h;

                b3,cr3,offset,0,pi/4,l_45, e, g, b, h;
                m3,cr3,offset,0,0,l, e, g, b, h;
                t3,cr3,offset,0,-pi/4,l_45, e, g, b, h;

                b3,m3,offset,0,pi/2,l, e, g, b, h;
                m3,t3,offset,0,pi/2,l, e, g, b, h;

   %face 4
                cf4,b4,offset,pi/2,-pi/4,l_45, e, g, b, h;
                cf4,m4,offset,pi/2,0,l, e, g, b, h;
                cf4,t4,offset,pi/2,pi/4,l_45, e, g, b, h;

                b4,cb4,offset,pi/2,pi/4,l_45, e, g, b, h;
                m4,cb4,offset,pi/2,0,l, e, g, b, h;
                t4,cb4,offset,pi/2,-pi/4,l_45, e, g, b, h;

                b4,m4,offset,pi/2,pi/2,l, e, g, b, h;
                m4,t4,offset,pi/2,pi/2,l, e, g, b, h;

   %face 5 (top)

                f5,l5,-pi/2+offset,0,3*pi/4,l_45, e, g, b, h;
                f5,m5,-pi/2+offset,0,pi/2,l, e, g, b, h;
                f5,r5,-pi/2+offset,0,pi/4,l_45, e, g, b, h;

                l5,bk5,-pi/2+offset,0,pi/4,l_45, e, g, b, h;
                m5,bk5,-pi/2+offset,0,pi/2,l, e, g, b, h;
                r5,bk5,-pi/2+offset,0,3*pi/4,l_45, e, g, b, h;

                l5,m5,-pi/2+offset,0,0,l, e, g, b, h;
                m5,r5,-pi/2+offset,0,0,l, e, g, b, h;

%face 6 (bottom)
                f6,l6,-pi/2+offset,0,3*pi/4,l_45, e, g, b, h;
                f6,m6,-pi/2+offset,0,pi/2,l, e, g, b, h;
                f6,r6,-pi/2+offset,0,pi/4,l_45, e, g, b, h;

                l6,bk6,-pi/2+offset,0,pi/4,l_45, e, g, b, h;
                m6,bk6,-pi/2+offset,0,pi/2,l, e, g, b, h;
                r6,bk6,-pi/2+offset,0,3*pi/4,l_45, e, g, b, h;

                l6,m6,-pi/2+offset,0,0,l, e, g, b, h;
                m6,r6,-pi/2+offset,0,0,l, e, g, b, h;
    ]; 

end

function conn_config = find_conn_nodes_x(nodes, r,c,s, id_arr)
    %assumes you checked that there was a neighbor in x (cols)
    % r, c, s
    vox_id1 = id_arr(r,c,s);
    vox_id2 = id_arr(r,c+1,s);
    off1 = get_offset(nodes, vox_id1);
    off2 = get_offset(nodes, vox_id2);
    new_id1 = @(node) node + off1;
    new_id2 = @(node) node + off2;
    
    id_pairs = [new_id1(nodes.b2),new_id2(nodes.b4);
                new_id1(nodes.cf2),new_id2(nodes.cf4);
                new_id1(nodes.cb2),new_id2(nodes.cb4);
                new_id1(nodes.t2),new_id2(nodes.t4)];

    rotations = [pi/2,0,0;
                 0,0,0;
                 0,0,0;
                 pi/2,0,0];

    conn_config = [id_pairs, rotations];
end

function conn_config = find_conn_nodes_y(nodes, r,c,s, id_arr)
    %assumes you checked that there was a neighbor in y (slices)
    vox_id1 = id_arr(r,c,s);
    vox_id2 = id_arr(r,c,s+1);
    off1 = get_offset(nodes, vox_id1);
    off2 = get_offset(nodes, vox_id2);
    new_id1 = @(node) node + off1;
    new_id2 = @(node) node + off2;
    
    id_pairs = [new_id1(nodes.b3),new_id2(nodes.b1);
                new_id1(nodes.cl3),new_id2(nodes.cl1);
                new_id1(nodes.cr3),new_id2(nodes.cr1);
                new_id1(nodes.t3),new_id2(nodes.t1)];

    rotations = [pi/2,pi/2,0;
                 0,0,pi/2;
                 0,0,pi/2;
                 pi/2,pi/2,0];
    
    conn_config = [id_pairs, rotations];
end

function conn_config = find_conn_nodes_z(nodes, r,c,s, id_arr)
    %assumes you checked that there was a neighbor in x (cols)
    vox_id1 = id_arr(r,c,s);
    vox_id2 = id_arr(r+1,c,s);
    off1 = get_offset(nodes, vox_id1);
    off2 = get_offset(nodes, vox_id2);
    new_id1 = @(node) node + off1;
    new_id2 = @(node) node + off2;
    
    id_pairs = [new_id1(nodes.b1),new_id2(nodes.t1);
                new_id1(nodes.b2),new_id2(nodes.t2);
                new_id1(nodes.b3),new_id2(nodes.t3);
                new_id1(nodes.b4),new_id2(nodes.t4)];

    rotations = [0,-pi/2,0;
                 -pi/2,0,pi/2;
                 0,-pi/2,0;
                 -pi/2,0,pi/2];
    conn_config = [id_pairs, rotations];
end

function conn_config = get_connector_edges(conn_map, id_arr, nodes)
    
    conn_config = [];
    if any(conn_map.x,'all')
        for id = id_arr(logical(conn_map.x))'
            idx = find(id_arr==id, 1);            % linear index
            [rx,cx,sx] = ind2sub(size(id_arr), idx);
            conn_config = [conn_config; find_conn_nodes_x(nodes,rx,cx,sx,id_arr)];
        end
    end

    if any(conn_map.y,'all')
        for id = id_arr(logical(conn_map.y))'
            idx = find(id_arr==id, 1);            % linear index
            [ry,cy,sy] = ind2sub(size(id_arr), idx);
            conn_config = [conn_config; find_conn_nodes_y(nodes,ry,cy,sy,id_arr)];
        end
    end
    if any(conn_map.z,'all')
        for id = id_arr(logical(conn_map.z))'
            idx = find(id_arr==id, 1);            % linear index
            [rz,cz,sz] = ind2sub(size(id_arr), idx);
            conn_config = [conn_config; find_conn_nodes_z(nodes,rz,cz,sz,id_arr)];
        end
    end
    
    %connector parameters:
    b = 1.6e-3;% width in m
    h = 6e-3; % height in m

    e = 25e9; % youngs modulus (25-30 GPA)in Pa
    g = 10e9; %shear modulus estimate in Pa
    l = 3e-3; % length of connector in m
    param = [l, e, g, b, h];

    %append parameters
    param_arr = repmat(param, size(conn_config, 1), 1, 1);
    conn_config = [conn_config, param_arr];
end

