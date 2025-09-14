%define nodes


function n = define_nodes(A,r,c,s)
    %id_arr is a matrix containing voxel ids
    id_arr = init_arr(A);

    %face 1 (front face)
    base.b1 = 1; %bottom
    base.cl1 = 2; %corner left
    base.m1 = 3; %middle
    base.cr1 = 4; %corner right
    base.t1 = 5; %top

    %face 2 (right face)
    base.b2 = 6; %bottom
    base.cf2 = base.cr1; %corner front
    base.m2 = 7; %middle
    base.cb2 = 8; %corner back
    base.t2 = 9; %top

    %face 3 (back face)
    base.b3 = 10; %bottom
    base.cl3 = 11; %corner left (viewed from the front)
    base.m3 = 12; %middle
    base.cr3 = base.cb2; %corner right
    base.t3 = 13; %top

    %face 4 (left face)
    base.b4 = 14; %bottom
    base.cf4 = base.cl1; %corner front
    base.m4 = 15; %middle
    base.cb4 = base.cl3; %corner back
    base.t4 = 16; %top

    %face 5 (top face)
    base.f5 = base.t1; %front
    base.l5 = base.t4; %left
    base.m5 = 17; %middle
    base.r5 = base.t2; %right
    base.bk5 = base.t3; %back

    %face 6 (bottom face)
    base.f6 = base.b1; %front
    base.l6 = base.b4; %left
    base.m6 = 18; %middle
    base.r6 = base.b2; %right
    base.bk6 = base.b3; %back

    if ~A(r,c,s), error('No voxel at (%d,%d,%d).',r,c,s); end
    off = get_offset(base, id_arr(r,c,s));
    f = fieldnames(base);
    n = base;
    for k = 1:numel(f), n.(f{k}) = n.(f{k}) + off; end
end

function offset = get_offset(nodes, vox_id)
    all_vals = struct2array(nodes);         % dump all field values into one array
    nodes_per_vox = max(all_vals);  % assumes values are given sequentially
    offset = (vox_id - 1)*nodes_per_vox;% create offset
end