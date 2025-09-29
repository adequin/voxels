function lat = get_lattice(A,conn_l, conn_b, conn_h, conn_e, conn_g, voxel_size, chamfer_size)
    %creates struct with voxel ids
    %defines connector dimensions in stiffness matrix
    %A is a logical matrix indicating voxel positions

    %default connector values
    l = 3e-3; % length in m
    b = 20e-3;% width in m
    h = 1.6e-3; % height in m
    e = 30e9; % youngs modulus (25-30 GPA)in Pa
    g = 10e9; %shear modulus estimate in Pa

    %default voxel 
    c = 59.243e-3;%m cube edge offset
    v = 140.088e-3;%m "cube" edge length

    if nargin<8
        chamfer_size = c;
    end
    if nargin<7
        voxel_size = v;
    end
    if  nargin<6
        conn_g = g;
    end
    if nargin<5
        conn_e = e;
    end
    if nargin<4
        conn_h = h;
    end
    if  nargin<3
        conn_b = b;
    end
    if nargin<2
        conn_l = l;
    end

    A = logical(A);
    lat.A = A;
    lat.dims = size(A);
    lat.id = zeros(lat.dims,'double'); % 
    lat.id(A) = double(1:nnz(A)); % replace true entries with sequential values. nnz-> number of non zero values
    
    lat.conn_l = conn_l; %connector length
    lat.conn_b = conn_b;
    lat.conn_h = conn_h;
    lat.conn_e = conn_e;
    lat.conn_g = conn_g;
    lat.voxel_size = voxel_size;
    lat.chamfer = chamfer_size; 
end