function K_global = compute_matrix(beam_data_3D)
    % example K_global row: [node1, node2, local axis, local y axis, local z axis, L,E,G,B,H]

    % --- global DOF size ---
    num_nodes = max(max(beam_data_3D(:,1:2)));
    nDOF = 6 * num_nodes;

    % --- preallocate triplets: 144 entries per beam (12x12) ---
    nElem = size(beam_data_3D,1);
    I_all = zeros(144*nElem,1);
    J_all = zeros(144*nElem,1);
    V_all = zeros(144*nElem,1);
    ptr = 0;

    for i = 1:nElem
        % UNPACK VOX CONFIG: [node1, node2, ex, ey, ez, length, young's mod, shear mod, width, height]
        node1  = beam_data_3D(i,1);
        node2  = beam_data_3D(i,2);
        ex = beam_data_3D(i,3:5); %row vectors
        ey = beam_data_3D(i,6:8);
        ez = beam_data_3D(i,9:11);
        L      = beam_data_3D(i,12);
        E      = beam_data_3D(i,13);
        G      = beam_data_3D(i,14);
        B      = beam_data_3D(i,15);
        H      = beam_data_3D(i,16);
        A = B*H;
        Iy = H*B^3/12;
        Iz = B*H^3/12; % area moment of inertia about z
        % Ix = Iy+Iz;
        J  = B*H^3*(1/3 - 0.21*(H/B)*(1 - (H^4)/(12*B^4)));  % Saint-Venant approx

        % local stiffness
        K_local = beam_stiffness_matrix_3D(E, G, A, J, Iy, Iz, L);

        % rotation (12x12 from your 3x3)
        R = [ex', ey', ez'];  R = blkdiag(R, R, R, R); % should be 12x12
        K_beam = R * K_local * R.';                    % element in GLOBAL coords

        % triplets for this element
        dof = [6*node1-5:6*node1, 6*node2-5:6*node2];    % 1x12
        [JJ,II] = meshgrid(dof, dof);                    % 12x12 each

        % write into preallocated slots
        I_all(ptr+1:ptr+144) = II(:);
        J_all(ptr+1:ptr+144) = JJ(:);
        V_all(ptr+1:ptr+144) = K_beam(:);
        ptr = ptr + 144;
    end

    % --- build global sparse once ---
    K_global = sparse(I_all, J_all, V_all, nDOF, nDOF);
    K_global = (K_global + K_global.')/2;  % enforce symmetry (numerical safety)
end

function K_local = beam_stiffness_matrix_3D(E, G, A, Ix, Iy, Iz, L)
    % Returns the 12x12 local stiffness matrix for a 3D beam element
    K_node1 = [
        E*A/L,         0,            0,           0,        0,          0;
        0,             12*E*Iz/L^3,  0,           0,        0,          6*E*Iz/L^2
        0,             0,            12*E*Iy/L^3, 0,        -6*E*Iy/L^2,0;
        0,             0,            0,           G*Ix/L,   0,          0;
        0,             0,            -6*E*Iy/L^2, 0,        4*E*Iy/L,   0;
        0,             6*E*Iz/L^2,   0,           0,        0,          4*E*Iz/L
    ];
    K_node2 = [
        E*A/L,         0,            0,           0,        0,          0;
        0,             12*E*Iz/L^3,  0,           0,        0,          -6*E*Iz/L^2
        0,             0,            12*E*Iy/L^3, 0,        6*E*Iy/L^2,0;
        0,             0,            0,           G*Ix/L,   0,          0;
        0,             0,            6*E*Iy/L^2,  0,        4*E*Iy/L,   0;
        0,             -6*E*Iz/L^2,  0,           0,        0,          4*E*Iz/L
    ];
    K_12= [
        -E*A/L,        0,            0,           0,        0,          0;
        0,             -12*E*Iz/L^3, 0,           0,        0,          6*E*Iz/L^2
        0,             0,            -12*E*Iy/L^3,0,       -6*E*Iy/L^2, 0;
        0,             0,            0,          -G*Ix/L,   0,          0;
        0,             0,            6*E*Iy/L^2,  0,        2*E*Iy/L,   0;
        0,             -6*E*Iz/L^2,  0,           0,        0,          2*E*Iz/L
    ];
    K_21= [
        -E*A/L,        0,            0,           0,        0,          0;
        0,             -12*E*Iz/L^3, 0,           0,        0,          -6*E*Iz/L^2
        0,             0,            -12*E*Iy/L^3,0,        6*E*Iy/L^2,0;
        0,             0,            0,          -G*Ix/L,  0,          0;
        0,             0,            -6*E*Iy/L^2, 0,        2*E*Iy/L,   0;
        0,             6*E*Iz/L^2,   0,           0,        0,          2*E*Iz/L
    ];
    
    K_local = [K_node1, K_12;
               K_21, K_node2];
    
end
