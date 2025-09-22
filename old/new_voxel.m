% Voxel simplified model stiffness matrix
% Notes: 
% solve symmetry issue, probably some rotated beams
% individual dof constraints
% convert strain gauge data to force inputs
% show constrained nodes and strain gauge data in visualizer
% add 


% Validation 2D tests
% - [o]  test single non rotated beam in 2D
% - [o]  Test single rotated beam in 2D -> rotations cause small errors on
% the order of e-16
% - [o]  test parallel non rotated beam in 2D
% - [o]  test parallel non rotated beam of different lengths in 2D
% - [o]  test single pi/4 rotated beam in 2D
% - [o]  test parallel rotated beam in 2D - if both fixed, doubles
% stiffness, if one fixed, translations caused by node 1 rotation are seen
% - [o]  test combined rotated + non rotated beam in 2D
% - [o]  test serial non rotated beam in 2D
% - [o]  test serial rotated beam in 2D

%Validation 3D tests
% - [o]  test single rotated beam in 3D (y-axis): figured out an issue with
% the rotation matrices, was rotating the wrong way
% - []  test single doubly rotated beam in 3D (y,z-axis) - ok this actually
% matters.
% - [o]  test single non rotated beam with 3D force vector - validated with
% same relative force applied to doubly rotated beam and then converted displacements 
% - [o]  test parallel doubly rotated beam in 3D (just validating ratios 
% between two 45deg beams and single axial beam; should be a factor of 2)

%%% implemented x rotations and fliped rotation orders to be xyz %%%
% - [o]  test parallel structure with a beam separating them - reorganized
% force matrix to be like displacement output
% - [o]  test serial doubly rotated beam in 3D -> rotation order is applied
% in reverse direction z->y->x global frame

%% debugging lack of symmetry
% - adding vector display -> all green arrows point in the right direction
% - find y stiffness values for each face ->
% - main issue was still rotation order and direction, vector display does
% not necessarily match stiffness matrix

%% 
% clear; clc;
b = 1.6e-3;% width in m
h = 6e-3; % height in m

e = 25e9; %youngs modulus (25-30 GPA)in Pa
g = 10e9; %shear modulus estimate in Pa
a = b*h; % in m^2

iy = h*b^3/12;% area moment of inertia about y
iz = b*h^3/12; % area moment of inertia about z
ix = iy+iz;% geometric value, but try the empirical version if results are off
l = 47.686e-3; %m
l_45 = l*sqrt(2);
F_value = 1000; %Define load condition: 1 Newton in y1 direction

%%name nodes%%
%face 1 (front face)
b1 = 1; %bottom
cl1 = 2; %corner left
m1 = 3; %middle
cr1 = 4; %corner right
t1 = 5; %top

%face 2 (right face)
b2 = 6; %bottom
cf2 = cr1; %corner front
m2 = 7; %middle
cb2 = 8; %corner back
t2 = 9; %top

%face 3 (back face)
b3 = 10; %bottom
cl3 = 11; %corner left (viewed from the front)
m3 = 12; %middle
cr3 = cb2; %corner right
t3 = 13; %top

%face 4 (left face)
b4 = 14; %bottom
cf4 = cl1; %corner front
m4 = 15; %middle
cb4 = cl3; %corner back
t4 = 16; %top

%face 5 (top face)
f5 = t1; %front
l5 = t4; %left
m5 = 17; %middle
r5 = t2; %right
bk5 = t3; %back

%face 6 (bottom face)
f6 = b1; %front
l6 = b4; %left
m6 = 18; %middle
r6 = b2; %right
bk6 = b3; %back

offset = pi/2;
% rotations happen right to left
beam_config = [
   % all faces
   %face 1
                b1,cl1,offset,0,3*pi/4, l_45
                b1,m1,offset,0,pi/2, l
                b1,cr1,offset,0,pi/4,l_45
                cl1,t1,offset,0,pi/4,l_45
                m1,t1,offset,0,pi/2, l
                cr1,t1,offset,0,3*pi/4, l_45
                cl1,m1,offset,0,0, l
                m1,cr1,offset,0,0, l
    % %face 2
                cf2,b2,offset,pi/2,-pi/4,l_45
                cf2,m2,offset,pi/2,0, l
                cf2,t2,offset,pi/2,pi/4,l_45

                b2,cb2,offset,pi/2,pi/4,l_45
                m2,cb2,offset,pi/2,0, l
                t2,cb2,offset,pi/2,-pi/4,l_45

                b2,m2,offset,pi/2,pi/2, l
                m2,t2,offset,pi/2,pi/2, l
   %            
   % %face 3
                cl3,b3,offset,0,-pi/4,l_45
                cl3,m3,offset,0,0,l
                cl3,t3,offset,0,pi/4,l_45

                b3,cr3,offset,0,pi/4,l_45
                m3,cr3,offset,0,0,l
                t3,cr3,offset,0,-pi/4,l_45

                b3,m3,offset,0,pi/2,l
                m3,t3,offset,0,pi/2,l

   %   %face 4
                cf4,b4,offset,pi/2,-pi/4,l_45
                cf4,m4,offset,pi/2,0,l
                cf4,t4,offset,pi/2,pi/4,l_45

                b4,cb4,offset,pi/2,pi/4,l_45
                m4,cb4,offset,pi/2,0,l
                t4,cb4,offset,pi/2,-pi/4,l_45

                b4,m4,offset,pi/2,pi/2,l
                m4,t4,offset,pi/2,pi/2,l

      %face 5 (top)

                f5,l5,-pi/2+offset,0,3*pi/4,l_45
                f5,m5,-pi/2+offset,0,pi/2,l
                f5,r5,-pi/2+offset,0,pi/4,l_45

                l5,bk5,-pi/2+offset,0,pi/4,l_45
                m5,bk5,-pi/2+offset,0,pi/2,l
                r5,bk5,-pi/2+offset,0,3*pi/4,l_45

                l5,m5,-pi/2+offset,0,0,l
                m5,r5,-pi/2+offset,0,0,l
        %face 6 (bottom)
                f6,l6,-pi/2+offset,0,3*pi/4,l_45
                f6,m6,-pi/2+offset,0,pi/2,l
                f6,r6,-pi/2+offset,0,pi/4,l_45

                l6,bk6,-pi/2+offset,0,pi/4,l_45
                m6,bk6,-pi/2+offset,0,pi/2,l
                r6,bk6,-pi/2+offset,0,3*pi/4,l_45

                l6,m6,-pi/2+offset,0,0,l
                m6,r6,-pi/2+offset,0,0,l
    ]; 

% beam_config1 = [1,2,0,pi/2,pi/2,l];
% beam_config2 = [1,2,0,0,pi/2,l];

%compute
K_actual = assemble_stiffness_3D(beam_config, e, g, a, ix, iy, iz);% Define stiffness matrix

%Define applied Force
F_input = zeros(length(K_actual)/6,6); %index using F(node, dofs)
F_input(t1,1:3) = [0,0,-F_value];
F_input(t2,1:3) = [0,0,-F_value];
F_input(t3,1:3) = [0,0,-F_value];
F_input(t4,1:3) = [0,0,-F_value];
% F_input(m5,1:3) = [0,-F_value,0];

%Define constraints
constrained_nodes = [b1,b2,b3,b4];

Displacements = find_displacements(K_actual, F_input, constrained_nodes);

% % -- after you have beam_config and Displacements -----------------
coords = beamCoords(beam_config);
plotBeamDeform(beam_config, coords, Displacements, 20);  % scale=50×
% plotBeamDeform(beam_config, coords, zeros(size(coords,1),6), 1); % no deformation
% showBeamAxes(beam_config, coords, 0.12); 

%find individual stiffnesses -> K = F/x


% reference axial displacement
d = F_value*l/(e*a);

function x_full = find_displacements(K_actual, F, constrained_nodes)
    % finds numerical displacement values
    % Inputs:
    % k_assembly - 3D Stiffness matrix
    % loads - Force vector with applied loads
    % constraints - fixed degrees of freedom
    % Outputs: 
    % 6 by n matrix where each column is a node
    % index using (node, dof)
    F = reshape(F',[],1);
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
% 3D General Stiffness Matrix Assembly
function K_global = assemble_stiffness_3D(beam_data_3D, E, G, A, Ix, Iy, Iz)
    % Inputs:
    % num_nodes - Total number of nodes in the structure
    % beam_data - Matrix where each row is [node1, node2, angles (radians)]
    % E - Young's modulus
    % G - Shear modulus
    % A - Cross-sectional area
    % Ix- Polar second moment of area
    % Iy - Second moment of area about y
    % Iz - Second moment of area about z

    % Initialize global stiffness matrix
    num_nodes = max(max(beam_data_3D(:,1:2)));
    K_global = zeros(6 * num_nodes);

    % Loop through each beam
    % disp(size(beam_data_3D, 1))
    for i = 1:size(beam_data_3D, 1)
        % Extract beam data
        node1 = beam_data_3D(i, 1);
        node2 = beam_data_3D(i, 2);
        anglex = beam_data_3D(i, 3);
        angley = beam_data_3D(i, 4);
        anglez = beam_data_3D(i, 5);
        L      = beam_data_3D(i, 6);
        % Compute local stiffness matrix for the beam
        K_local = beam_stiffness_matrix_3D(E, G, A, Ix, Iy, Iz, L);

        % Compute rotation matrix
        R = rotation_matrix_3D(anglex, angley, anglez);

        % Transform local stiffness matrix to global coordinates
        K_global_beam = R * K_local * R';

        % Assemble into the global stiffness matrix
        dof_indices = [6*node1-5:6*node1, 6*node2-5:6*node2];
        for ii = 1:12
            for jj = 1:12
                K_global(dof_indices(ii), dof_indices(jj)) = ...
                    K_global(dof_indices(ii), dof_indices(jj)) + K_global_beam(ii, jj);
            end
        end
    end
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

function R = rotation_matrix_3D(anglex, angley, anglez)
    % Returns the 12x12 rotation matrix for a beam rotated by a an angle
    % around z and y
    % use the negative because we are rotating the beam vectors from local to
    % global, not actually rotating the beam itself, which is in the
    % opposite direction to the beam rotation
    cz = cos(anglez);
    sz = sin(anglez);

    % Rotation matrix for 3D with moments
    Rz = [
        cz, -sz,  0;
        sz,  cz,  0;
         0,   0,  1;
   
    ];
    
    cy = cos(angley);
    sy = sin(angley);

    Ry = [
        cy,  0,  sy;
        0,   1,   0;
       -sy,  0,  cy 
    ];

    cx = cos(anglex);
    sx = sin(anglex);

    Rx = [
        1,   0,   0;
        0,   cx,  -sx;
        0,   sx,  cx
    ];

    Rxyz = Rx*Ry*Rz; %the order is flipped since matrix operations are right to left
    R = blkdiag(Rxyz, Rxyz, Rxyz, Rxyz);
    % e1 = -Rxyz' * [1;0;0];  % local x
    % e2 = -Rxyz' * [0;1;0];  % local y
    % e3 = -Rxyz' * [0;0;1];  % local z
    % disp('Local axes in global frame:');
    % disp([e1, e2, e3]);
end

% Display beam_config
function nodeCoords = beamCoords(beam_config)
% Returns an array  [x y z]  for every node index that appears.
%
% beam_config cols: [node1 node2 anglex angley anglez L]

    ex = [1;0;0];                           % local x unit
    maxNode = max( beam_config(:,1:2) ,[],'all');
    nodeCoords = nan(maxNode,3);            % pre‑fill NaNs
    nodeCoords(1,:) = [0 0 0];              % start at origin

    unresolved = true;
    while unresolved
        unresolved = false;
        for i = 1:size(beam_config,1)
            n1 = beam_config(i,1);      n2 = beam_config(i,2);
            ax = beam_config(i,3);      ay = beam_config(i,4);
            az = beam_config(i,5);      L  = beam_config(i,6);

            % rotation matrix (intrinsic X‑Y‑Z, matches your code)
            Rx = [1 0 0; 0 cos(ax) -sin(ax); 0 sin(ax) cos(ax)];
            Ry = [cos(ay) 0 sin(ay); 0 1 0; -sin(ay) 0 cos(ay)];
            Rz = [cos(az) -sin(az) 0; sin(az) cos(az) 0; 0 0 1];
            R  = Rx*Ry*Rz;
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

function plotBeamDeform(beam_config, nodeCoords, U, scale)
% Plots undeformed (gray) and deformed (blue) beams.
% U : N×6 displacement matrix, translational DOFs in columns 1..3

    if nargin < 4, scale = 1; end
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


