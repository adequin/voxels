function u = estimate_with_sensors(lat, K, beam_n1, beam_n2, readings)
% K        : (ndof x ndof) global stiffness
% beam_n1/2     : (nb x 1) node IDs of gauged vertical members
% readings    : (nb x 1) measured axial strains on those members (need to calculate or calibrate from resistor readings)
    
    [bc_idx, bc_val, ~] = estimate_topface(lat, beam_n1, beam_n2, readings);
    assert(numel(unique(bc_idx)) == numel(bc_idx),"duplicates");
    [u, ~] = solve_with_dirichlet(K, bc_idx, bc_val);
end

function [bc_idx, bc_val, fit] = estimate_topface(lat, beam_n1, beam_n2, readings, weights)

    % beam_n1/2     : (nb x 1) node IDs of gauged vertical members
    % readings    : (nb x 1) measured axial strains on those members
    % weights       : (optional nb x 1) confidence weights for each gauge
    
    % RETURNS
    % bc_idx        : DOF indices (global) to constrain (z DOFs of top face)
    % bc_val        : prescribed values for those DOFs (negative = compression)
    % fit           : struct with fields: dz, thx, thy, rms
     
    nb = length(beam_n1);
    e = readings;

    assert(length(beam_n2)==nb & length(readings)==nb);
    assert(size(lat.A,2)==1 & size(lat.A,3)==1,"assembly should be vertical");

    max_height = max(sum(lat.A,1)); %sums along rows, assume assembly has no gaps
    H = max_height * lat.voxel_size + (max_height-1)*lat.conn_l; %finds total height

    p1 = zeros(nb,3);
    p2 = zeros(nb,3);
    pmid = zeros(nb,3);
    Li = H*ones(nb,1);  % structure height, assumes strain is the same everywhere vertically.
    
    for i = 1:nb
        p1(i,:) = get_xyz(beam_n1(i));
        p2(i,:) = get_xyz(beam_n2(i));
        pmid(i,:) = (p1(i,:)+p2(i,:))*0.5;
        % Li(i) = norm(p2(i,:) - p1(i,:));
    end

    xi = pmid(:,1);
    yi = pmid(:,2);

    % Design matrix A * p = eps, where p = [dz; thx; thy]
    % eps_i ≈ (dz + thx*yi - thy*xi) / Li
    A = [ ones(nb,1),  yi,     -xi ];
    A = A ./ Li;                % divide each column by Li elementwise

    %weight stuff
    if nargin >= 5 && ~isempty(weights)
        W = diag(weights(:));
        Aw = W^(1/2) * A;
        ew = W^(1/2) * e;
    else
        Aw = A; ew = e;
    end

    ew = reshape(ew, [],1);
    p = Aw\ew;
    dz  = p(1);
    thx = p(2);
    thy = p(3);

    % Build top-face displacement BCs: uz(x,y) = (dz + thx*y - thy*x)
    n_top = define_nodes(lat,1,1,1); %assumes nx1x1 column
    n_bot = define_nodes(lat,size(lat.A,1),1,1); %assumes nx1x1 column
    top = [];
    bot = []; 
    XY = [];

    for name = n_bot.face_nodes
        bot = [bot, n_bot.f6.(name)];
        top = [top, n_top.f5.(name)]; 
        XY = [XY; get_xyz(n_top.f5.(name))]; %we dont care about z so this should work
    end

    % ux  = zeros(size(top)); % not used here
    % uy  = zeros(size(top)); % not used here
    uz  = (dz + thx*XY(:,2) - thy*XY(:,1));

    % Map node IDs to global DOF indices: assuming DOF order [ux, uy, uz, rx, ry, rz] per node.
    % If your ordering differs, change the stride and offset.
    dof_per_node = 6;
    uz_offset    = 3;  % ux=1, uy=2, uz=3 in this convention
    top_z = (top-1)*dof_per_node + uz_offset;
    bot_all = [];
    for node = bot
        bot_all = [bot_all, 6*node-5:6*node];
    end
    bc_idx = [top_z, bot_all];
    bc_val = zeros(numel(bc_idx),1);
    bc_val(1:numel(top)) = uz;

    % Diagnostics
    r    = (A * (size(A,2)==1)*dz + (size(A,2)~=1)*([ones(nb,1) yi -xi]*[dz; thx; thy])./Li) - e;
    fit.dz   = dz;  fit.thx  = thx;   fit.thy = thy;
    fit.rms  = sqrt(mean(r.^2));
end