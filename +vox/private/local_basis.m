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
    tol = 1e-9;
    % assert( abs(dot(ex,ey))<tol && abs(dot(ey,ez))<tol && abs(dot(ez,ex))<tol );
    % assert( abs(norm(ex)-1)<tol && abs(norm(ey)-1)<tol && abs(norm(ez)-1)<tol );

end