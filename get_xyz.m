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

    % assert(isequal(size(coords), [1 3]))
end