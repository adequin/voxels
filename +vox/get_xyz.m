function coords = get_xyz(id, lat, map)
%returns relative coordinates from voxel center. add lat input to specify
%global coordinates

    if nargin < 3
        if nargin >= 2
            map = coord_map(lat.voxel_size, lat.chamfer);
        else
            map = coord_map(); % cached default
        end
    end

    id_loc = mod(id-1,42)+1; %shift back to base node coords (need to account for mod(42,42)->0)
    if id_loc<=36
        edge  = ceil(id_loc/3);
        slot  = id_loc - 3*(edge-1);
        coords = reshape(map.edge_points(edge,slot,:),1,[]);

    elseif id_loc <=42
        coords = map.face_points(id_loc-36,:);
    else
        error('error mapping id = %d',id_loc);
    end

    if nargin>=2 && ~isempty(lat)
        vox_id = floor((id-1)/42)+1;
        [r,c,s] = ind2sub(size(lat.id), find(lat.id==vox_id));
        offset = (lat.voxel_size + lat.conn_l)*[c-1,s-1,-(r-1)];
        coords = coords+offset;
    end

    % assert(isequal(size(coords), [1 3]))
end