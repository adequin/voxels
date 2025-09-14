function id_arr = init_arr(A)
    %creates struct with voxel ids
    %A is a logical matrix indicating voxel positions
    A = logical(A);
    id_arr = zeros(size(A),'double'); % 
    id_arr(A) = double(1:nnz(A)); % replace true entries with sequential values. nnz-> number of non zero values
end