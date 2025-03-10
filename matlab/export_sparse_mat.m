function export_sparse_mat(A, filename)
    % Get the nonzero elements and their indices
    [rows, cols, values] = find(A);
    
    % Open the file for writing
    f = fopen(filename, 'w');
    
    % Write number of nonzero elements
    fprintf(f, '%d\n', nnz(A));
    
    % Write triplets
    for i = 1:length(values)
        fprintf(f, '%d,%d,%g\n', rows(i), cols(i), values(i));
    end
    
    % Close the file
    fclose(f);
end
