for num_dims = 2:2
    dim = 256;
    final_dim = dim^num_dims
    while final_dim <= 32768
        % [dim], [dim dim], [dim dim dim]
        dims = repmat(dim, 1, num_dims);
        % {'DD'}, {'DD', 'DD'}, {'DD', 'DD', 'DD'}
        bcs = repmat({'DD'}, 1, num_dims);

        % Compute the Laplacian matrix
        [~, ~, A] = laplacian(dims, bcs);
        A = -A;
        A_exp = expm(A);

        first_row_exp = A_exp(1, :);

        % Generate filenames
        laplacian_filename = sprintf('laplacian_%dd_%d.csv', num_dims, final_dim);
        exp_filename = sprintf('laplacian_%dd_%d_exp.csv', num_dims, final_dim);

        % Save matrices
        export_sparse_mat(A, laplacian_filename);
        export_sparse_mat(first_row_exp, exp_filename);

        dim = dim * 2;
        final_dim = dim^num_dims
    end
end
