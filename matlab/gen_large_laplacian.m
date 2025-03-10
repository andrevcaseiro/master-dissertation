num_dims = 3
dim = 2^6
final_dim = dim^num_dims

% [dim], [dim dim], [dim dim dim]
dims = repmat(dim, 1, num_dims);
% {'DD'}, {'DD', 'DD'}, {'DD', 'DD', 'DD'}
bcs = repmat({'DD'}, 1, num_dims);

% Compute the Laplacian matrix
[~, ~, A] = laplacian(dims, bcs);
A = -A;

% Generate filename
laplacian_filename = sprintf('laplacian_%dd_%d.csv', num_dims, final_dim);

export_sparse_mat(A, laplacian_filename);
