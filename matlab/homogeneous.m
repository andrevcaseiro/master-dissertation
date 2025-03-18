[~, ~, A] = laplacian([4 4 4], {'DD', 'DD', 'DD'});
A = -A;

expA = expm(A);

x_0 = zeros(size(A,1), 1);
x_0(1) = 1;

b = zeros(size(A,1), 1);
b(1) = 1;

res = (expA - eye(size(A))) * inv(A)
res = expA * x_0 + res * b
