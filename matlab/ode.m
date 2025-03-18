A = [-2, 1; 1, -1];
x_0 = [1; 0];
b = [2; 0];

expA = expm(A);

res = expA * x_0 + (expA - eye(size(A))) / A * b

