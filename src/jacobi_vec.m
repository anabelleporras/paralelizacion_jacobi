function [x, iter, res_hist] = jacobi_vec(A, b, x0, tol, max_iter)
%JACOBI_VEC Método de Jacobi en forma vectorizada (versión de referencia).
%
%   [x, iter, res_hist] = jacobi_vec(A, b, x0, tol, max_iter)
%
%   Mismas entradas y salidas que jacobi_seq. Usa la forma equivalente
%     x^(k+1) = x^(k) + D^{-1} (b - A x^(k)),
%   que evita construir L + U (no crea otra matriz n x n en memoria).

d = diag(A);
x = x0;
norma_b = norm(b);
res_hist = zeros(max_iter, 1);
iter = 0;

for k = 1:max_iter
    x = x + (b - A * x) ./ d;

    R = norm(A * x - b) / norma_b;
    iter = k;
    res_hist(k) = R;
    if R < tol
        break;
    end
end

res_hist = res_hist(1:iter);
end