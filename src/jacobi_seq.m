function [x, iter, res_hist] = jacobi_seq(A, b, x0, tol, max_iter)
%JACOBI_SEQ Método de Jacobi secuencial, con un ciclo sobre las componentes.
%
%   [x, iter, res_hist] = jacobi_seq(A, b, x0, tol, max_iter)
%
%   Entradas:
%     A        - matriz n x n con diagonal no nula
%     b        - lado derecho (n x 1)
%     x0       - aproximación inicial (n x 1)
%     tol      - tolerancia sobre el residuo relativo R_k = ||Ax - b||_2 / ||b||_2
%     max_iter - número máximo de iteraciones
%
%   Salidas:
%     x        - última aproximación
%     iter     - iteraciones realizadas
%     res_hist - historial de R_k (iter x 1)
%
%   Cada componente se calcula como
%     x_i^(k+1) = ( b_i - sum_{j~=i} a_ij x_j^(k) ) / a_ii
%   usando el producto punto de la columna i de A.' (memoria contigua en
%   MATLAB). Se usa el MISMO núcleo de cálculo que en jacobi_par, para que la
%   comparación solo mida el efecto de paralelizar.

n = length(b);
At = A.';                 % columnas de At = filas de A (acceso contiguo)
d = diag(A);
x = x0;
x_nuevo = zeros(n, 1);
norma_b = norm(b);
res_hist = zeros(max_iter, 1);
iter = 0;

for k = 1:max_iter
    for i = 1:n
        x_nuevo(i) = (b(i) - At(:, i).' * x + d(i) * x(i)) / d(i);
    end
    x = x_nuevo;

    % Criterio de parada: idéntico en las tres versiones
    R = norm(A * x - b) / norma_b;
    iter = k;
    res_hist(k) = R;
    if R < tol
        break;
    end
end

res_hist = res_hist(1:iter);
end