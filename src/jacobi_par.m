function [x, iter, res_hist] = jacobi_par(A, b, x0, tol, max_iter)
%JACOBI_PAR Método de Jacobi paralelizado por componentes con parfor.
%
%   [x, iter, res_hist] = jacobi_par(A, b, x0, tol, max_iter)
%
%   Mismas entradas y salidas que jacobi_seq. Requiere un pool abierto,
%   preferiblemente de hilos:  parpool("Threads", p)
%
%   Las componentes de x^(k+1) solo dependen de x^(k), así que el ciclo sobre
%   i se reparte entre los workers. El cálculo del residuo y la revisión del
%   criterio de parada se hacen en el cliente, fuera del parfor: esa es la
%   fracción secuencial del algoritmo (ley de Amdahl).

n = length(b);
At = A.';
d = diag(A);
x = x0;
norma_b = norm(b);
res_hist = zeros(max_iter, 1);
iter = 0;

for k = 1:max_iter
    x_nuevo = zeros(n, 1);
    parfor i = 1:n
        % At(:,i), b(i), d(i): variables "sliced" (cada worker recibe su parte)
        % x: variable "broadcast" (todos los workers leen x^(k) completo)
        x_nuevo(i) = (b(i) - At(:, i).' * x + d(i) * x(i)) / d(i);
    end
    x = x_nuevo;

    R = norm(A * x - b) / norma_b;
    iter = k;
    res_hist(k) = R;
    if R < tol
        break;
    end
end

res_hist = res_hist(1:iter);
end