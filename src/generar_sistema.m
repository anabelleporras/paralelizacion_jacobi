function [A, b, x_exacta] = generar_sistema(n, semilla)
%GENERAR_SISTEMA Construye un sistema Ax = b con dominancia diagonal estricta.
%
%   [A, b, x_exacta] = generar_sistema(n, semilla)
%
%   Entradas:
%     n       - dimensión del sistema
%     semilla - semilla del generador aleatorio (reproducibilidad)
%
%   Salidas:
%     A        - matriz densa n x n; fuera de la diagonal tiene entradas en
%                [0,1) y en la diagonal rand + n, de modo que |a_ii| > sum_{j~=i} |a_ij|
%     b        - lado derecho, construido como b = A * x_exacta
%     x_exacta - solución conocida (vector de unos), útil para validar

rng(semilla);
A = rand(n) + n * eye(n);
x_exacta = ones(n, 1);
b = A * x_exacta;
end