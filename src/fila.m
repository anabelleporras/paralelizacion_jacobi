function T = fila(n, version, p, t, it, res, dif)
%FILA Una fila de la tabla de resultados.
T = table(n, string(version), p, t, it, res, dif, 'VariableNames', ...
    {'n', 'version', 'p', 'tiempo_s', 'iteraciones', 'residuo_final', 'dif_vs_seq'});
end