function [t, x, it, res] = medir(f, reps)
%MEDIR Ejecuta f() reps veces y devuelve la mediana del tiempo.
tiempos = zeros(reps, 1);
for r = 1:reps
    tic;
    [x, it, res] = f();
    tiempos(r) = toc;
end
t = median(tiempos);
end