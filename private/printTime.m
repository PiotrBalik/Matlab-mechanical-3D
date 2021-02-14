function t_m=printTime
t_m = toc; % measure time
fprintf('\nTime elapsed: %2.0fm %2.0fs\n', floor(t_m/60), t_m-60*floor(t_m/60) );
end
