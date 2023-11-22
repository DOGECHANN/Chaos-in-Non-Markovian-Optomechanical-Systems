function updateSurface(s, d)
s.ZData(d(1)) = d(2);
drawnow('limitrate');
end