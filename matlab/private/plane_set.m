function planes = plane_set(M)
%PLANE_SET Normal and slip vectors for both nodal planes of every event.
[planes.normal1, planes.slip1] = plane_vectors(M.strike_1, M.dip_1, M.rake_1);
[planes.normal2, planes.slip2] = plane_vectors(M.strike_2, M.dip_2, M.rake_2);
end
