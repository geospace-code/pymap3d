function helper_vdist(lat, lon, N)

  addons = matlab.addons.installedAddons();

  M1 = repmat(lat, N, 1);
  M2 = repmat(lon, N, 1);
  L1 = rand(N,1);
  L2 = rand(N,1);

  if any(addons.Name == "Mapping Toolbox")
    disp("Using Mapping Toolbox distance()")
    f = @() distance(lat, lon, L1, L2);
  elseif ~isempty(what("matmap3d"))
    disp("Using matmap3d.vdist()")
    f = @() matmap3d.vdist(M1, M2, L1, L2);
  else
    error("Matlab Mapping Toolbox is not installed")
  end

  t = timeit(f);
  disp(t)

end
