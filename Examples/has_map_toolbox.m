function has_map = has_map_toolbox()
  if usejava('jvm')
    addons = matlab.addons.installedAddons();

    has_map = any(contains(addons.Name, 'Mapping Toolbox'));
  else
    has_map = ~isempty(ver("map"));
  end
end
