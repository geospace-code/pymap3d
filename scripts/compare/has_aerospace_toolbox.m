function has_map = has_aerospace_toolbox()
  if usejava('jvm')
    addons = matlab.addons.installedAddons();

    has_map = any(contains(addons.Name, 'Aerospace Toolbox'));
  else
    has_map = ~isempty(ver("aero"));
  end
end
