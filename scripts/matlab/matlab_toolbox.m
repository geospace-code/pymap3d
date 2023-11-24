function h = matlab_toolbox()

h = struct(mapping=has_mapping(), aerospace=has_aerospace());

end


function has_map = has_mapping()
if usejava('jvm')
  addons = matlab.addons.installedAddons();

  has_map = any(contains(addons.Name, 'Mapping Toolbox'));
else
  has_map = ~isempty(ver("map"));
end
end


function has_map = has_aerospace()
if usejava('jvm')
  addons = matlab.addons.installedAddons();

  has_map = any(contains(addons.Name, 'Aerospace Toolbox'));
else
  has_map = ~isempty(ver("aero"));
end
end
