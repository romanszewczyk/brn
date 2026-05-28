function r = pkg_ready(t)
  r = 0;              % default: not loaded
  pkg_list = pkg("list");

  % Check if installed
  installed = false;
  for k = 1:length(pkg_list)
    if strcmp(pkg_list{k}.name, t)
      installed = true;
      break
    end
  end

  % Try installing if missing
  if ~installed
    try
      pkg("install", "-forge", t);
      installed = true;
    catch
      installed = false;
    end
  end

  % Try loading
  if installed
    try
      pkg("load", t);
      r = 1;
    catch
      r = 0;
    end
  end
end

