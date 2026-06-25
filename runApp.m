function app = runApp
%RUNAPP Launch the Optimal Filter Sets MATLAB app from the project root.

root = fileparts(mfilename('fullpath'));
appDir = fullfile(root, 'filterSetApp');
if ~isfolder(appDir)
    appDir = fullfile(root, 'FilterSetApp');
end
addpath(appDir);

if nargout > 0
    app = filterDesignerApp(root);
else
    filterDesignerApp(root);
end
end
