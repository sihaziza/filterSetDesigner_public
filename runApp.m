function app = runApp
%RUNAPP Launch the Optimal Filter Sets MATLAB app from the project root.

root = fileparts(mfilename('fullpath'));
appDir = fullfile(root, 'FilterSetApp');
addpath(appDir);

if nargout > 0
    app = OptimalFilterApp(root);
else
    OptimalFilterApp(root);
end
end
