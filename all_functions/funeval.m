% Vikas Pejaver
% August 2015

% Custom function to handle the issue with 'eval' for decision
% trees - based on suggestions from stackoverflow forum

function [varargout] = funeval(func_handle, func_path, varargin)

% Global variables
global CURRDIR;

% Move to function directory
cd(func_path)

[varargout{1:nargout}] = func_handle(varargin{:});

% Return to working directory
cd(CURRDIR);

return
