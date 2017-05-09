function varargout = my_call(value)

% fprintf('%d ',value)
fprintf('You clicked: %s \t', get(gcbo,'String'))  % "gcbo" is the handle of the object whose callback is executing