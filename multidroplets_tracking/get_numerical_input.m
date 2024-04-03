function output=get_numerical_input(message)
% This function asks for a numerical value from uses with message
% It outputs a single numerical value as given by the user.

% Assert message is a character array
if ischar(message)
    prompt = message;
else
    prompt = char(message);
end

% Accept any type of user input
try 
    user_input=input(prompt);
catch
end
output=user_input;