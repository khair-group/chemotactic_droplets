function output=get_textual_input(message,acceptables)
% This function asks for textual input from user with 'message'
% It outputs a textual value as given by the user.
% If a list of acceptable inputs is applied, it keeps asking until the
% answer is provided is in the list of acceptables
% Acceptables gives a list of values that are considered acceptable

if nargin>0
    % Assert message is a character array
    if ischar(message)
        prompt = message;
    else
        prompt = char(message);
    end
end


if nargin>1 && ~isempty(acceptables)
    % Accept only textual user input given in list acceptables
    while 1==1
        try
            user_input=input(prompt,'s');
        catch
        end
        if contains(acceptables,user_input)
            valid_input=user_input;
            break
        else
            disp(['Please only provide response with a character in the following list: ' acceptables '.'])
        end
    end
else
    valid_input=input(prompt, 's');
end

output=string(valid_input);