function string_x_u = create_x_and_u_string(x,u)
% Given a state-input pair, this function produces the string with the name
% of the pair. This is used in the package to identify a state-input pair.

% If the second argument is empty, the function returns a list with x
% elements only.

temp_x = compose('%.4f',x);

if ~isempty(u)
    temp_u = compose('%.4f',u);
    string_x_u = strcat('X:',temp_x{:},'|','U:',temp_u{:});
else % If 2nd arg is empty, create the string for the x element only
    string_x_u = strcat('X:',temp_x{:});
end

end

