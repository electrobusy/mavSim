% ---
% Function used to bound a certain value within certain limits - acts as a
% "saturation function" 
% ---
function value = bound(value, min_value, max_value)

if (value > max_value)
    value = max_value;
elseif (value < min_value)
    value = min_value;
end

end