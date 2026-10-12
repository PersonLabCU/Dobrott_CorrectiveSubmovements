%% Function to Adjust Length Mismatch
function output = enforce_length(data, expected_length)
    actual_length = size(data, 2);
    
    if actual_length < expected_length
        warning("Data length (%d) shorter than expected length (%d). Padding with NaNs.", actual_length, expected_length);
        output = [data, nan(1, expected_length - actual_length)];
    elseif actual_length > expected_length
        warning("Data length (%d) longer than expected length (%d). Trimming data.", actual_length, expected_length);
        output = data(:, 1:expected_length);
    else
        output = data;
    end
end