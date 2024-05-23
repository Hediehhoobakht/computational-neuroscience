function [new_vector] = upsampling(old_vector, old_dt, new_dt)
    % Calculate the upsampling factor
    scale_ratio = round(old_dt/new_dt);

    % Calculate the length of the upsampled vector
    length_new = length(old_vector) * scale_ratio;

    % Create the upsampled vector
    new_vector = zeros(1, length_new);

    % Upsample by repeating each value from the old vector
    for k = 1:length(old_vector)
        new_vector((k-1)*scale_ratio+1:k*scale_ratio) = old_vector(k);
    end
end
