function [new_vector] = expandbin(old_vector, old_dt, new_dt)
length_old = length(old_vector);
scale_ratio = round(new_dt/old_dt);
length_new = round(length_old/scale_ratio);
new_vector = zeros(1,length_new);
tsteps = 50;
for k = 1:length(new_vector)
new_vector(k) = mean(old_vector((k-1)*tsteps+1:k*tsteps));
end
