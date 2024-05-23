function [I_app_n,v_n] = n_Iapp(I_0, t,leak_potential,C,R,dt,threshold,reset_potential,sigma_I)
    v_n =zeros(size(t));
    v_n(1) = leak_potential;
    I_app_n = I_0*ones(size(t));
    noise_vec = randn(size(t))*sigma_I*sqrt(dt);
    for i = 2:length(t)
        % Update membrane potential with Forward Euler
        v_n(i) = (v_n(i-1) + dt * ((leak_potential-v_n(i-1)) / R + I_app_n(i)) / C)+noise_vec(i);

        % Check for and reset spike
        if v_n(i) >= threshold
            v_n(i) = reset_potential;
        end
    end
end