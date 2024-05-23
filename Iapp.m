function [I_app,v] = Iapp(I_0, t,leak_potential,C,R,dt,threshold,reset_potential)
    v =zeros(size(t));
    v(1) = leak_potential;
    I_app = I_0*ones(size(t));
    for i = 2:length(t)
        % Update membrane potential with Forward Euler
        v(i) = v(i-1) + dt * ((leak_potential-v(i-1)) / R + I_app(i)) / C;

        % Check for and reset spike
        if v(i) >= threshold
            v(i) = reset_potential;
        end
    end
end