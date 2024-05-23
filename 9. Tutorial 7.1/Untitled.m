close all
clear
clc

r_max=100;
theta_E=-5;
theta_I=0;
alpha_E=0.05;
alpha_I=1;
W_EE=2;
W_EI=2.5;
W_IE=-2.5;
W_II=-2;
t_max=3;
r_E(1)=50;
r_I(1)=50;
dt=0.1e-3;
t=0:dt:t_max;
Istim = zeros(size(t));
ind = find (t > 1 & t < 2);
Istim(ind) = 20;
Ibase_E = 0;
Ibase_I = 0;
tau_E=5e-3;
tau_I=5e-3;
I_E_app = ones(size(t))*Ibase_E;
I_I_app = ones(size(t))*Ibase_I+Istim;
% I_E(1)=0
% I_I(1)=0
for n = 2:length(t)
    I_E(n)=W_EE*r_E(n-1)+W_IE*r_I(n-1)+I_E_app(n-1);
    I_I(n)=W_EI*r_E(n-1)+W_II*r_I(n-1)+I_I_app(n-1);
    r_E(n) =  r_E(n-1)+(dt/tau_E)*(-r_E(n-1)+alpha_E*((I_E(n)-theta_E)^2).*sign(I_E(n)-theta_E));
    allowed_vals = find((r_E >= 0 ).*( r_E <= r_max ));
    

    r_I(n) =  r_I(n-1)+(dt/tau_I)*(-r_I(n-1)+alpha_I*((I_I(n)-theta_I)));
    allowed_vals = find((r_I>= 0 ).*( r_I <= r_max ));
    

    
end

f1=figure(1);
subplot(3,1,1);
 plot(t(allowed_vals),r_E(allowed_vals),'k:');
%plot(t,r_E)
subplot(3,1,2);
plot(t,r_I)
subplot(3,1,3);
plot(r_E,r_I)
I_E(1)=0;
I_I(1)=0;
I_E_base(1)=25;
I_I_base(1)=15;
r_E(1)=0;
r_I(1)=0;
for n = 2:length(t)
    I_E_app(n)=I_E_base;
    if t(n)>1 && t(n)<=2

            I_I_app(n)=I_I_base+20;
    else

            I_I_app(n)=I_I_base;
     end

    I_E(n)=W_EE*r_E(n-1)+W_IE*r_I(n-1)+I_E_app(n-1);
    I_I(n)=W_EI*r_E(n-1)+W_II*r_I(n-1)+I_I_app(n-1);
    r_E(n) =  r_E(n-1)+(dt/tau_E)*(-r_E(n-1)+alpha_E*((I_E(n)-theta_E)^2).*sign(I_E(n)-theta_E));
    
    if r_E(n) < 0
            r_E(n) = 0; % If r_E_temp is negative, set r_E to 0
    elseif r_E(n) > r_max
        r_E(n) = r_max; % If r_E_temp exceeds 10, set r_E to 10
    else
        r_E(n) = r_E(n); % Otherwise, keep the calculated value
    end

    r_I(n) =  r_I(n-1)+(dt/tau_I)*(-r_I(n-1)+alpha_I*((I_I(n)-theta_I)));
    if r_I(n) < 0
            r_I(n) = 0; % If r_E_temp is negative, set r_E to 0
    elseif r_I(n) > r_max
        r_I(n) = r_max; % If r_E_temp exceeds 10, set r_E to 10
    else
        r_I(n) = r_I(n); % Otherwise, keep the calculated value
    end

    
end

f2=figure(2);
subplot(3,1,1);
plot(t,r_E)
subplot(3,1,2);
plot(t,r_I)
subplot(3,1,3);
plot(r_E,r_I)
% hold on
% plot(t,r_2)
% xlabel('time');
% ylabel('Firing rate');
% title('part 1');
% legend('r1','r2')
% saveas(f1, sprintf('1.png'));

tau_E=2e-3;
tau_I=10e-3;
dt=0.1e-3;
t=0:dt:t_max;
I_E(1)=0;
I_I(1)=0;
I_E_base(1)=0;
I_I_base(1)=0;
r_E(1)=0;
r_I(1)=0;
for n = 2:length(t)
    I_E_app(n)=I_E_base;
    if t(n)>1 && t(n)<=2

            I_I_app(n)=I_I_base+20;
    else

            I_I_app(n)=I_I_base;
    end
    I_E(n)=W_EE*r_E(n-1)+W_IE*r_I(n-1)+I_E_app(n-1);
    I_I(n)=W_EI*r_E(n-1)+W_II*r_I(n-1)+I_I_app(n-1);
    r_E(n) =  r_E(n-1)+(dt/tau_E)*(-r_E(n-1)+alpha_E*((I_E(n)-theta_E)^2).*sign(I_E(n)-theta_E));
    if r_E(n) < 0
            r_E(n) = 0; % If r_E_temp is negative, set r_E to 0
    elseif r_E(n) > r_max
        r_E(n) = r_max; % If r_E_temp exceeds 10, set r_E to 10
    else
        r_E(n) = r_E(n); % Otherwise, keep the calculated value
    end

    r_I(n) =  r_I(n-1)+(dt/tau_I)*(-r_I(n-1)+alpha_I*((I_I(n)-theta_I)));
    if r_I(n) < 0
            r_I(n) = 0; % If r_E_temp is negative, set r_E to 0
    elseif r_I(n) > r_max
        r_I(n) = r_max; % If r_E_temp exceeds 10, set r_E to 10
    else
        r_I(n) = r_I(n); % Otherwise, keep the calculated value
    end

end

f3=figure(3);
subplot(3,1,1);
plot(t,r_E)
subplot(3,1,2);
plot(t,r_I)
subplot(3,1,3);
plot(r_E,r_I)
I_E(1)=0;
I_I(1)=0;
I_E_base(1)=25;
I_I_base(1)=15;
r_E(1)=0;
r_I(1)=0;
for n = 2:length(t)
    I_E_app(n)=I_E_base;
    if t(n)>1 && t(n)<=2

            I_I_app(n)=I_I_base+20;
    else

            I_I_app(n)=I_I_base;
     end

    r_E(n) =  r_E(n-1)+(dt/tau_E)*(-r_E(n-1)+alpha_E*((I_E(n-1)-theta_E)^2).*sign(I_E(n-1)-theta_E));
    if r_E(n) < 0
            r_E(n) = 0; % If r_E_temp is negative, set r_E to 0
    elseif r_E(n) > r_max
        r_E(n) = r_max; % If r_E_temp exceeds 10, set r_E to 10
    else
        r_E(n) = r_E(n); % Otherwise, keep the calculated value
    end

    r_I(n) =  r_I(n-1)+(dt/tau_I)*(-r_I(n-1)+alpha_I*((I_I(n-1)-theta_I)));
    if r_I(n) < 0
            r_I(n) = 0; % If r_E_temp is negative, set r_E to 0
    elseif r_I(n) > r_max
        r_I(n) = r_max; % If r_E_temp exceeds 10, set r_E to 10
    else
        r_I(n) = r_I(n); % Otherwise, keep the calculated value
    end

    I_E(n)=W_EE*r_E(n)+W_IE*r_I(n)+I_E_app(n);
    I_I(n)=W_EI*r_E(n)+W_II*r_I(n)+I_I_app(n);
end

f4=figure(4);
subplot(3,1,1);
plot(t,r_E)
subplot(3,1,2);
plot(t,r_I)
subplot(3,1,3);
plot(r_E,r_I)