function vel = dd_vel_profile(V0, Re, R, switching)
    if switching == 0
        vel = vasitas_vortex(V0, Re, R);
    elseif switching == 1
        vel = rankine_vortex(V0, Re, R);
    elseif switching == 2
        vel = burgers_vortex(V0, Re, R);
    elseif switching == 3
        vel = hollow_core_vortex(V0, Re, R);
    else
        error('Invalid switching value for velocity profile.');
    end
end

function vel = rankine_vortex(V0, Re, R)
    vel = zeros(size(R));
    vel(abs(R) <= Re) = V0 * R(abs(R) <= Re) / Re; 
    vel(abs(R) > Re) = V0 * Re ./ R(abs(R) > Re); 
end

function vel = burgers_vortex(V0, Re, R)
    a = 1;
    vt = 1;
    vel = V0 * Re ./ R .* (1 - exp(-a * R.^2 / (4 * vt)));  
end

function vel = vasitas_vortex(V0, Re, R)
    rbar = abs(R) / Re; 
    a = 2;
    n = 2;
    vel = V0 * R ./ (Re * ( (rbar.^(a * n)) + 1 ).^(1 / n)); 
end

function vel = hollow_core_vortex(V0, Re, R)
    vel = zeros(size(R));
    vel(abs(R) > Re) = V0 * Re ./ R(abs(R) > Re);
end


function vel = tan_vel_profile(H,switching)
    g0 = 9.81;
    Ro = 0.375;
    deltaT = 8.0;
    T0 = 305;
    b0 = (g0 * deltaT) / T0;
    K = 0.0;
    if switching == 0 
        K = pi / 2.0;
    elseif switching == 1 
        K = 1.0;
    elseif switching == 2 
        K = 1.702;
    elseif switching == 3 
        K = 0.5;
    else 
        error('Invalid switching value for velocity profile.');
    end 
   coeff = b0 / (K - (Ro^2 / 2.0)); 
   vel = (coeff * H)^0.5; 
end

r = -10:0.01:10;
Re = 7/4;
H = 14.7;

figure
hold on;

for switching = 0:3
    V0 = tan_vel_profile(H, switching);
    plot(r, dd_vel_profile(V0, Re, r, switching));
end
legend('Vatistas','Rankine','Burgers','Hollow-Core')
xlabel('$R$ [m]');
ylabel('$V_\theta$ [ms\textsuperscript{-1}]');
legend show;
grid on;

