clear
syms f
%constants%
gravity = 9.81;
rho = 998;
diameter = 7.94/1000;
viscosity = 1.0016/1000;
roughness = 0.0015/1000;
length = 32/100;
width = 26/100;
tube_length = 20/100;
angle = 1/150;
t_joint_diameter = 11.1125/1000;
t_joint_length = 4/100;

%initial conditions%
time = 0;
height = 10/100;
interval = 1;
velocity_out = 0;

%initialize arrays for plots%
velocity_plot = [];
mass_flow_plot = [];
Re_plot = [];
ff_plot = [];
height_plot = [];
time_plot = [];

while height>=0.02
    friction_factor_lower = 0.01;
    friction_factor_upper = 0.015;
    while abs(friction_factor_upper-friction_factor_lower)>= 0.0001
        friction_factor_lower=friction_factor_upper;
        
        %Velocity out of the tube, sans T-joint
        velocity_out = sqrt((2*gravity*height+2*gravity*tube_length/150)/(1+tube_length*friction_factor_lower/diameter+0.5));
        
        %Velocity out of the tube, with T-joint
%         velocity_out = sqrt((2*gravity*height+2*gravity*tube_length/150)/(1+tube_length*friction_factor_lower/diameter+0.5+(1-(diameter/t_joint_diameter))+(friction_factor_lower*t_joint_length)/t_joint_diameter));
                
        Re = (rho*velocity_out*diameter)/viscosity;
        
        if Re > 2300
            friction_factor_eqn = 1/sqrt(f) == -2*log10(roughness/(diameter*3.7)+2.51/(Re*sqrt(f)));
            friction_factor_upper = double(solve(friction_factor_eqn, f));
        else
            friction_factor_eqn = f == 64/Re;
            friction_factor_upper = double(solve(friction_factor_eqn, f));
        end
    end
    volumetric_flow_rate = velocity_out*(pi*diameter^2)/4;
    mass_flow = rho*(pi*(diameter/2)^2)*velocity_out;
    
    time = time+interval;
    height = height-((volumetric_flow_rate/(width*length))*interval);
    disp("volumentric_flow_rate: "+volumetric_flow_rate);
    disp("mass_flow: "+mass_flow);
    disp("velocity_out: "+velocity_out);
    disp(height + " at time " + time);
    disp("Re: " + Re);
    disp("friction_factor: " + friction_factor_upper)
    
    %update plot arrays
    velocity_plot = [velocity_plot, velocity_out];
    mass_flow_plot = [mass_flow_plot, mass_flow];
    Re_plot = [Re_plot, Re];
    ff_plot = [ff_plot, friction_factor_lower];
    height_plot = [height_plot, height];
    time_plot = [time_plot, time];
    
end

%Plots%
disp("time to drain: "+time);
figure();
plot(time_plot, velocity_plot);
title("Velocity vs Time")
xlabel("Time (s)")
ylabel("Velocity (m/s)")

figure();
plot(time_plot, Re_plot);
title("Re vs Time")
xlabel("Time (s)")
ylabel("Re")

figure();
plot(time_plot, mass_flow_plot);
title("Mass Flow Rate vs Time")
xlabel("Time (s)")
ylabel("Mass Flow Rate (kg/s)")

figure();
plot(time_plot, height_plot);
title("Height of Water vs Time")
xlabel("Time (s)")
ylabel("Height of Water (m)")

figure();
plot(ff_plot, velocity_plot);
title("Velocity vs Friction Factor")
xlabel("Friction Factor")
ylabel("Velocity (m/s)")
