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
tube_length = 40/100;
angle = 1/150;
t_joint_diameter = 11.1125/1000;
t_joint_length = 4/100;

%initial conditions%
time = 0;
height = 8/100;
interval = 2;
velocity_out = 0;

while height>=0
    friction_factor_lower = 0.01;
    friction_factor_upper = 0.015;
    while abs(friction_factor_upper-friction_factor_lower)>= 0.0001
        friction_factor_lower=friction_factor_upper;
        %Velocity out of the tube, sans T-joint
        %velocity_out = sqrt((2*gravity*height+2*gravity*tube_length/150)/(1+tube_length*friction_factor_lower/diameter+0.42*(1-diameter^2/(2*width*length/(width+length))^2)));
        %Velocity out of the tube, with T-joint
        velocity_out = sqrt((2*gravity*height+2*gravity*tube_length/150)/(1+tube_length*friction_factor_lower/diameter+0.42*(1-diameter^2/(2*width*length/(width+length))^2)+(1-(diameter/t_joint_diameter))+(friction_factor_lower*t_joint_length)/t_joint_diameter));
        Re = (rho*velocity_out*diameter)/viscosity;
        
        colebrook_eqn = 1/sqrt(f) == -2*log(roughness/(diameter*3.7)+2.51/(Re*sqrt(f)));
        friction_factor_upper = double(solve(colebrook_eqn, f));
    end
    volumetric_flow_rate = velocity_out*(pi*diameter^2)/4;
    
    time = time+interval;
    height = height-((volumetric_flow_rate/(width*length))*interval);
    disp("volumentric_flow_rate: "+volumetric_flow_rate);
    disp("velocity_out: "+velocity_out);
    disp(height + " at time " + time);
end

disp("time to drain: "+time);
