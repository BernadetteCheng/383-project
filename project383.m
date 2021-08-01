%constants%
gravity = 9.81;
rho = 998;
diameter = 7.94/1000;
viscosity = 1.0016e-3;
roughness = 0.0015/1000;
length = 32/1000;
width = 26/1000;
tube_length = 10/100;
angle = 1/150;

%initial conditions%
time = 0;
height = 8/100;
friction_factor_upper = 0.01;
friction_factor_lower = 0.02;
interval = 1;

while height>=0
    velocity_out = 0;
    while (friction_factor_lower-friction_factor_upper)>= 0.0001
        velocity_out = sqrt((2*gravity*height+2*gravity*tube_length/150)/(1+tube_length*ffactor/diameter+0.42*(1-diameter^2/(2*width*length/(width+length)^2))));
        Re = (rho*velocity_out*diameter)/viscosity;
        
        colebrook_eqn = 1/sqrt(f) == -2*log(roughness/(diameter*3.7)+2.51/(Re*sqrt(f)));
        friction_factor = double(solve(eqn, f));
    end
    volumetric_flow_rate = velocity_out*(pi*diameter^2)/4;
    
    time = time+interval;
    height = height-(volumetric_flow_rate/width*length)*time;
end
