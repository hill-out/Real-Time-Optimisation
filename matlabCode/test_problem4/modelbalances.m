function [residuals] = modelbalances(x,u, theta)
% u = [Fa,Fb, Tr]
% x = [xa, xb, xc, xp, xe, xg]
% uncertain parameters: 
%theta(1) = kinetic constantl (plant value = 1.6599*1e6)
%theta(2) = exponent in a rate law (plant val = 1)
% This is the structurally perfect model

%constants
Fa = u(1); %kg/s
Ma = 1;
Mb = 1;
Mp = 1; %Marchetti assumes these
Mc = 2; %from stoichiometry
Me = 2;
Mg = 3;

W = 2105; %kg

%inputs
Fb = u(2); %kg/s
Tr = u(3); %celsius

%states
xa = x(1);
xb = x(2);
xc = x(3);
xp = x(4);
xe = x(5);
xg = x(6);



%rate constants
k1 = theta(1)*exp(-6666.7/(Tr + 273.15));
k2 = 7.2117*1e8*exp(-8333.3/(Tr + 273.15));
k3 = 2.6745*1e12*exp(-11111/(Tr + 273.15));

r1 = k1*xa*xb;
r2 = k2*xb*xc^(theta(2));
r3 = k3*xc*xp;

%material balance equations
residuals = [Fa - (Fa + Fb)*xa - W*r1,...
            Fb - (Fa +Fb)*xb - (Mb/Ma)*W*r1 - W*r2,...
            -(Fa + Fb)*xc + (Mc/Ma)*W*r1 - (Mc/Mb)*W*r2 - W*r3,...
            -(Fa + Fb)*xp + (Mp/Mb)*W*r2 - (Mp/Mc)*W*r3,...
            -(Fa +Fb)*xg + (Mg/Mc)*W*r3,...
            -xe + (Me/Mp)*xp + (Me/Mg)*xg];
             
        
        
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        