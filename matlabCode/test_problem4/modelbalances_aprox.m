function [residuals] = modelbalances_aprox(x,u,E_1, E_2)
% u = [Fa, Fb, Tr]
% x = [xa, xb, xp, xe, xg]
% uncertain parameters: E_1 and E_2
% This is Marchetti's two-reaction approximation (STRUCTURAL MISMATCH)

%% uncertain parameters
k_10 = 2.189*1e8;
k_20 = 4.310*1e13;
%E_1 = 8077.6;
%E_2 = 12438.5;


%constants
Ma = 1;
Mb = 1;
Mp = 1; %Marchetti assumes these
Me = 2;
Mg = 3;

W = 2105; %kg

%inputs
Fa = u(1); %kg/s
Fb = u(2); %kg/s
Tr = u(3); %celsius

%states
xa = x(1);
xb = x(2);
xp = x(3);
xe = x(4);
xg = x(5);



%rate constants
k1 = k_10*exp(-E_1/(Tr + 273.15));
k2 = k_20*exp(-E_2/(Tr + 273.15));

r1 = k1*xa*xb^2;
r2 = k2*xa*xb*xp;

%material balance equations
residuals = [Fa - (Fa + Fb)*xa - W*r1 - W*r2,...
            Fb - (Fa +Fb)*xb - (Mb/Ma)*2*W*r1 - (Mb/Ma)*W*r2,...
            -(Fa + Fb)*xp + (Mp/Ma)*W*r1 - (Mp/Ma)*W*r2,...
            -(Fa +Fb)*xe + (Me/Ma)*W*r1,...
            -xg + (Mg/Me)*xe + (Mg/Mp)*xp];
             
        
        
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        