function [residuals] = plantbalances(x,Xa,Fb,Fa)
% y = [xa (scaled), Fb]
% x = [Tr, xb, xc, xp, xe, xg]
%e_c = control error on xa (plant val = 0)

%constants       
Ma = 1;
Mb = 1;
Mp = 1; %Marchetti assumes these
Mc = 2; %from stoichiometry
Me = 2;
Mg = 3;

W = 2105; %kg

%inputs
xa = Xa/100; %the real x_a 



%states
Tr = x(1);
xb = x(2);
xc = x(3);
xp = x(4);
xe = x(5);
xg = x(6);



%rate constants
k1 = 1.6599*1e6*exp(-6666.7/(Tr + 273.15));
k2 = 7.2117*1e8*exp(-8333.3/(Tr + 273.15));
k3 = 2.6745*1e12*exp(-11111/(Tr + 273.15));

r1 = k1*xa*xb;
r2 = k2*xb*xc;
r3 = k3*xc*xp;


%material balance equations
residuals = [Fa - (Fa + Fb)*xa - W*r1,...
            Fb - (Fa +Fb)*xb - (Mb/Ma)*W*r1 - W*r2,...
            -(Fa + Fb)*xc + (Mc/Ma)*W*r1 - (Mc/Mb)*W*r2 - W*r3,...
            -(Fa + Fb)*xp + (Mp/Mb)*W*r2 - (Mp/Mc)*W*r3,...
            -(Fa +Fb)*xg + (Mg/Mc)*W*r3,...
            -xe + (Me/Mp)*xp + (Me/Mg)*xg];
             
      
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        