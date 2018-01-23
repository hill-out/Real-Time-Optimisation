function [out] = CSTRsymb(in)
% % Symbolic code to evaluate fixed matrices for NE and SOC - Also used for MU

syms k1 k2 u1 u2 V cain cbin w dH1 dH2 Qmax Dmax
u = [u1; u2];
theta = [k1;k2];
 
a = 2*k1*k2*V^2;
b = V*(2*k2+k1)*(u1+u2);
c = (u1+u2)^2 + (u1*cain-u2*cbin)*V*k1 ;
d = - u2*cbin*(u1+u2);
cb = 1/6/a*(36*c*b*a-108*d*a^2-8*b^3+12*3^(1/2)*(4*c^3*a-c^2*b^2-18*c*b*a*d+27*d^2*a^2+4*d*b^3)^(1/2)*a)^(1/3)-2/3*(3*c*a-b^2)/a/(36*c*b*a-108*d*a^2-8*b^3+12*3^(1/2)*(4*c^3*a-c^2*b^2-18*c*b*a*d+27*d^2*a^2+4*d*b^3)^(1/2)*a)^(1/3)-1/3*b/a;
ca = u1*cain / (V*k1*cb+u1+u2);
cc  = k1*ca*cb*V/(u1+u2);
cd = k2*cb^2*V/(u1+u2);
 
Phi = -cc^2*(u1+u2)^2/(u1*cain) + w*(u1^2 + u2^2);
Phi_u = jacobian(Phi,u);
Phi_uu = jacobian(Phi_u,u);
Phi_utheta = jacobian(Phi_u,theta);

H = [ca; cb; cc; cd];
H_u = jacobian(H,u);
H_theta = jacobian(H,theta);

Q = V*k1*ca*cb*dH1 + V*k2*cb^2*dH2;
G1 = Q/Qmax - 1;
G1_u = jacobian(G1,u);
G1_uu = jacobian(G1_u,u);
G1_utheta = jacobian(G1_u,theta);

D = cd/(ca+cb+cc+cd);
G2 = D/Dmax - 1;
G2_u = jacobian(G2,u);
G2_uu = jacobian(G2_u,u);
G2_utheta = jacobian(G2_u,theta); 

symb = [k1 k2 u1 u2 V cain cbin w dH1 dH2 Qmax Dmax];
out = {subs(Phi, symb, in);
    subs(Phi_u, symb, in);
    subs(Phi_uu, symb, in);
    subs(Phi_utheta, symb, in);
    subs(H, symb, in);
    subs(H_u, symb, in);
    subs(H_theta, symb, in);
    subs(G1, symb, in);
    subs(G1_u, symb, in);
    subs(G1_uu, symb, in);
    subs(G1_utheta, symb, in);
    subs(G2, symb, in);
    subs(G2_u, symb, in);
    subs(G2_uu, symb, in);
    subs(G2_utheta, symb, in)};

end



