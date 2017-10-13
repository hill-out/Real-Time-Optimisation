% solves the steady state outlets of an elementary reaction CSTR

function [cSol] = CSTR(reactOrder, kVal, c0, volFlow, V)

%check function variables (unfinished)
kSize = size(kVal);
rSize = size(reactOrder);
cSize = size(c0);
vSize = size(volFlow);

if kSize(1) ~= rSize(1) || rSize(2) ~= cSize(2) || cSize(2) ~= vSize(2)
    error('Function input dimensional mismatch')
end

%calc variables
v0 = sum(volFlow);
nReact = size(reactOrder);

%solve steady state
options = optimset('Display','off');
cSol = fsolve(@massBal,c0,options);

    function [F] = massBal(c)
        %evaluates the species balance
        r = rateCalc(c);
        F = volFlow.*c0/V+r-v0*c/V;
    end

    function [r] = rateCalc(c)
        %evaluates the rate of reaction
        c = repmat(c,nReact(1),1);
        
        elemOrder = reactOrder;
        elemOrder(elemOrder>0)=0;
        baseM = c.^(-elemOrder);
        rReact = kVal.*prod(baseM,2);
        r = rReact'*reactOrder;
    end


end