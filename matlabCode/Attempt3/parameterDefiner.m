function [out] = parameterDefiner(para)
% -------------------------------------------------------------------------
% This function takes the input and searches the database style set-up to
% get the required parameters as out
%
% para      - [string]      - required parameters:
%                               'm' = model
%                               'p' = plant
%                               'c' = contraints
% 
% out       - [cell]        - output of the parameters required
%
% -------------------------------------------------------------------------

% check inputs
if ~isa(para,'char')
    error('Expected "para" to be type ''char'' but got type ''%s''',class(para))
end

% initialise out
out = cell(1,numel(para));

% run for each letter in para
for i = 1:numel(para)
    if para(i) == 'm' % para(i) is for model
        
        % model specific parameters
        m.k = [0.75, 1.5];              % k value for model [L/(mol min)]
        m.c_in = [2; 1.5; 0; 0];        % conc_in for model [mol/L]
        
        % parameters for both model and plant
        m.s = [-1,  0;
               -1, -2;
                1,  0;
                0,  1];    % stoichiometry
        m.V = 500;         % volume of reactor [L]
        m.dH = [3.5, 1.5]; % heat of reaction
        m.o = -m.s;
        m.o(m.o<0) = 0; % assume elementary
        
        % constraints
        m.Qmax = 110;
        m.Dmax = 0.1;
        % add to out
        out{i} = m;
        
    elseif para(i) == 'p' % para(i) is for plant
        
        % plant specific parameters
        p.k = [1.4, 0.4];              % k value for model [L/(mol min)]
        p.c_in = [2.5; 1.5; 0; 0];        % conc_in for model [mol/L]
        
        % parameters for both model and plant
        p.s = [-1,  0;
               -1, -2;
                1,  0;
                0,  1];    % stoichiometry
        p.V = 500;         % volume of reactor [L]
        p.dH = [3.5, 1.5]; % heat of reaction
        p.o = -p.s;
        p.o(p.o<0) = 0; % assume elementary
        
        % constraints
        p.Qmax = 110;
        p.Dmax = 0.1;
        
        % add to out
        out{i} = p;
        
    end
end
        
        