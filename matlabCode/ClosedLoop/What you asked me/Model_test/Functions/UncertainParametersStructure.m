function ThetaUncertain = UncertainParametersStructure(test)
% test == 1 => good parameters
%      ~= 1 =>  mismatch

ThetaUncertain = struct();

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WO reactor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5 reactions
%
% A +  B -> C      (1)
% C +  B -> P + E  (2)
% P +  C -> G      (3)
% A + 2B -> P + E  (4)
% A +  B + P -> G      (5)
%
if strcmp(test,'Plant') == 1
    %     % Plant 2
% ThetaUncertain.A1 = 1.6599*1e6;   % [1/h]
% ThetaUncertain.B1 = 6666.6666667; % [°K]
% ThetaUncertain.A2 = 7.2117*1e8;  % [1/h]
% ThetaUncertain.B2 = 8333.3333333; % [°K]
% ThetaUncertain.A3 = 2.6745*1e12;  % [1/h]
% ThetaUncertain.B3 = 11111.111111; % [°K]
% 
% ThetaUncertain.A4 = 0;  % [1/h]
% ThetaUncertain.B4 = 0; % [°K]
% ThetaUncertain.A5 = 0;  % [1/h]
% ThetaUncertain.B5 = 0; % [°K]
    
%     Plant
ThetaUncertain.A1 = 5.9755*1e9;   % [1/h]
ThetaUncertain.B1 = 6666.6666667; % [°K]
ThetaUncertain.A2 = 2.5962*1e12;  % [1/h]
ThetaUncertain.B2 = 8333.3333333; % [°K]
ThetaUncertain.A3 = 9.6283*1e15;  % [1/h]
ThetaUncertain.B3 = 11111.111111; % [°K]

ThetaUncertain.A4 = 0;  % [1/h]
ThetaUncertain.B4 = 1; % [°K]
ThetaUncertain.A5 = 0;  % [1/h]
ThetaUncertain.B5 = 1; % [°K]

else
    % 
%     % Plant 2 Mismatch
% ThetaUncertain.A1 = 0;   % [1/h]
% ThetaUncertain.B1 = 1; % [°K]
% ThetaUncertain.A2 = 0;  % [1/h]
% ThetaUncertain.B2 = 1; % [°K]
% ThetaUncertain.A3 = 0;  % [1/h]
% ThetaUncertain.B3 = 1; % [°K]
% 
% ThetaUncertain.A4 = 5.0728*1e13;  % [1/h] 
% ThetaUncertain.B4 = 8077.6; % [°K]
% ThetaUncertain.A5 = 1.4155*1e16;  % [1/h]
% ThetaUncertain.B5 = 12438.5; % [°K]
    
% %         % Model 3  Parametric mismatch
% ThetaUncertain.A1 = 6*1e9;   % [1/h]
% ThetaUncertain.B1 = 7000; % [°K]
% ThetaUncertain.A2 = 2.7*1e12;  % [1/h]
% ThetaUncertain.B2 = 8000; % [°K]
% ThetaUncertain.A3 = 9.5*1e15;  % [1/h]
% ThetaUncertain.B3 = 11000; % [°K]
% 
% ThetaUncertain.A4 = 0;  % [1/h]
% ThetaUncertain.B4 = 0; % [°K]
% ThetaUncertain.A5 = 0;  % [1/h]
% ThetaUncertain.B5 = 0; % [°K]

%         % Model 3  Parametric mismatch
ThetaUncertain.A1 = 0;   % [1/h]
ThetaUncertain.B1 = 1; % [°K]
ThetaUncertain.A2 = 0;  % [1/h]
ThetaUncertain.B2 = 1; % [°K]
ThetaUncertain.A3 = 0;  % [1/h]
ThetaUncertain.B3 = 1; % [°K]

% Model 1
% ThetaUncertain.A4 = 5.0728*1e10;  % [1/h] 
% ThetaUncertain.B4 = 7050; % [°K]
% ThetaUncertain.A5 = 8.4155*1e11;  % [1/h]
% ThetaUncertain.B5 = 8500; % [°K]

% Model 2
ThetaUncertain.A4 = 5.0728*1e10;  % [1/h] 
ThetaUncertain.B4 = 7250; % [°K]
ThetaUncertain.A5 = 8.4155*1e11;  % [1/h]
ThetaUncertain.B5 = 8320; % [°K]

end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Distillation column
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(test,'Plant') == 1
    % Plant
ThetaUncertain.eff = 0.1;
else
   % Model
% Model 1
ThetaUncertain.eff = 0.01;
% Model 2
% ThetaUncertain.eff = 0.5;

end

end