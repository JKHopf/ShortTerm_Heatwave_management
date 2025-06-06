function mngt = ParaMngt_Implicit_v6_v1(scenario)

% set default scenarios off
mngt.fish = "N"; 
mngt.culling = "N"; 
mngt.restore= "N"; 
mngt.time = NaN;
mngt.length = NaN; 
mngt.degree = NaN;

% set time and length for any scenario (if one occurs)

if contains(scenario,{'fish', 'cull', 'rest'}) 

    % timing of action start
    % relative to start of disturbance, in terms of seasons (timesteps)
    % before = -x, during = 0, after = x-1 (yr)
    mngt.time = -4; % [-4,0,4]; %   (-8:2:20); %     

    % length of mngt action (timesteps)
    mngt.length = 12; % 20; % [0,4,12,20]; %  0:4:20; %   

end


% if the scenario contains any of the strings then the action will be
% implimented

if contains(scenario,'fish') 
    % Temporarily reducing fishing pressure (or temp MPA)
    mngt.fish = "Y"; 

    % temporary fishing pressure (baseline for SH = 0.05 per season)
    mngt.degreeF = 0; % [0,0.05]; % 0:0.005:0.05; %  

end

if contains(scenario,'cull') 
    % Urchin Removal
    % = reducing urchin biomass by removing a set amount 
    mngt.culling = "Y";  

    % biomass removed in a season (kg.ha)
    % avg pre-disturbance urchin biomass, get value from:
        % model was run with disturbance & no mngt action, to get biomass
        % time-step before disturbance
        % load("Model outputs\MngtScenarios\Implicitv6b_MngtScen_v1_none_D0.01_20250305.mat")
        % urchinA_avg_pre = 1.7*10^4
    mngt.degreeC = 1*1.7*10^4; %   [0,0.01,0.05,0.1,0.25,0.5,0.75,1].*1.7*10^4; %   

end
    
if contains(scenario,'rest') 
    % Kelp restoration  
    % = reseeding juveniles into the population
    mngt.restore = "Y"; 

    % biomass of recruits added per season (kg.ha)
    % avg pre-disturbance juv + adult standing kelp biomass (fished state)
    % get value from:
        % model was run with disturbance & no mngt action, to get biomass
        % time-step before disturbance
        % load("Model outputs\MngtScenarios\Implicitv6b_MngtScen_v1_none_D0.01_20250305.mat")
        % calculate kelpJA_avg_pre = kelpJ_avg_pre + kelpA_avg_pre;
        % = 1.43*10^5 kg for D=0.5;  1.63*10^5 for D=0.01
    mngt.degreeR = 1*1.43*10^5; %  [0,0.01,0.05,0.1,0.25,0.5,0.75,1].*1.43.*10^5; %

end


end 