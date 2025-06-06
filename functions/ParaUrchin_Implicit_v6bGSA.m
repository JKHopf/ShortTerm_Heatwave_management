function urchin = ParaUrchin_Implicit_v6bGSA(tmax)

% Sets para values for urchin portion of the model.
% Relevant to PredUrchinKelp_ImplicitCC_v6.m


% Recruitment
    % **constant incoming settlers (mean value; inc fecundity & dispersal)
    % This value makes all the difference to what state the model is in! 
    % - set to 0 to turn off urchins in model
    % - threshold urchin biomass value = 2.64x10^4 kg.ha, above which urchins dominate) get this in non-MPA scenario with RU = 3.605*10^5  
    % - 1*10^6 gives urchin values seen at Anacapa East, with sheephead present (6.7x10^3 kg.ha), but urchins dominate when sheephead are
    %   under fishing pressure
    % Set so that kg.60m^2 is within the range seen in the Channel Islands
    % (~60 kg.60m2 outside MPAs, ~24 kg.60m2 in MPAs) 3*10^5;
    urchin.RU = 10^1 + (10^8-10^1).*rand(1);

    % temporal stdv of recruitment (noise) 0.621
    urchin.RUstdv = rand(1);

    % Reproduction timing
    % vector that decitates if reproduction occurs in that time step
    % urchins reproduce 90% March-July (March-May = spring)
    % urchin.RTu = repmat([0.25 0.25 0.25 0.25],1,tmax/4);
    urchin.RTu = repmat([0.05 0.54 0.36 0.05],1,tmax/4);

% Growth    
    % proportion maturing from juvenile to adults
    % 4 seasons in a year * 2 = 2 yrs maturation 1/(4*2); 
    urchin.gJ = rand(1);

% Mortlaity rates (instantaneous)
    % natural mortality
        % juveniles 0.1;
        LA = log(0.001); LB = log(2);
        urchin.MJ = exp(LA + (LB-LA) * rand(1)); 
            % strength of DD 
            % smaller = stronger DD effects (survival increases slower with
            % inc adult biomass) 0.00001;
            LA = log(0.0000001); LB = log(0.01);
            urchin.alpha = exp(LA + (LB-LA) * rand(1));
            
            % variance in adult urchin densitites 
            % (for use in scale transition to landscape)
            % value from PISCO_kelp-urchin-sheephead_data.Rmd 16113; % 
            urchin.alphavar = 10^1 + (10^7-10^1).*rand(1);

        % hiding adults 0.1;
        LA = log(0.0001); LB = log(2);
        urchin.MH = exp(LA + (LB-LA) * rand(1));
        % exposed adults 0.1;
        urchin.ME = exp(LA + (LB-LA) * rand(1));

    % predation rate on urchins by predators (lower = less mortality)
        % exposed 0.013;
        LA = log(0.00001); LB = log(2);
        urchin.PE = exp(LA + (LB-LA) * rand(1));
        % hiding urchin.PE*0.5; 0.0065
        urchin.PH = exp(LA + (LB-LA) * rand(1));


    % fishing rate
        urchin.F = 0;   

    % discount rate for juv survival
    % PLD in days 65;
    urchin.PLD = 10 + (300-10).*rand(1);
    urchin.tau = 1-urchin.PLD/91;


% **behavioural switching function (hiding <> exposed)
    % inflection point (drift density at which urchin hiding:exposed = 1:1)
    urchin.w1 = 1; % 
    % slope around inflection point 
    % (larger values = steeper inflection = sharper transition)
    % (smaller values = flatter line = urchins out all the time) 0.5; % 
    urchin.w2 = rand(1);

% **step function (turns predation/fishing of exposed on/off) 
% min density of all standing kelp (juvs + adults)
% 0 = events never off (predation on)
% large number = events always off (predation off) 1170;% 0;% 10^10; %  
    urchin.kmin = 1 + (10^10-1).*rand(1);

end