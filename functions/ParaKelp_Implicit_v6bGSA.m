function kelp = ParaKelp_Implicit_v6bGSA(tmax)

% Sets para values for kelp portion of the model.
% Relevant to PredUrchinKelp_ImplicitCC_v6.m
% r = a + (b-a).*rand(1) (a,b)

% Recruitment
    % yearly incomming settlers per kg adult standing kelp
    % (zoospore production and successful fertilization and settlement of spore)
    % mean value 4*10^4
    kelp.RK = 10^1 + (10^8-10^1).*rand(1); 

    % temporal stdev (noise)
    % normalised 0.389
    kelp.RKstdv = 0 + (0.9-0).*rand(1); 

    % **strength of denisity dependence (shading by local adults)
    % smaller = weaker effect = higher numbers
    % Set so that kg.60m^2 is within the range seen in the Channel Islands
    % (~ 700-800kg.60m2) 0.00001
    kelp.mu = 0.5*10^-5 + (2*10^-5 - 0.5*10^-5).*rand(1); %  0.1*10^-5 + (10*10^-5 - 0.1*10^-5).*rand(1); % 

        % the relative per capita effect on juvenile survival by adults and juveniles
        % ddD = 0 for ricker (inter-cohort DD, adults->)
        % ddD = 1 for beverton-holt (intra-cohort DD, juvs <->)
        % we'll assume that one adult has the same effect as one juvenile on the survival for giant kelp
        kelp.ddD = 0 + (0.5-0).*rand(1);

        % spatial variance in kelp densitites 
        % (for use in scale transition to landscape)
        % value from PISCO_kelp-urchin-sheephead_data.Rmd 189090;
        kelp.muvar = 0 + (10^10-0).*rand(1);

    % Reproduction timing
    % vector that decitates if, and how much, reproduction occurs in that time step
    % [winter, spring, summer, autumn]
    % summer & autum peaks
    kelp.RTk = repmat([0.1 0.1 0.4 0.4],1,tmax/4);
    % kelp.RTk = repmat([0.25 0.25 0.25 0.25],1,tmax/4); 
    
% Growth/Maturation    
    % growth rate season to season 6.825
    kelp.g = 1 + (20-1).*rand(1); 

% Mortality/survival
    % change in standing kelp biomass 
    % can including seasonal variation here [winter spring summer autum]
    % kelp.lambda = repmat([1 1 0.8 0.9],1,tmax/4);
    kelp.lambda = repmat([1 1 1 1],1,tmax/4);  % 1; 
 
    % standing (juvenile + adult) kelp retension 0.5688;
    kelp.rS = rand(1);     

    % drift production [0,1]
    % (standing kelp mortality)
    % low = more kelp, high = more drift
    kelp.c = rand(1);  

    % ++ drift retention [0,1]
    % low = less retained, high = more retained
    kelp.rD = rand(1);

    % ++ decomposition [0,1]
    % low = slower decomp, high = fast decomp
    kelp.d = rand(1); 

% Grazing
    % attack rates (of urchin stage j on kelp stage i)
    % (p.c. grazing mortality in absence of conspecifics)
    kelp.ajae = rand(1);  %0.5
    kelp.adc = rand(1);  %1

    kelp.aij = [ 0  kelp.ajae;
                 0  kelp.ajae;
                 kelp.adc  0 ];

    % max prey consumed (1/handling time) 
    % (of urchin stage j on kelp stage i)    
    kelp.h = 1 + (10-1).*rand(1); % 2.985
    kelp.bhij = [ 0  kelp.h;  
                  0  kelp.h;
                  kelp.h  0 ];

    % not seasonal
    % kelp.hij = repmat(hij.*reshape([1 1 1 1],1,1,4),1,1,1,tmax/4);

    % including seasonal variation
    % [winter spring summer autum]
    kelp.hij = repmat(kelp.bhij.*reshape([1 0.9 1.15 1.2],1,1,4),1,1,1,tmax/4);
 

end