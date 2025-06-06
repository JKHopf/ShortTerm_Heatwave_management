% PredUrchinKelp_ImplicitCC_v6b_MngtScen_GSA.m
% Nov 2024
% 
% Jess Hopf
% jess.k.hopf@gmail.com 

% Toolboxes requried:
%   - Statistics and Machine Learning Toolbox
%   - Optimization Toolbox

% GLOBAL Sensitvity Analysis for PredUrchinKelp_ImplicitCC_v6b_MngtScen_v1.m


% Top-level stuff --------------------

    % clear all
    clear
    
    % add paths
    addpath('functions\', 'functions\cbrewer')
    addpath('Model outputs\')

    % set random number
    % rng(1)

    tic
        
% ------- Model parameters -------------

% runtimes
% runs winter, spring, summer, autmum
% disturbance happens in yr 20-23
T1 = 30*4; %                pre for fished predator
T2 = 35*4; % 25*4; %        for kelp-urchin runs 

% number of replicates
RR = 10000;%    2;%  

% pre-assign variables
GSA_mat = NaN(10000,26);

% Which management scenario to run over?
mngt_scen = 'none'; 
          

parfor sa = 1:10^4

    % load in parameter values
    tmax = 10^4;
    kelp = ParaKelp_Implicit_v6bGSA(tmax);
    urchin = ParaUrchin_Implicit_v6bGSA(tmax);
    % pred = ParaPred_Implicit_v6(tmax);
    
    
    % --- Disturbance -------
    dist = struct();
    
    % single acute disturbance 
        % length of disturbance (years)
        dist.lngth = 2; %
    
        % time in model when disturbance happen  (20 years pre run is added to this)
        % 1 & 2 = winter and spring 3 & 4 = summer & autum
            dist.yrs = NaN; % no disturbance scenario
    
       % How do vital rates change during the disturbance (heatwave) 
        % kelp recruitment
        dist.RK = kelp.RK/7; %  0; % 
        % change in kelp biomass
        dist.lambda = kelp.lambda .* repmat([1 1 0.5 0.5],1,tmax/4); %  * 0.5; %  
        % urchin grazing rates
        dist.hij = repmat(kelp.bhij.*reshape([1.15 1.05 1.2 1.3],1,1,4),1,1,1,tmax/4); % kelp.hij * 1.2; % * 1.4; % 
      
    % get vector values
    mngt = ParaMngt_Implicit_v6_v1(mngt_scen);
    
    % degree of mngt action
    % if contains(mngt_scen,'fish')
    %    pred.fishF = mngt.degreeF(1); end
    if contains(mngt_scen,'cull') 
       urchin.culln = mngt.degreeC(1); end
    if contains(mngt_scen,'rest') 
       kelp.restn = mngt.degreeR(1);  end   
    
    % pred.fish = mngt.fish;
    urchin.culling = mngt.culling;
    kelp.restore = mngt.restore;


    % Intial conditions
        % KELP 🌿
        % 1.17*10^5 kg.ha = time-average kelp biomass at Anacapa East MPA (in kelp forest)
        % start with high drift so that it give the system a chance to be in a kelp state
        kt0 = [1.17*10^5,1.17*10^5,1.17*10^6]; %  [0,1,0]; %  [0,0,0]; %  
        
        % URCHINS 🟣
        % 6.7x10^3 kg.ha = time-average urchin biomass at Anacapa East MPA
        % 2.64x10^4 kg.ha = threshold urchin biomass value (above which urchins dominate in data)
        % start with low numbers to give system a chance to be in a kelp state
        ut0 = [0,0,0];%  [urchin.RU,urchin.RU,0]; %


    % run pred model over time
    % note T1 is to get stable state and T2 is the same as the
    % urchin-kelp run
    % [nt,nb] = run_Predator_Implicit_v6(pred, T1+T2, RR, ones(pred.meshno,1), dist);

    % biomass predating on urchins 
    % (turn off for GSA and set as a para to test)
    % PBE(:,h,i,j,:) = sum(nb(pred.Lgraze_ind:end,(end-T2+1):end,:),1);
    PBEsa = 1 + (1000-1).*rand(1);
    PBE = repmat(PBEsa,T2, RR);

    % urchin-kelp section
    % note that input predator biomass needs to be size T2 x RR)
    [kt2,ut2,~,RK_noise] = run_UrchinKelp_Implicit_v6b(kelp, urchin, T2, RR, ...
                     kt0, ut0, PBE, dist);

    % calculate kelp average (mean over last 1 year of run)
    kelp_avg = mean(kt2(2,(end-4):end,:));

    % proportion of sims with kelp persistence at sample time
    kelp_pers_t = nnz(kt2(2,88+8*4+1,:)>0);

    % save info
    GSA_mat(sa,:) = [kelp_pers_t,PBEsa,...
        kelp.RK,kelp.RKstdv,kelp.mu,kelp.ddD,kelp.muvar,kelp.g,kelp.rS,...
        kelp.c,kelp.rD,kelp.d,kelp.ajae,kelp.adc,kelp.h,...
        urchin.RU,urchin.RUstdv,urchin.gJ,urchin.MJ,...
        urchin.MH,urchin.ME,urchin.PE,urchin.PH, ...
        urchin.PLD,urchin.w2,urchin.kmin];

end  

toc

% 
% save(['Implicitv6a_MngtScen_GSA_',mngt_scen,'_',datestr(now, 'yyyymmdd'),'.mat'], '-v7.3')
% writematrix(GSA_mat, ['Implicitv6b_MngtScen_GSA_',mngt_scen,'_',datestr(now, 'yyyymmdd'),'.csv'])
% 
% histogram(GSA_mat(:,1)./RR)