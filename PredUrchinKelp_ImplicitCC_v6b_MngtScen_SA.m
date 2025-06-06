% PredUrchinKelp_ImplicitCC_v6b_MngtScen_SA.m
% Nov 2024
% 
% Jess Hopf
% jess.k.hopf@gmail.com 

% Toolboxes requried:
%   - Statistics and Machine Learning Toolbox
%   - Optimization Toolbox

% Sensitvity analysis for PredUrchinKelp_ImplicitCC_v6b_MngtScen_v1.m


% Top-level stuff --------------------

    % clear all
    clear
    
    % add paths
    addpath('functions\', 'functions\cbrewer')
    addpath('Model outputs\')

    % set random number
    rng(1)

    tic
        
% ------- Model parameters -------------

% runtimes
% runs winter, spring, summer, autmum
% disturbance happens in yr 20-23
T1 = 30*4; %                pre for fished predator
T2 = 35*4; % 25*4; %        for kelp-urchin runs 

% number of replicates
RR = 10000;% 2;%   


% load in parameter values
tmax = 1000;
kelp = ParaKelp_Implicit_v6b(tmax);
urchin = ParaUrchin_Implicit_v6(tmax);
pred = ParaPred_Implicit_v6(tmax);


% --- Disturbance -------

% single acute disturbance 
    % length of disturbance (years)
    dist.lngth = 2; %

    % time in model when disturbance happen  (20 years pre run is added to this)
    % 1 & 2 = winter and spring 3 & 4 = summer & autum
        % dist.yrs = (20*4) + repmat(1:4,1,dist.lngth) + repelem(((1:dist.lngth)-1)*4,4); % disturbance is on
        dist.yrs = NaN; % no disturbance scenario

   % How do vital rates change during the disturbance (heatwave) 
    % kelp recruitment
    dist.RK = kelp.RK/7; %  0; % 
    % change in kelp biomass
    dist.lambda = kelp.lambda .* repmat([1 1 0.5 0.5],1,tmax/4); %  * 0.5; %  
    % urchin grazing rates
    dist.hij = repmat(kelp.bhij.*reshape([1.15 1.05 1.2 1.3],1,1,4),1,1,1,tmax/4); % kelp.hij * 1.2; % * 1.4; % 


% ---------- Run models -----------------

% --- Management scenarios ------

% Which management scenario to run over?
mngt_scen = 'none'; %  'culling'; % 'fishing'; %  'restoration'; %   
    
% get vector values
mngt = ParaMngt_Implicit_v6_v1(mngt_scen);

% degree of mngt action
if contains(mngt_scen,'fish')
   pred.fishF = mngt.degreeF(1); end
if contains(mngt_scen,'cull') 
   urchin.culln = mngt.degreeC(1); end
if contains(mngt_scen,'rest') 
   kelp.restn = mngt.degreeR(1);  end   

pred.fish = mngt.fish;
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


% pre-assign variables
deglngth = structfun(@numel,mngt);  
deglngth = deglngth(end);

kts = NaN(3, T2+1, length(mngt.time), length(mngt.length), deglngth,RR);
uts = kts;
kelp_avg = NaN(length(mngt.time), length(mngt.length), deglngth, RR);
PBE = NaN(T2, length(mngt.time), length(mngt.length), deglngth, RR);  

% Run over range of values for SA parameter 
paravalSA = 0.0065.*[0.5,1,2];    % urchin.PH
% paravalSA = [0.25,0.5688,1];       % kelp.rS
% paravalSA = [0.75,0.9,0.95];      % kelp.c
% paravalSA = 6.825.*[0.5,1,2];     % kelp.g
% paravalSA = [255,355,455];      % PBE, through sheephead recruits

for h = 1:length(paravalSA)
    
    % chose which parameter
    urchin.PH = paravalSA(h);
    % kelp.rS = paravalSA(h);
    % kelp.c = paravalSA(h);
    % kelp.g = paravalSA(h);
    % pred.Rmean = paravalSA(h);


    % Run over timing of mngt action
    for i = 1:length(mngt.time)
    
        if contains(mngt_scen,'fish')
            pred.fishtime = mngt.time(i)+T1;    end
        if contains(mngt_scen,'cull') 
            urchin.culltime = mngt.time(i);     end  
        if contains(mngt_scen,'rest') 
            kelp.resttime = mngt.time(i);       end
    
    
        % Run over length of mngt action
        for j = 1:length(mngt.length)
        
            if contains(mngt_scen,'fish')
                pred.fishlgth = mngt.length(j);     end
            if contains(mngt_scen,'cull') 
                urchin.culllgth = mngt.length(j);   end  
            if contains(mngt_scen,'rest') 
                kelp.restlgth = mngt.length(j);     end
                    
     
            % run pred model over time
            % note T1 is to get stable state and T2 is the same as the
            % urchin-kelp run
            [nt,nb] = run_Predator_Implicit_v6(pred, T1+T2, RR, ones(pred.meshno,1), dist);

            % biomass predating on urchins
            PBE(:,h,i,j,:) = sum(nb(pred.Lgraze_ind:end,(end-T2+1):end,:),1);
    
            % urchin-kelp section
            % note that input predator biomass needs to be size T2 x RR)
            [kt2,ut2,~,RK_noise] = run_UrchinKelp_Implicit_v6b(kelp, urchin, T2, RR, ...
                             kt0, ut0, squeeze(PBE(:,h,i,j,:)), dist);
        
            % outputs
            kts(:,:,h,i,j,:) = kt2;
            uts(:,:,h,i,j,:) = ut2;
        
            % calculate kelp average (mean over last 1 year of run)
            kelp_avg(h,i,j,:) = mean(kt2(2,(end-4):end,:));

            
    
        end    
    end
end

% proportion of sims with kelp persistence over time
kelp_pers_t = squeeze(sum(kts(2,88+8*4+1,:,:,:,:)>0,6));

reshape(kelp_pers_t,[],1)

% figure(98324)
% hold on



     
% 
% save(['Implicitv6a_MngtScen_v1_',mngt_scen,'_',datestr(now, 'yyyymmdd'),'.mat'],... % '_degreesLRG.mat'],... % 
%       "kelp_pers_t", "T1", "T2", "kelp", "urchin", "pred", 'dist', 'mngt', 'mngt_scen', '-v7.3')

% for baseline no disturbance
    % baseline_persit_overtime = kelp_pers_t(:,:,:,:,1);
    % save(['Implicitv6a_MngtScen_v1_nodist_',datestr(now, 'yyyymmdd'),'.mat'],... % '.mat'],... %
    %  "baseline_persit_overtime", '-v7.3')

% for disturbance no mngt action
    % % baseline
    % baseline_persit_overtime = kelp_pers_t(:,:,:,:,1);
    % 
    % % mean values before disturbance
    % persisting_runs = find(squeeze(prod(kts(2,1:(dist.yrs(1)-1),:,:,:,:)>0,2))>0)';
    % kelpJ_avg_pre = mean(kt2(1,(dist.yrs(1)-40):(dist.yrs(1)-1),persisting_runs),[2,3]);
    % kelpA_avg_pre = mean(kt2(2,(dist.yrs(1)-40):(dist.yrs(1)-1),persisting_runs),[2,3]);
    % urchinA_avg_pre = mean(sum(ut2(2:3,(dist.yrs(1)-40):(dist.yrs(1)-1),persisting_runs),1),[2,3]);
    % save(['Implicitv6a_MngtScen_v1_',mngt_scen,'_',datestr(now, 'yyyymmdd'),'.mat'],... % '.mat'],... %
    %       "kelpJ_avg_pre","kelpA_avg_pre","urchinA_avg_pre", "baseline_persit_overtime", '-v7.3')
    % 


