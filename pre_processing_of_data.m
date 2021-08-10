%% pre_processing_of_data script
D_t                = D_2010(Time_Start+1:Time_Length+Time_Start);
W_on_t             = W_e_on_2010(Time_Start+1:Time_Length+Time_Start);
W_of_t             = W_e_of_2010(Time_Start+1:Time_Length+Time_Start); 
S_t                = S_2010(Time_Start+1:Time_Length+Time_Start); %(LV)
%------------------------------------------------------------------------------------------------
% Reserve requirements to meet 3.5 std of demand-wind-solar forecast errors
lambda_up          = 1.5;                % Upwards spinning reserve
lambda_dn          = 1.5;                % Downwards spinning reserve
lambda_up1         = 3.5 - lambda_up;    % Standing reserve
lambda_dn1         = 3.5 - lambda_dn;    % Standing reserve
%
W_on_t_Original    = W_on_t;
W_of_t_Original    = W_of_t;
%
S_t_large          = S_t*S_large;        % Large scale PV (> 5 MW) generation (LV)
S_t_small          = S_t*(1-S_large);    % Small scale PV (< 5 MW) generation (LV)
S_t_Original_large = S_t_large;
S_t_Original_small = S_t_small;
%
D_t_Original       = D_t;
W_t                = W_on_t + W_of_t;
W_t_Original       = W_on_t + W_of_t;
S_t                = S_t_large + S_t_small; %(LV)
S_t_Original       = S_t;        
Net_D_t_Original   = D_t - W_t_Original - S_t_Original; %(LV)
DEMAND             = Net_D_t_Original;
%------------------------------------------------------------------------------------------------
NG               = size(gen_data,1);   % no. of generators
NT               = size(DEMAND,1);     % number of time periods (hours)
%------------------------------------------------------------------------------------------------
index_ref        = length(D_t);
index_1          = 1;
index_2          = index_ref;
min_period       = index_1;
max_period       = index_2;
window           = max_period - min_period + 1;
period_t1        = 0;
period_t2        = 0;
max_array_size   = max_period - min_period + 1;
%------------------------------------------------------------------------------------------------
%% Minimum allowable demand - minimum synchronised thermal generation to provide intertia and reserve requirements
Tot_c_t=zeros(1,max_array_size); % All the necessary curtailment to maintain D_MIN_t (LV)
S_c_inert=zeros(1,max_array_size);
S_c_am_inert=zeros(1,max_array_size);
W_c_t_Original=zeros(1,max_array_size); 
D_MIN_t = K_WIND_INERTIA*W_t + MIN_LOAD; % to curtail wind and solar output for min load
W_c_t1=zeros(1,max_array_size);% (LV)
W_c_t2=zeros(1,max_array_size);% (LV)
W_on_c_inert=zeros(1,max_array_size);
W_on_c_am_inert=zeros(1,max_array_size);
W_of_c_inert=zeros(1,max_array_size);
W_of_c_am_inert=zeros(1,max_array_size);
S_c_t_small=zeros(1,max_array_size); % (LV)
S_c_inert_s=zeros(1,max_array_size);% (LV)
S_c_am_inert_s=zeros(1,max_array_size);% (LV)
% Previous curtailment (to maintain required levels of inertia + reserve)
W_on_c_feas=zeros(1,max_array_size);
W_on_c_am_feas=zeros(1,max_array_size);
W_of_c_feas=zeros(1,max_array_size);
W_of_c_am_feas=zeros(1,max_array_size);
S_c_feas=zeros(1,max_array_size);
S_c_am_feas=zeros(1,max_array_size);
S_c_feas_s=zeros(1,max_array_size);
S_c_am_feas_s=zeros(1,max_array_size);

for i = min_period:max_period
    if DEMAND(i) < D_MIN_t(i)
        Tot_c_t(i) = D_MIN_t(i) - DEMAND(i); % total amount to be curtailed
        S_t_large(i)= S_t_large(i) - min(S_t_large(i),Tot_c_t(i)); % Firstly, large scale solar is curtailed (LV)
        S_c_inert(i)=(S_t_Original_large(i)>S_t_large(i));% count how many times large solar is curtailed for inertia requirement (logical)
        S_c_am_inert(i) = (S_t_Original_large(i)>S_t_large(i)).*((S_t_Original_large(i)-S_t_large(i)));%how much large solar is curtailed for inertia requirement
        if Tot_c_t(i) > S_t_Original_large(i)
            W_c_t1(i) = Tot_c_t(i) - S_t_Original_large(i);
            W_on_t(i) = W_on_t(i) - min(W_c_t1(i),W_on_t(i)); % Secondly, onshore wind is curtailed and up to its available level
            W_on_c_inert(i)=(W_on_t_Original(i)>W_on_t(i)); % count how many times offshore wind is curtailed for inertia requirement
            W_on_c_am_inert(i) = (W_on_t_Original(i)>W_on_t(i)).*((W_on_t_Original(i)-W_on_t(i)));% how much onshore wind is curtailed for inertia requirement
            if Tot_c_t(i) > (W_on_t_Original(i) + S_t_Original_large(i))
                W_c_t2(i) = W_c_t1(i) - W_on_t_Original(i);
                W_of_t(i) = W_of_t(i) - min(W_c_t2(i),W_of_t(i)); % Thirdly, offshore wind is curtailed
                W_of_c_inert(i)=(W_of_t_Original(i)>W_of_t(i));% count how many times onshore wind is curtailed for inertia requirement
                W_of_c_am_inert(i) = (W_of_t_Original(i)>W_of_t(i)).*((W_of_t_Original(i)-W_of_t(i)));% how much offshore wind is curtailed for inertia requirement
                if Tot_c_t(i) > (W_t_Original(i)+ S_t_Original_large) % Part added for solar curtailment (LV)
                    S_c_t_small(i)= W_c_t2(i) - W_of_t_Original(i); % Fourthly, small scale solar is curtailed (LV)
                    S_t_small(i)= S_t_small(i) - min(S_t_small(i),S_c_t_small(i));
                    S_c_inert_s(i)=(S_t_Original_small(i)>S_t_small(i));% count how many times small solar is curtailed for inertia requirement
                    S_c_am_inert_s(i) = (S_t_Original_small(i)>S_t_small(i)).*((S_t_Original_small(i)-S_t_small(i)));% how much small solar is curtailed for inertia requirement
                end
            end
        end
    end
end
%
Tot_c_t     = Tot_c_t'; % (LV)
S_c_t_small = S_c_t_small'; % (LV)
S_t         = S_t_large + S_t_small; % (LV)
W_t         = W_on_t + W_of_t; % (LV)
W_c_t       = W_c_t1' - S_c_t_small; % (LV)
S_c_t_large = Tot_c_t - W_c_t - S_c_t_small; % (LV)
S_c_t       = S_c_t_large + S_c_t_small; % (LV)
DEMAND      = DEMAND + Tot_c_t; % (LV)
%
R_SP_UP_MAX = ones(NT,1)*P_MAX; %credible maximum generator at each time
R_SP_DN_MAX = zeros(NT,1); 
%
R_SP_UP_D_t = zeros(NT,1);
R_SP_DN_D_t = zeros(NT,1);
R_SP_UP_W_t = zeros(NT,1);
R_SP_DN_W_t = zeros(NT,1);
R_SP_UP_S_t = zeros(NT,1);
R_SP_DN_S_t = zeros(NT,1);
%
x = min_period:max_period;
R_SP_UP_D_t(x) = D_u*D_t(x); % Demand forecast uncertainty (std) is D_u*100% of forecast demand
R_SP_DN_D_t(x) = D_u*D_t(x); % Demand forecast uncertainty (std) is D_u*100% of forecast demand
R_SP_UP_W_t(x) = min(W_u*W_t(x),max(0,W_t(x) - W_c_t(x))); % Wind uncertainty is 10% of forecast wind or less if curtailed or zero if D_MIN_t>=D_t
R_SP_DN_W_t(x) = W_u*W_t(x);  % Wind uncertainty is W_u*100% of forecast wind
R_SP_DN_S_t(x) = S_u*S_t(x);  % Solar forecast uncertainty is S_u*100% of forecast solar
% R_SP_UP_S_t(x) = min(S_u*S_t(x),max(0,S_t(x) - S_c_t(x))); 
R_SP_DN_S_t(x) = S_u*S_t(x); % Solar forecast uncertainty is S_u*100% of forecast solar
%
R_SP_UP_D_W_S_t = ((R_SP_UP_D_t).^2 + (R_SP_UP_W_t).^2 + (R_SP_UP_S_t).^2).^0.5;
R_SP_DN_D_W_S_t = ((R_SP_DN_D_t).^2 + (R_SP_DN_W_t).^2 + (R_SP_DN_S_t).^2).^0.5;
%
RES_UP = R_SP_UP_MAX + lambda_up*R_SP_UP_D_W_S_t;
RES_UP(RES_UP<0) = 0;        
RES_DN = R_SP_DN_MAX + lambda_dn*R_SP_DN_D_W_S_t;
RES_DN(RES_DN<0) = 0;        
%
GMAX_0=GMAX;
GMIN_0=GMIN;
%----------------------------------------------------------------------------------------------------
%% check up the availability of data for certain cases:
if (DISPATCH_METHOD == 2 || DISPATCH_METHOD == 3) && (any(isnan(GNLC)) || any(isnan(GFC)) || any(isnan(GINC)))
    STR = ['To use linear cost model, you must provide data for NO LOAD COSTS,'...
        'FUEL COSTS and INCREMENTAL COSTS.'];   % there are no data for quick dispatch method,
    msgbox(STR,'DATA AVAILABILITY CHECK FAILURE!','warn');       % write a message
    return
elseif DISPATCH_METHOD == 1 && (any(isnan(COEF_A)) || any(isnan(COEF_B)) || any(isnan(COEF_C)))
    STR = ['To use quadratic cost model, you must provide data for the cost coefficients:,'...
        'COEFF_A (ï¿½), COEFF_B (ï¿½/MWh) and COEFF_C (ï¿½/MW^2h).'];   % there are no data for quick dispatch method,
    msgbox(STR,'DATA AVAILABILITY CHECK FAILURE!','warn');       % write a message
    return
end
if MIN_UP_DOWN_TIME_FLAG == 1 && (any(isnan(GMINUP)) || any(isnan(GMINDOWN)))
    STR = ['To use minimum up and down time constraints, you must provide data for GMINUP'...
        ' and GMINDOWN.'];
    msgbox(STR,'DATA AVAILABILITY CHECK FAILURE!','warn');
    return
end
if RAMP_UP_DOWN_FLAG == 1 && (any(isnan(GRAMPUP)) || any(isnan(GRAMPDOWN)))
    STR = ['To use rump constraints, you must provide data for GRAMPUP '...
        'and GRAMPDOWN.'];
    msgbox(STR,'DATA AVAILABILITY CHECK FAILURE!','warn');
    return
end
if COMPLETE_ENUMERATION_FLAG == 0 && (any(isnan(GNLC)) || any(isnan(GFC)) || any(isnan(GINC)))
    STR = ['To use priority list, you must provide data for NO LOAD COSTS,'...
        'FUEL COSTS and INCREMENTAL COSTS.'];
    msgbox(STR,'DATA AVAILABILITY CHECK FAILURE!','warn');
    return
end
if START_UP_COST_METHOD == 2 && (any(isnan(GSH)) || any(isnan(GCSTIME)) )
    STR = ['To use cold/hot start up cost method, you must provide data for GSH '...
        'and GCSTIME.'];
    msgbox(STR,'DATA AVAILABILITY CHECK FAILURE!','warn');
    return
elseif START_UP_COST_METHOD == 3 && (any(isnan(GSH)) || any(isnan(TAU)) )
    STR = ['To use exponential start up cost method, you must provide data for GSH '...
        'and TAU.'];
    msgbox(STR,'DATA AVAILABILITY CHECK FAILURE!','warn');
    return
end
%------------------------------------------------------------------------------------------------
%% enabling / disabling some of the constraints
if MIN_UP_DOWN_TIME_FLAG == 0               % minimum up and down time enabled/disabled
    GMINUP(:)   = 1;                        % if disabled, all min up and down times
    GMINDOWN(:) = 1;                        % are set to 1
end
%
if RAMP_UP_DOWN_FLAG == 0                   % ramping constraints enabled/disabled
    GRAMPUP(:)   = Inf;                     % if disabled, ramp rates are
    GRAMPDOWN(:) = Inf;                     % set to a large number
end
%
if COMPLETE_ENUMERATION_FLAG == 0           % use priority list...
    [GMAXlst,GMINlst,LIST_STATES,GEN_ORDER] = priority_list(GNLC,GFC,GMAX,GMIN,GINC,NG,ECO2_g,CoC_g,Y_MAX_CAPT_g,STK_g,VOM_g,CC_VOM_g,CC_TRANS_g,STACKING_OPTION);
elseif COMPLETE_ENUMERATION_FLAG == 1       % ...or complete enumeration (consisting of all possible combinations) or... 
    [GMAXlst,GMINlst,LIST_STATES,GEN_ORDER] = complete_enumeration(GNLC,GFC,GMAX,GMIN,GINC,NG,ECO2_g,CoC_g,Y_MAX_CAPT_g,STK_g,...
        VOM_g,CC_VOM_g,CC_TRANS_g,STACKING_OPTION);
else                                        % ...flexible enumeration
    [GMAXlst,GMINlst,LIST_STATES,GEN_ORDER] = dynamic_flex(GNLC,GFC,GMAX,GMIN,GINC,NG,ECO2_g,CoC_g,Y_MAX_CAPT_g,STK_g,...
        VOM_g,CC_VOM_g,CC_TRANS_g,STACKING_OPTION);
end
%
if RESERVE_FLAG == 0
    RES_UP = zeros(size(DEMAND));           % if reserve not required,
    RES_DN = zeros(size(DEMAND));           % set it to zero.
end
%
if START_UP_COST_METHOD == 3                % if start-up costs are exponential
    ALPHA = GSH;                            % then define ALPHA
    BETA  = GSC;                            % and BETA
else
    ALPHA = NaN*ones(NG,1);                 % otherwise, just define the names for variables
    BETA  = NaN*ones(NG,1);                 % since they will be passed to the functions
end
%------------------------------------------------------------------------------------------------
%% Determines the initial status (ON/OFF = 1/0) for each generator, based on the input data (GSTATINI).
% (GSTATINI contains the number of hours that a generator was ON/OFF before the 1st time step)
% INI_STATE [NG x 1] - initial states of generators (1-commited, 0-not commited)
% INI_STATE_NUM      - position (column) of vector INI_STATE in the list of states
INI_STATE = (GSTATINI > 0);
[I, INI_STATE_NUM]= ismember(INI_STATE',LIST_STATES','rows');
HOUR_inf = zeros(length(max_array_size),1);