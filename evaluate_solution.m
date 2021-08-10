function evaluate_solution(NT,BEST_PATH,LIST_STATES,GMIN,GMAX,D_t,W_t,W_on_t,W_of_t,W_c_t,S_t,S_c_t,W_t_Original,W_on_t_Original,W_of_t_Original,...
    S_t_Original,Tot_c_t,Net_D_t_Original,DEMAND,Price,RES_UP,RES_DN,GEN_ORDER,GNLC,GFC,GINC,m_A_g,m_B_g,m_C_g,c_A_g,c_B_g,c_C_g,GSC,SU_FUEL_COLD_g,SU_CO2_COLD_g,...
    GSDC,SD_FUEL_g,SD_CO2_g,INI_STATE,NG,GMINUP,GMINDOWN,...
    GRAMPUP,GRAMPDOWN,COEF_A,COEF_B,COEF_C,CC_FIXED_g,CC_OP_g,DISPATCH_METHOD,FLEXIBLE_CO2_CAPTURE_OPTION,GSTATINI,...
    GMINUP_Original,GMINDOWN_Original,START_UP_COST_METHOD,GCSTIME,GSH,ALPHA,BETA,TAU,ECO2_g,CoC_g,FLEX_g,Y_MAX_CAPT_g,...
    STK_g,VOM_g,P_MIN_CAPT_g,P_MAX_CAPT_g,Y_CAPT_g_t,u_CAPT_g_t,P_MIN_CAPT_g_t,...
    P_MAX_CAPT_g_t,CC_VOM_g,CC_CSOLV_g,CC_SOLVD_g,CC_SOLVD_TH_g,CC_TRANS_g,RAMP_UP_g,RAMP_DN_g,GMIN_0,GMAX_0,RELAXATION,EXP_COST_OPTION,STACKING_OPTION,RAMP_COSTS_INCLUDED,...
    On_subsidy,Off_subsidy,Uplift_c,Uplift_b,Uplift_VOLL,Uplift_SUM_G,Uplift_k,D_MIN_t,W_on_c_inert,W_of_c_inert,S_c_inert,S_c_inert_s,W_on_c_am_inert,W_of_c_am_inert,S_c_am_inert,...
    S_c_am_inert_s,W_on_c_feas,W_on_c_am_feas,W_of_c_feas,W_of_c_am_feas,S_c_feas,S_c_am_feas,S_c_feas_s,S_c_am_feas_s,HOUR_inf)
% --------------------------------------------------------------------------------------------------------------
% For the given set BEST_PATH of the states in each time step, calculates production costs, transition costs
% total costs etc.
% In the end, this function calls function to print the results in a tabulated form (optional).
% This function follows pretty much the same fashion as the main procedure and will not be described in details.
% The only difference is that now the optimal path is known, so there is no need for searching.
% For each state in the path this procedure calculates the costs and finally calls the printing routine.
% ---------------------------------------------------------------------------------------------------------------
GEN_START_SHUT_COST1  = zeros(NG,NT);
CO2_SU_SD1            = zeros(NG,NT);
FUEL_SU_SD1           = zeros(NG,NT);
COUNT_SU_SD1          = zeros(NG,NT);
GEN_PRODUCTION1       = zeros(NG,NT);
CURRENT_STATE1        = zeros(NG,NT);
MW_RAMP1              = zeros(NG,NT);
X_t_11                = zeros(NG,NT);
%
FUEL_COMP1            = zeros(NG,NT);
PROD_COST1            = zeros(NG,NT);
COST_FUEL1            = zeros(NG,NT);
COST_CO21             = zeros(NG,NT);
COST_GVOM1            = zeros(NG,NT);
COST_RAMP1            = zeros(NG,NT);
CO2_PROD1             = zeros(NG,NT);
CO2_CAPT1             = zeros(NG,NT);
FCOST1                = zeros(NT,1);
%
GENERATING_COST1      = zeros(NT,1);
GEN_PRODUCTION        = zeros(NG,1);
P_CAPT                = zeros(NG,1);
P_R_SP_UP_t=zeros(NT,NG);
P_R_SP_DN_t=zeros(NT,NG);
P_CAPT1=zeros(NG,NT);
%
X  = GSTATINI;
for HOUR = 1:NT
    PREV_STATES_NUM = BEST_PATH(HOUR);
    FEASIBLE_STATES_NUM = BEST_PATH(HOUR+1);
    X_PREV = X;
    if HOUR==1 && PREV_STATES_NUM == 0
        PREVIOUS_STATE = INI_STATE;
    else
        PREVIOUS_STATE = LIST_STATES(:,PREV_STATES_NUM);
        for G = 1:length(PREVIOUS_STATE)
            g = GEN_ORDER(G);
            if X_PREV(g) > 0
                PREVIOUS_STATE(g) = 1;
            end
        end
    end
    CURRENT_STATE  = LIST_STATES(:,FEASIBLE_STATES_NUM);
    PRODUCTION_PREV = GEN_PRODUCTION;
    %    
    [GEN_PRODUCTION,PROD_COST,P_CAPT,FUEL_COMP] = production(CURRENT_STATE,u_CAPT_g_t,PREVIOUS_STATE,GMIN,GMAX,D_t,W_t,W_on_t,W_of_t,S_t,DEMAND,RES_UP,RES_DN,HOUR,...
        GNLC,GFC,GINC,m_A_g,m_B_g,m_C_g,c_A_g,c_B_g,c_C_g,NG,GRAMPUP,GRAMPDOWN,PRODUCTION_PREV,GEN_ORDER,COEF_A,COEF_B,COEF_C,DISPATCH_METHOD,ECO2_g,CoC_g,Y_MAX_CAPT_g,STK_g,...
        VOM_g,Price,CC_FIXED_g,CC_OP_g,P_MIN_CAPT_g,P_MAX_CAPT_g,P_MIN_CAPT_g_t,P_MAX_CAPT_g_t,Y_CAPT_g_t,CC_VOM_g,CC_CSOLV_g,CC_SOLVD_g,CC_SOLVD_TH_g,CC_TRANS_g,GMIN_0,GMAX_0);
    %
    % Check to see if net demand equals the sum of thermal generators
    % Ramps that are violated are assigned to flexible generators
    GMINUP = GMINUP_Original;
    GMINDOWN = GMINDOWN_Original;
    %
    STATE_DIFF = CURRENT_STATE - PREVIOUS_STATE;
    [X,SUCCESS] = check_up_down_time(CURRENT_STATE,PREVIOUS_STATE,X_PREV,GMINUP,GMINDOWN,NG,RELAXATION);
    %
    if RELAXATION==1
        if SUCCESS == 0
            GMINUP(:)   = 0;
            GMINDOWN(:) = 0;
            [X,~]       = check_up_down_time(CURRENT_STATE,PREVIOUS_STATE,X_PREV,GMINUP,GMINDOWN,NG,RELAXATION);
        end
    end
    %
    GMINUP = GMINUP_Original;
    GMINDOWN = GMINDOWN_Original;
    %
    P_RAMP = zeros(NG,1);
    RAMP_VIOLATIONS = DEMAND(HOUR) - sum(GEN_PRODUCTION);
    if RAMP_VIOLATIONS > 1
        for G = 1:length(CURRENT_STATE)
            g = GEN_ORDER(G);
            P_RAMP(g) = P_RAMP(g) + min((GMAX(g)-(GEN_PRODUCTION(g)+P_CAPT(g)).*CURRENT_STATE(g)).*FLEX_g(g),RAMP_VIOLATIONS);
            RAMP_VIOLATIONS = RAMP_VIOLATIONS - min((GMAX(g)-(GEN_PRODUCTION(g)+P_CAPT(g)).*CURRENT_STATE(g)).*FLEX_g(g),RAMP_VIOLATIONS);
            if P_RAMP(g) > 0
                CURRENT_STATE(g) = 1;
                STATE_DIFF(g) = CURRENT_STATE(g) - PREVIOUS_STATE(g);
                if (X_PREV(g) >= 1) && (CURRENT_STATE(g) == 1)
                    X(g) = X_PREV(g) + 1;
                elseif (X_PREV(g) <= -GMINDOWN(g)) && (CURRENT_STATE(g) == 1)
                    X(g) = 1;
                elseif (X_PREV(g) <= -1) && (CURRENT_STATE(g) == 0)
                    X(g) = X_PREV(g) - 1;
                elseif (X_PREV(g) >= GMINUP(g)) && (CURRENT_STATE(g) == 0)
                    X(g) = -1;
                end
            end
            GEN_PRODUCTION(g) = GEN_PRODUCTION(g) + P_RAMP(g);
        end
    end
    %
    FUEL_COMP = FUEL_COMP + (1./GINC).*P_RAMP.*CURRENT_STATE;
    COST_FUEL = GFC.*FUEL_COMP;
    COST_CO2 = ECO2_g.*CoC_g.*(1-Y_CAPT_g_t(HOUR,:)').*FUEL_COMP;
    CO2_PROD = ECO2_g.*(1-Y_CAPT_g_t(HOUR,:)').*FUEL_COMP;
    CO2_CAPT = ECO2_g.*(Y_CAPT_g_t(HOUR,:)').*FUEL_COMP;
    COST_GVOM = VOM_g.*GEN_PRODUCTION;
    COST_CC_VOM = CO2_CAPT.*CC_VOM_g;
    COST_CC_SOLV = GMAX.*(1./GINC).*ECO2_g.*Y_MAX_CAPT_g.*CC_CSOLV_g.*(((CC_SOLVD_g - CC_SOLVD_TH_g).*Y_CAPT_g_t(HOUR,:)') + CC_SOLVD_TH_g);
    COST_CC_TRANS = CO2_CAPT.*CC_TRANS_g;
    PROD_COST = COST_FUEL + COST_CO2 + COST_GVOM + COST_CC_VOM + COST_CC_SOLV + COST_CC_TRANS;
    %
    if START_UP_COST_METHOD == 1   % start-up costs are constant and equal to cold start cost
        GEN_START_SHUT_COST = (STATE_DIFF > 0) .* (GSC + SU_FUEL_COLD_g.*GFC + SU_CO2_COLD_g.*CoC_g.*(1-Y_CAPT_g_t(HOUR,:)'));
    elseif START_UP_COST_METHOD == 2 % cold/hot start-up costs
        GEN_START_SHUT_COST = ((STATE_DIFF > 0) & (-X_PREV >= (GMINDOWN + GCSTIME))) .* (GSC...
            + SU_FUEL_COLD_g.*GFC + SU_CO2_COLD_g.*CoC_g.*(1-Y_CAPT_g_t(HOUR,:)'));  % cold start-up cost
        GEN_START_SHUT_COST = GEN_START_SHUT_COST + ((STATE_DIFF > 0) & (-X_PREV <  (GMINDOWN + GCSTIME))) .* (GSH +...
            +SU_FUEL_COLD_g.*GFC + SU_CO2_COLD_g.*CoC_g.*(1-Y_CAPT_g_t(HOUR,:)'));  % hot start-up cost
    else
        if  EXP_COST_OPTION==0 % all costs inside the brackets
            GEN_START_SHUT_COST = (STATE_DIFF > 0) .* ((GSC + SU_FUEL_COLD_g.*GFC + SU_CO2_COLD_g.*CoC_g.*(1-Y_CAPT_g_t(HOUR,:)')) .* ...
                (1-exp(X_PREV ./ TAU)));
        elseif EXP_COST_OPTION==1 % % start-up fixed cost outside of the brackets
            GEN_START_SHUT_COST = (STATE_DIFF > 0) .* (GSC + ((SU_FUEL_COLD_g.*GFC + SU_CO2_COLD_g.*CoC_g.*(1-Y_CAPT_g_t(HOUR,:)'))) .* ...
                (1-exp(X_PREV ./ TAU)));
        elseif EXP_COST_OPTION==2 % "cold" alpha-beta approach
            GEN_START_SHUT_COST = (STATE_DIFF > 0) .* (ALPHA + (BETA + SU_FUEL_COLD_g.*GFC + SU_CO2_COLD_g.*CoC_g.*(1-Y_CAPT_g_t(HOUR,:)')) .* ...
                (1-exp(X_PREV./ TAU)));
        else % "hot" alpha-beta approach
            GEN_START_SHUT_COST = (STATE_DIFF > 0) .* (ALPHA + SU_FUEL_COLD_g.*GFC + SU_CO2_COLD_g.*CoC_g.*(1-Y_CAPT_g_t(HOUR,:)') + BETA .* ...
                (1-exp(X_PREV./ TAU)));
        end
    end
    %exponential start-up fuel
    FUEL_SU_SD = (STATE_DIFF > 0) .* ((SU_FUEL_COLD_g) .* (1-exp(X_PREV ./ TAU)));
    CO2_SU_SD = (STATE_DIFF > 0) .* ((SU_CO2_COLD_g.*(1-Y_CAPT_g_t(HOUR,:)')) .* (1-exp(X_PREV ./ TAU)));
    COUNT_SU_SD = (STATE_DIFF > 0) .* X_PREV;
    %
    GEN_START_SHUT_COST = GEN_START_SHUT_COST + (STATE_DIFF < 0 ) .* (GSDC + + SD_FUEL_g.*GFC + SD_CO2_g.*CoC_g.*(1-Y_CAPT_g_t(HOUR,:)'));       % shut down cost
    CO2_SU_SD = CO2_SU_SD + (STATE_DIFF  < 0) .* (SD_CO2_g.*(1-Y_CAPT_g_t(HOUR,:)'));
    %shut-down fuel
    FUEL_SU_SD = FUEL_SU_SD + (STATE_DIFF  < 0) .* (SD_FUEL_g);
    %
    GEN_START_SHUT_COST(isnan(GEN_START_SHUT_COST)) = 0;
    FUEL_SU_SD(isnan(FUEL_SU_SD)) = 0;
    CO2_SU_SD(isnan(CO2_SU_SD)) = 0;
    %
    COST_RAMP = (STATE_DIFF == 0) .* ((GEN_PRODUCTION.*CURRENT_STATE - PRODUCTION_PREV.*PREVIOUS_STATE) > 0) .* abs(GEN_PRODUCTION.*CURRENT_STATE - PRODUCTION_PREV.*PREVIOUS_STATE).*RAMP_UP_g + ...
        (STATE_DIFF == 0) .* ((PRODUCTION_PREV.*PREVIOUS_STATE - GEN_PRODUCTION.*CURRENT_STATE)> 0) .* abs(PRODUCTION_PREV.*PREVIOUS_STATE - GEN_PRODUCTION.*CURRENT_STATE).*RAMP_DN_g;
    MW_RAMP = (STATE_DIFF == 0) .* ((GEN_PRODUCTION.*CURRENT_STATE - PRODUCTION_PREV.*PREVIOUS_STATE) > 0) .* (GEN_PRODUCTION.*CURRENT_STATE - PRODUCTION_PREV.*PREVIOUS_STATE) + ...
        (STATE_DIFF == 0) .* ((PRODUCTION_PREV.*PREVIOUS_STATE - GEN_PRODUCTION.*CURRENT_STATE)> 0) .* (PRODUCTION_PREV.*PREVIOUS_STATE - GEN_PRODUCTION.*CURRENT_STATE);
    %
            if RAMP_COSTS_INCLUDED==1
                if HOUR == 1
                    TOTAL_COST = sum(PROD_COST) + sum(GEN_START_SHUT_COST) + sum(COST_RAMP);
                else
                    TOTAL_COST = sum(PROD_COST) + sum(GEN_START_SHUT_COST) + sum(COST_RAMP) + FCOST1(HOUR-1);
                end
            else
                if HOUR == 1
                    TOTAL_COST = sum(PROD_COST) + sum(GEN_START_SHUT_COST);
                else
                    TOTAL_COST = sum(PROD_COST) + sum(GEN_START_SHUT_COST) + FCOST1(HOUR-1);
                end
            end
    %
    FCOST1(HOUR)                = TOTAL_COST;
    GEN_PRODUCTION1(:,HOUR)     = GEN_PRODUCTION;
    P_CAPT1(:,HOUR)             = P_CAPT;
    CURRENT_STATE1(:,HOUR)      = CURRENT_STATE;
    MW_RAMP1(:,HOUR)            = MW_RAMP;
    X_t_11(:,HOUR)              = X_PREV;
    FUEL_COMP1(:,HOUR)          = FUEL_COMP;
    PROD_COST1(:,HOUR)          = PROD_COST;
    COST_FUEL1(:,HOUR)          = COST_FUEL;
    COST_CO21(:,HOUR)           = COST_CO2;
    CO2_PROD1(:,HOUR)           = CO2_PROD;
    CO2_CAPT1(:,HOUR)           = CO2_CAPT;
    COST_GVOM1(:,HOUR)          = COST_GVOM;
    COST_RAMP1(:,HOUR)          = COST_RAMP;
    COST_CC_VOM1(:,HOUR)        = COST_CC_VOM;
    COST_CC_SOLV1(:,HOUR)       = COST_CC_SOLV;
    COST_CC_TRANS1(:,HOUR)      = COST_CC_TRANS;
    %
    GEN_START_SHUT_COST1(:,HOUR)= GEN_START_SHUT_COST;
    CO2_SU_SD1(:,HOUR) = CO2_SU_SD;
    FUEL_SU_SD1(:,HOUR) = FUEL_SU_SD;
    COUNT_SU_SD1(:,HOUR) =  COUNT_SU_SD;
    %
    GENERATING_COST1(HOUR)      = sum(PROD_COST);
    %
end   % HOUR = 1:NT
GEN_START_SHUT_COST_TOTAL = sum(GEN_START_SHUT_COST1).';
%
if FLEXIBLE_CO2_CAPTURE_OPTION == 1
    P_R_SP_UP_t_NC = (ones(NT,1)*GMAX' - ones(NT,1)*CC_FIXED_g').*CURRENT_STATE1' - GEN_PRODUCTION1';
else
    P_R_SP_UP_t_NC = (ones(NT,1)*GMAX').*CURRENT_STATE1' - GEN_PRODUCTION1' - P_CAPT1';
end
%
P_R_SP_UP_t_RC = ((ones(NT,1)*GRAMPUP') + P_CAPT1' - ones(NT,1)*CC_FIXED_g').*CURRENT_STATE1';
P_R_SP_DN_t_NC = GEN_PRODUCTION1' - (ones(NT,1)*GMIN' - P_MIN_CAPT_g_t).*CURRENT_STATE1';
P_R_SP_DN_t_RC = (ones(NT,1)*GRAMPDOWN').*CURRENT_STATE1';
%
for i = 1:NT
    for j = 1:NG
        P_R_SP_UP_t(i,j) = min([P_R_SP_UP_t_NC(i,j),P_R_SP_UP_t_RC(i,j)]);
        P_R_SP_DN_t(i,j) = min([P_R_SP_DN_t_NC(i,j),P_R_SP_DN_t_RC(i,j)]);
    end
end
%
RESERVE_UP_MAX_NC = sum(P_R_SP_UP_t_NC,2);
RESERVE_DN_MAX_NC = sum(P_R_SP_DN_t_NC,2);
%
RESERVE_UP_MAX = sum(P_R_SP_UP_t,2);
RESERVE_DN_MAX = sum(P_R_SP_DN_t,2);
%
Status = CURRENT_STATE1';
Starts = GEN_START_SHUT_COST1';
%
Days = NT/24;
Op_hours = zeros(Days,NG);
COST_GEN_START_SHUT = zeros(Days,NG);
l = 1:24:NT;
%
for k = 1:Days
    for j = 1:NG
        Op_hours(k,j) = sum(Status(l(k):k*24,j));
        COST_GEN_START_SHUT(k,j) = sum(Starts(l(k):k*24,j))./(max(2,Op_hours(k,j)).*GMAX_0(j));
    end
end
COST_GEN_START_SHUT(isnan(COST_GEN_START_SHUT)) = 0;
COST_GEN_START_SHUT = kron(COST_GEN_START_SHUT,ones(24,1));
%% Marginal calculations
Merit_Order  = cumsum(GMAX);
%
if STACKING_OPTION==2 % short-run marginal cost (SRMC)
    MERIT_COST       = 0*GNLC + GFC.*(1./GINC) + VOM_g + ECO2_g.*(1./GINC).*CoC_g.*(1-Y_MAX_CAPT_g) + ECO2_g.*Y_MAX_CAPT_g.*(1./GINC).*(CC_VOM_g + CC_TRANS_g) - STK_g;
else %average full-load cost
    MERIT_COST = (0*GNLC + GFC.*GMAX.*GINC/1000)./GMAX;
end
%
Marginal_Generator = CURRENT_STATE1';
Marginal_Power = GEN_PRODUCTION1';
Marginal_Fuel = FUEL_COMP1';
%
MC = zeros(NT,1);
SC = zeros(NT,1);
UL = zeros(NT,1);
UL_e = zeros(NT,1);
UL_max = zeros(NT,1);
CR = zeros(NT,1);
%
for HOUR = 1:NT
    MG1 = find(Marginal_Generator(HOUR,:));
    MG = MG1(end);
    P_MARGINAL1 = Marginal_Power(HOUR,:);
    F_MARGINAL1 = Marginal_Fuel(HOUR,:);
    P_MARGINAL = P_MARGINAL1(MG);
    F_MARGINAL = F_MARGINAL1(MG);
    %
    if STACKING_OPTION==2 % short-run marginal cost (SRMC)
        MC(HOUR) = (0*GNLC(MG) + GFC(MG).*F_MARGINAL + VOM_g(MG) + ECO2_g(MG).*F_MARGINAL.*CoC_g(MG).*(1-Y_MAX_CAPT_g(MG)) + ECO2_g(MG).*Y_MAX_CAPT_g(MG).*F_MARGINAL.*(CC_VOM_g(MG)...
            + CC_TRANS_g(MG)) - STK_g(MG))./P_MARGINAL;
    else %average full-load cost
        MC(HOUR) = ((0*GNLC(MG) + GFC(MG).*GMAX(MG).*(1/F_MARGINAL)/1000)./GMAX(MG))./P_MARGINAL;
    end
    %
    SC(HOUR) = COST_GEN_START_SHUT(HOUR,MG);
    UL(HOUR) = ((Uplift_k*DEMAND(HOUR))/(Uplift_SUM_G-(DEMAND(HOUR))));
    UL_e(HOUR) = Uplift_VOLL*exp(-Uplift_b*((Uplift_SUM_G-DEMAND(HOUR))/Uplift_SUM_G));
    UL_max(HOUR) = max(UL(HOUR),UL_e(HOUR));
    CR(HOUR) = -On_subsidy*exp(Uplift_c*((W_t_Original(HOUR)-D_t(HOUR)))/D_t(HOUR));% I guess On_subsidy=50 ï¿½/MWh is the subsidy level for onshore wind (VA)
    Price(HOUR) =  MC(HOUR) + SC(HOUR) + UL_max(HOUR) + CR(HOUR);
end
Price(Price<-Off_subsidy) = -Off_subsidy;
%% Output generation
v9_GENERATION = GEN_PRODUCTION1';
v9_P_CAPT = P_CAPT1';
v9_P_NET = GEN_PRODUCTION1' - P_CAPT1';
v9_STATUS = CURRENT_STATE1';
v9_COUNT_SU_SD = COUNT_SU_SD1';
v9_MW_RAMP = MW_RAMP1';
v9_X_PREV = X_t_11';
%
v9_GEN_START_SHUT_COST = GEN_START_SHUT_COST1';
v9_RESERVE_UP = P_R_SP_UP_t;
v9_RESERVE_DN = P_R_SP_DN_t;
%
v9_FUEL_SU_SD = FUEL_SU_SD1';
v9_FUEL_COMP = FUEL_COMP1';
v9_FUEL_ALL = v9_FUEL_SU_SD + v9_FUEL_COMP;
%
v9_PROD_COST = PROD_COST1';
v9_COST_FUEL = COST_FUEL1';
v9_COST_CO2 = COST_CO21';
v9_COST_GVOM = COST_GVOM1';
v9_COST_RAMP = COST_RAMP1';
%
v9_CO2_SU_SD = CO2_SU_SD1';
v9_CO2_PROD = CO2_PROD1';
v9_CO2_CAPT = CO2_CAPT1';
v9_CO2_ALL = v9_CO2_SU_SD + v9_CO2_PROD;
v9_CO2_EMISSION_INTENSITY1 = sum(v9_CO2_ALL(:))/sum(D_t(:));
v9_COST_SUM_START_STOP = sum(v9_GEN_START_SHUT_COST,2);
v9_COST_SUM_PROD = sum(v9_PROD_COST,2);
v9_COST_SUM_FUEL = sum(v9_COST_FUEL,2);
v9_COST_SUM_CO2 = sum(v9_COST_CO2,2);
v9_COST_SUM_RAMP = sum(v9_COST_RAMP,2);
%
v9_CO2_SUM_START_STOP = sum(v9_CO2_SU_SD,2);
v9_CO2_SUM_PROD = sum(v9_CO2_PROD,2);
v9_CO2_SUM_CAPT = sum(v9_CO2_CAPT,2);
v9_CO2_SUM_ALL = sum(v9_CO2_ALL,2);
%
if RAMP_COSTS_INCLUDED==1
    v9_COST_ALL = v9_GEN_START_SHUT_COST + v9_PROD_COST + v9_COST_RAMP;
    v9_TOTAL_COST = v9_COST_SUM_START_STOP + v9_COST_SUM_PROD + v9_COST_SUM_RAMP;
else
    v9_COST_ALL = v9_GEN_START_SHUT_COST + v9_PROD_COST;
    v9_TOTAL_COST = v9_COST_SUM_START_STOP + v9_COST_SUM_PROD;
end
%
v9_COST_SUM_ALL = sum(v9_COST_ALL,2);
%
v9_count_curt_on_inert          = W_on_c_inert';
v9_count_curt_of_inert          = W_of_c_inert';
v9_count_curt_s_inert_l         = S_c_inert';
v9_count_curt_s_inert_s         = S_c_inert_s';
v9_am_curt_on_inert             = W_on_c_am_inert';
v9_am_curt_of_inert             = W_of_c_am_inert';
v9_am_curt_s_inert_l            = S_c_am_inert';
v9_am_curt_s_inert_s            = S_c_am_inert_s';
v9_count_curt_on_feas           = W_on_c_feas';
v9_count_curt_of_feas           = W_of_c_feas';
v9_count_curt_s_feas_l          = S_c_feas';
v9_count_curt_s_feas_s          = S_c_feas_s';
v9_am_curt_on_feas              = W_on_c_am_feas';
v9_am_curt_of_feas              = W_of_c_am_feas';
v9_am_curt_s_feas_l             = S_c_am_feas';
v9_am_curt_s_feas_s             = S_c_am_feas_s';
if RELAXATION==1
    v9_infeasible_hours         = find(HOUR_inf);
    csvwrite('v9_infeasible_hours.csv',v9_infeasible_hours);
end
%
v9_SYSTEM = [Price,v9_COST_SUM_ALL,v9_COST_SUM_START_STOP,v9_COST_SUM_PROD,v9_COST_SUM_FUEL,v9_COST_SUM_CO2,v9_COST_SUM_RAMP,v9_CO2_SUM_ALL,...
    v9_CO2_SUM_START_STOP,v9_CO2_SUM_PROD,v9_CO2_SUM_CAPT,D_t,W_t_Original,W_on_t_Original,W_of_t_Original,S_t_Original,Net_D_t_Original,...
    W_t,W_on_t,W_of_t,S_t,DEMAND,Tot_c_t,S_c_t,W_c_t,v9_count_curt_s_inert_l,v9_count_curt_on_inert,v9_count_curt_of_inert,v9_count_curt_s_inert_s,v9_am_curt_s_inert_l,...
    v9_am_curt_on_inert,v9_am_curt_of_inert,v9_am_curt_s_inert_s,v9_count_curt_s_feas_l,v9_count_curt_on_feas,v9_count_curt_of_feas,v9_count_curt_s_feas_s,v9_am_curt_s_feas_l,v9_am_curt_on_feas,...
    v9_am_curt_of_feas,v9_am_curt_s_feas_s];
%
v9_SYSTEM_h = {'Price (?/MWhe)','Total cost (?/h)','Start/stop cost (?/h)','Production cost (?/h)','Fuel cost (?/h)','CO2 cost (?/h)','Ramping Cost (?/h)',...
    'Total CO2 emissions (tCO2/h)','Start/stop emissions (tCO2/h)','Production CO2 emissions (tCO2/h)','Captured CO2 emissions (tCO2/h)','D_t (MWe)',...
    'W_t original (MWe)','W_t onshore original (MWe)','W_t offshore original (MWe)','S_t solar original (MWe)','Net demand origninal (MWe)',...
    'W_t (MWe)','W_t onshore (MWe)','W_t offshore (MWe)','S_t solar (MWe)',...
    'Net demand (MWe)','Initial curtailment (MWe)','Solar curtailment (MWe)','Wind curtailment (MWe)','Number of large solar curtailment for inertia',...
    'Number of onshore wind curtailment for inertia','Number of offshore wind curtailment for inertia','Number of small solar curtailment for inertia','Amount of large solar curtailment for inertia',...
    'Amount of onshore wind curtailment for inertia','Amount of offshore wind curtailment for inertia','Amount of small solar curtailment for inertia','Number of large solar curtailment for feasibility',...
    'Number of onshore wind curtailment for feasibility','Number of offshore wind curtailment for feasibility','Number of small solar curtailment for feasibility','Amount of large solar curtailment for feasibility',...
    'Amount of onshore wind curtailment for feasibility','Amount of offshore wind curtailment for feasibility','Amount of small solar curtailment for feasibility'};
csvwrite_with_headers('v9_SYSTEM.csv',v9_SYSTEM,v9_SYSTEM_h)
%
% print_results(BEST_PATH,LIST_STATES,INI_STATE,NT,NG,GMIN,GMAX,DEMAND,FCOST1,GENERATING_COST1,GEN_PRODUCTION1,PROD_COST1,GEN_START_SHUT_COST1,DETAIL_PRINT_FLAG);
v9_GENERATOR_h = {'NUC 1-2','NUC 3-4','NUC 5-6','NUC 7-8','CCGT+PCC 1','CCGT+PCC 2','CCGT+PCC 3','CCGT+PCC 4','CCGT 1','CCGT 2','CCGT 3','CCGT 4','CCGT 5',...
    'CCGT 6','CCGT 7','CCGT 8','CCGT 9','CCGT 10','CCGT 11','CCGT 12','CCGT 13','CCGT 14','CCGT 15','CCGT 16','CCGT 17','CCGT 18','CCGT 19','CCGT 20',...
    'CCGT 21','CCGT 22','CCGT 23','CCGT 24','CCGT 25','CCGT 26','CCGT 27','CCGT 28','CCGT 29','CCGT 30','CCGT 31','CCGT 32','CCGT 33','CCGT 34','CCGT 35',...
    'CCGT 36','CCGT 37','CCGT 38','CCGT 39','CCGT 40','OCGT 1-2','OCGT 3-4','OCGT 5-6','OCGT 7-8','OCGT 9-10','OCGT 11-12','OCGT 13-14','OCGT 15-16',...
    'OCGT 17-18','OCGT 19-20','OCGT 21-22','OCGT 23-24','OCGT 25-26','OCGT 27-28','OCGT 29-30','OCGT 31-32','OCGT 33-34','OCGT 35-36','OCGT 37-38','OCGT 39-40'};
%
csvwrite_with_headers('v9_GENERATION.csv',v9_GENERATION,v9_GENERATOR_h)
csvwrite_with_headers('v9_P_CAPT.csv',v9_P_CAPT,v9_GENERATOR_h)
csvwrite_with_headers('v9_P_NET.csv',v9_P_NET,v9_GENERATOR_h)
csvwrite_with_headers('v9_STATUS.csv',v9_STATUS,v9_GENERATOR_h)
csvwrite_with_headers('v9_COUNT_SU_SD.csv',v9_COUNT_SU_SD,v9_GENERATOR_h)
csvwrite_with_headers('v9_MW_RAMP.csv',v9_MW_RAMP,v9_GENERATOR_h)
csvwrite_with_headers('v9_X_PREV.csv',v9_X_PREV,v9_GENERATOR_h)
csvwrite_with_headers('v9_FUEL_COMP.csv',v9_FUEL_COMP,v9_GENERATOR_h)
csvwrite_with_headers('v9_FUEL_SU_SD.csv',v9_FUEL_SU_SD,v9_GENERATOR_h)
csvwrite_with_headers('v9_FUEL_ALL.csv',v9_FUEL_ALL,v9_GENERATOR_h)
csvwrite_with_headers('v9_COST_ALL.csv',v9_COST_ALL,v9_GENERATOR_h)
csvwrite_with_headers('v9_GEN_START_SHUT_COST.csv',v9_GEN_START_SHUT_COST,v9_GENERATOR_h)
csvwrite_with_headers('v9_PROD_COST.csv',v9_PROD_COST,v9_GENERATOR_h)
csvwrite_with_headers('v9_COST_FUEL.csv',v9_COST_FUEL,v9_GENERATOR_h)
csvwrite_with_headers('v9_COST_CO2.csv',v9_COST_CO2,v9_GENERATOR_h)
csvwrite_with_headers('v9_COST_GVOM.csv',v9_COST_GVOM,v9_GENERATOR_h)
csvwrite_with_headers('v9_COST_RAMP.csv',v9_COST_RAMP,v9_GENERATOR_h)
csvwrite_with_headers('v9_CO2_ALL.csv',v9_CO2_ALL,v9_GENERATOR_h)
csvwrite_with_headers('v9_CO2_SU_SD.csv',v9_CO2_SU_SD,v9_GENERATOR_h)
csvwrite_with_headers('v9_CO2_PROD.csv',v9_CO2_PROD,v9_GENERATOR_h)
csvwrite_with_headers('v9_CO2_CAPT.csv',v9_CO2_CAPT,v9_GENERATOR_h)
%
csvwrite('v9_CO2_EMISSION_INTENSITY.csv',v9_CO2_EMISSION_INTENSITY1);   %tCO2/MWh
csvwrite('v9_BEST_PATH.csv',BEST_PATH);
csvwrite('v9_FCOST1.csv',FCOST1);
csvwrite('v9_TOTAL_COST.csv',v9_TOTAL_COST);
csvwrite('v9_RESERVE_UP.csv',v9_RESERVE_UP);
csvwrite('v9_RESERVE_DN.csv',v9_RESERVE_DN);
end