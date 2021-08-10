% script: relax UT and DT
for K = 1: length(PREV_STATES_NUM)
    PREVIOUS_STATE = LIST_STATES(:,PREV_STATES_NUM(K));
    %
    GMINUP(:)   = 0;
    GMINDOWN(:) = 0;
    %
    [X,~] = check_up_down_time(CURRENT_STATE,PREVIOUS_STATE,X_PREV(:,K),GMINUP,GMINDOWN,NG,RELAXATION);
    %
    GMINUP = GMINUP_Original;
    GMINDOWN = GMINDOWN_Original;
    %
    STATE_DIFF = CURRENT_STATE - PREVIOUS_STATE;
    %
    % STATE_DIFF = 1  means unit is commited
    % STATE_DIFF = -1 means unit is decommited
    if START_UP_COST_METHOD == 1     % start-up costs are constant and equal to cold start costs
        GEN_START_SHUT_COST(:,K) = (STATE_DIFF > 0) .* (GSC + SU_FUEL_COLD_g.*GFC + SU_CO2_COLD_g.*CoC_g.*(1-Y_CAPT_g_t(HOUR,:)'));
    elseif START_UP_COST_METHOD == 2 % cold/hot start-up costs
        GEN_START_SHUT_COST(:,K) = ((STATE_DIFF > 0) & (-X_PREV(:,K) >= (GMINDOWN + GCSTIME))) .* (GSC...
            + SU_FUEL_COLD_g.*GFC + SU_CO2_COLD_g.*CoC_g.*(1-Y_CAPT_g_t(HOUR,:)'));  % cold start-up cost
        GEN_START_SHUT_COST(:,K) = GEN_START_SHUT_COST(:,K) + ((STATE_DIFF > 0) & (-X_PREV(:,K) <  (GMINDOWN + GCSTIME))) .* (GSH +...
            +SU_FUEL_COLD_g.*GFC + SU_CO2_COLD_g.*CoC_g.*(1-Y_CAPT_g_t(HOUR,:)'));  % hot start-up cost
    else                             % exponential start-up costs
        if  EXP_COST_OPTION==0 % all costs inside the brackets
            GEN_START_SHUT_COST(:,K) = (STATE_DIFF > 0) .* ((GSC + SU_FUEL_COLD_g.*GFC + SU_CO2_COLD_g.*CoC_g.*(1-Y_CAPT_g_t(HOUR,:)')) .* ...
                (1-exp(X_PREV(:,K) ./ TAU)));
        elseif EXP_COST_OPTION==1 % start-up fixed cost outside the brackets
            GEN_START_SHUT_COST(:,K) = (STATE_DIFF > 0) .* (GSC + ((SU_FUEL_COLD_g.*GFC + SU_CO2_COLD_g.*CoC_g.*(1-Y_CAPT_g_t(HOUR,:)'))) .* ...
                (1-exp(X_PREV(:,K) ./ TAU)));
        elseif EXP_COST_OPTION==2 % "cold" alpha-beta approach
            GEN_START_SHUT_COST(:,K) = (STATE_DIFF > 0) .* (ALPHA + (BETA + SU_FUEL_COLD_g.*GFC + SU_CO2_COLD_g.*CoC_g.*(1-Y_CAPT_g_t(HOUR,:)')) .* ...
                (1-exp(X_PREV(:,K)./ TAU)));
        else % "hot" alpha-beta approach
            GEN_START_SHUT_COST(:,K) = (STATE_DIFF > 0) .* (ALPHA + SU_FUEL_COLD_g.*GFC + SU_CO2_COLD_g.*CoC_g.*(1-Y_CAPT_g_t(HOUR,:)') + BETA .* ...
                (1-exp(X_PREV(:,K)./ TAU)));
        end
    end
    %
    GEN_START_SHUT_COST(:,K) = GEN_START_SHUT_COST(:,K) + (STATE_DIFF  < 0) .* (GSDC + SD_FUEL_g.*GFC + SD_CO2_g.*CoC_g.*(1-Y_CAPT_g_t(HOUR,:)'));   % shut down cost added
    GEN_START_SHUT_COST(isnan(GEN_START_SHUT_COST))=0; % Neglecting this line was the source of mischeduling in previous versions
    %
    [GENERATION,PROD_COST,~,~] = production(CURRENT_STATE,u_CAPT_g_t,PREVIOUS_STATE,GMIN,GMAX,D_t,W_t,W_on_t,W_of_t,S_t,DEMAND,RES_UP,RES_DN,HOUR,...
        GNLC,GFC,GINC,m_A_g,m_B_g,m_C_g,c_A_g,c_B_g,c_C_g,NG,GRAMPUP,GRAMPDOWN,PRODUCTION_PREV,GEN_ORDER,COEF_A,COEF_B,COEF_C,DISPATCH_METHOD,ECO2_g,CoC_g,Y_MAX_CAPT_g,STK_g,...
        VOM_g,Price,CC_FIXED_g,CC_OP_g,P_MIN_CAPT_g,P_MAX_CAPT_g,P_MIN_CAPT_g_t,P_MAX_CAPT_g_t,Y_CAPT_g_t,CC_VOM_g,CC_CSOLV_g,CC_SOLVD_g,CC_SOLVD_TH_g,CC_TRANS_g,GMIN_0,GMAX_0);
    %
    COST_RAMP = (STATE_DIFF == 0) .* ((GENERATION.*CURRENT_STATE - PRODUCTION_PREV(:,K).*PREVIOUS_STATE) > 0) .* abs(GENERATION.*CURRENT_STATE - PRODUCTION_PREV(:,K).*PREVIOUS_STATE).*...
        RAMP_UP_g + (STATE_DIFF == 0) .* ((PRODUCTION_PREV(:,K).*PREVIOUS_STATE - GENERATION.*CURRENT_STATE)> 0) .* abs(PRODUCTION_PREV(:,K).*PREVIOUS_STATE - ...
        GENERATION.*CURRENT_STATE).*RAMP_DN_g;
    %
    X_TEMP(:,K) = X;
    PRODUCTION_TEMP(:,K) = GENERATION;
    %
    if RAMP_COSTS_INCLUDED==1
        if HOUR == 1
            TOTAL_COST(K) = sum(PROD_COST) + sum(GEN_START_SHUT_COST(:,K)) + sum(COST_RAMP);
        else
            TOTAL_COST(K) = sum(PROD_COST) + sum(GEN_START_SHUT_COST(:,K)) + sum(COST_RAMP) + FC_PREV(K);
        end
    else
        if HOUR == 1
            TOTAL_COST(K) = sum(PROD_COST) + sum(GEN_START_SHUT_COST(:,K));
        else
            TOTAL_COST(K) = sum(PROD_COST) + sum(GEN_START_SHUT_COST(:,K)) + FC_PREV(K);
        end
    end
end % K
%
[MM,II] = sort(TOTAL_COST(TOTAL_COST ~= 0));
for K = 1:MN
    if isinf(MM(K)) == 0
        COUNTER = COUNTER +1;
        FC(COUNTER,1) = MM(K);
        TR(COUNTER,1:size(TR_PREV,2)) = TR_PREV(II(K),:);
        TR(COUNTER,end) = FEASIBLE_STATES_NUM(J);
        X_CURR(:,COUNTER) = X_TEMP(:,II(K));
        PRODUCTION_CURR(:,COUNTER) = PRODUCTION_TEMP(:,II(K));
    end
end