function [GENERATION,F_BASE,P_CAPT] = dispatch(CURRENT_STATE,u_CAPT_g_t,Y_CAPT_g_t,lb,ub,GMIN,GMAX,P_MIN_CAPT_g_t,P_MAX_CAPT_g_t,DEMAND,HOUR,GEN_ORDER,m_A_g,m_B_g,m_C_g,c_A_g,...
    c_B_g,c_C_g,DISPATCH_METHOD,CC_FIXED_g,CC_OP_g,ECO2_g,GINC,GMIN_0,GMAX_0)
% --------------------------------------------------------------------------------------------------------------
% For the given state, calculates the MW output for each commited generator
% Generators are dispatched in a merit order (first the least expensive, 
% last the most expensive)
% Note: GEN_ORDER is based on No Load Cost, Fuel Cost and Incremental costs.
% OUTPUT:
% GENERATION [NG x 1] - vector of power output for each generator
% ---------------------------------------------------------------------------------------------------------------
Y_CAPT_g_t = Y_CAPT_g_t(HOUR,:)';
GENERATION = GMIN.*CURRENT_STATE; % first set the output for each commited generator to their minimal stable generation
P_CAPT     = CURRENT_STATE.*P_MIN_CAPT_g_t(HOUR,:)';
P_CAPT_MIN = CURRENT_STATE.*P_MIN_CAPT_g_t(HOUR,:)';
LOAD = DEMAND(HOUR) - sum(GENERATION)+ sum(P_CAPT_MIN); % then reduce the load for the total minimal generation

if DISPATCH_METHOD == 0  % piece-wise quadratic approximation
    for K = 1:length(CURRENT_STATE)
        L = GEN_ORDER(K);
        GENERATION(L) = GENERATION(L) + min(GMAX(L)-GMIN(L),LOAD)*CURRENT_STATE(L);
        if GENERATION(L) == 0
            CURRENT_STATE(L) = 0;
        end
        if GENERATION(L) >= 0.8*GMAX_0(L)
            F_BASE(L) = (GENERATION(L)*m_C_g(L) + c_C_g(L))*CURRENT_STATE(L);
        elseif GENERATION(L) >= 0.6*GMAX_0(L)
            F_BASE(L) = (GENERATION(L)*m_B_g(L) + c_B_g(L))*CURRENT_STATE(L);
        else
            F_BASE(L) = (GENERATION(L)*m_A_g(L) + c_A_g(L))*CURRENT_STATE(L);
        end
        P_CAPT(L) = round(CC_FIXED_g(L) + (CC_OP_g(L).*(((GENERATION(L).*m_A_g(L) + c_A_g(L)).*CURRENT_STATE(L)).*ECO2_g(L)).*Y_CAPT_g_t(L))).*CURRENT_STATE(L);
        LOAD = LOAD - min(GMAX(L)-GMIN(L),LOAD)*CURRENT_STATE(L) + P_CAPT(L) - P_CAPT_MIN(L);
    end
else
    if DISPATCH_METHOD == 3 % quick dispatch linear
        for K = 1:length(CURRENT_STATE)                 % note that CURRENT_STATE is the feasible one, ie.  demand may be supplied by commited generators
            L = GEN_ORDER(K);                           % GEN_ORDER is the merit list of dispatching generators
            if GENERATION(L) == 0
                CURRENT_STATE(L) = 0;
            end
            GENERATION(L) = GENERATION(L) + min(GMAX(L)-GMIN(L),LOAD)*CURRENT_STATE(L); % increase the power of the next generator in the list either to their max. or to match the load
            P_CAPT(L) = round(CC_FIXED_g(L) + (CC_OP_g(L).*(((GENERATION(L).*(1./GINC(L)))*CURRENT_STATE(L)).*ECO2_g(L)).*Y_CAPT_g_t(L))).*CURRENT_STATE(L);
            LOAD = LOAD - min(GMAX(L)-GMIN(L),LOAD)*CURRENT_STATE(L) + P_CAPT(L) - P_CAPT_MIN(L); % whenever generation increases, load reduces, until they match
            F_BASE(L) = GENERATION(L).*(1./GINC(L));
        end
    end
end