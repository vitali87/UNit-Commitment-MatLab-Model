function DP_v9() %Vitali Avagyan, June 2017
%------------------------------------------------------------------------------------------------
% UNIT COMMITMENT BY DYNAMIC PROGRAMMING METHOD
%
% Program versions development:
% version v1    - basic DP algorithm with simple quick dispatch (linear, one segment generation curve)
%                 based on the book:
%                 A. J. Wood and B. F. Wollenberg, Power Generation Operation and Control,
%                 1984, John Wiley, New York
% version v2    - minimum up and down times taken into account: this is based on the following
%                 article (with slight modification since cooling time is not taken into account)
%                   C.Li, R.B.Johnson, and A.J.Svoboda: "A new unit commitment method",
%                   IEEE Transactions on Power Systems, Vol. 12, No. 1, 1997, pp. 113ï¿½119, 1997.
% version v3    - generators dispatch based on LINPROG,
%                   for faster processing linprog was replaced with cplex, however it is still slow(VA)
% version v4    - ramp up and down constraints taken into account, dispatching based on LINPROG
% version v4_b  - same as v4 with either simple quick dispatch or linprog based one
% version v5    - enhanced DP method (keep track of several predecessors, not only one)
%                 Version 5 is based on the article:
%                   W.J.Hobbs, et.al., ï¿½An enhanced dynamic programming approach for unit commitment,ï¿½
%                   IEEE Trans. Power Syst., Vol. 3, No. 3, pp. 1201-1205, August 1988.
% version v6    - deterministic spinning reserve added; also QUADPROG option included
% version v7    - shut down cost taken into account; start-up costs can now be cold/hot or
%                 exponential.
% version V8    - some basic check up on data availability added on.
% version v9    - two DP models were merged (Bruce (2015) + Stanojevic (2011)) and
%                 improved (bugs fixed) with enhanced features
% Program can use either:
%       - priority list of generators to be commited
%       - complete enumeration (all possible combinations)
%       - flexible enumaration (combinations for flexible plants?)
% Priority list is based on: 1) no load cost and incremental heat rate
%                               data or
%                            2) short-run marginal cost (SRMC)
% This program uses a forward DP technique and does not take into account
% time resolution different from 1h.
%
% Program features are controlled by a set of switches (flags) provided in
% the input data (DP_input_data_v9) script.
%
% THINGS THAT SHOULD BE IMPROVED/CHANGED IN THE CODE:
% - modify the code so it can deal with any time resolution (not only 1h)
% - since there is a relation between coefficients in linear cost methods and coefficients in
%   quadratic cost method, it is possible to estimate the former one (if not given)
%   based on the latter. This can be done off-line (using, for example, least-square method)
%   and it would offer more versatility to the program (eg. using quick dispatch or priority list
%   even if only quadratic coefficients are given).
% - etc.
%
% Author: Vitali Avagyan, June 2017
% Model was derived and improved from the following two models:
% 1. Alasdair Bruce (2015)
%    http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=7350242
% 2. Vladimir Stanojevic (2011)
%    https://uk.mathworks.com/matlabcentral/fileexchange/32073-unit-commitment-by-dynamic-programming-method
%------------------------------------------------------------------------------------------------
%%
close all;
clear all;                              % clear all varaibles
clc                                     % clean command window
tic                                     % initialise timer
warning off
%------------------------------------------------------------------------------------------------
% This is to add the matlab-cplex files into the directory (if CPLEX is
% installed) for a faster calculation with cplexlp instead of linprog.
% To do that, put LINPROG_OR_CPLEX=0 and add your cplex directory of matlab
% scripts and functions below, e.g.
addpath ('C:\Program Files\IBM\ILOG\CPLEX_Studio1271\cplex\matlab\x64_win64');
savepath;
%------------------------------------------------------------------------------------------------
%% These auxiliary files should be put in the same directory as DP_v9
% (master)file.
DP_input_data_v9;
pre_processing_of_data;
price_block;
CCS_bypass_or_not;
%------------------------------------------------------------------------------------------------
%% main loop of the program - search for optimal commitment by dynamic programming
%------------------------------------------------------------------------------------------------
% Here is the brief explanation of the algorithm: state is a unique
% combination of commited and non-commited generators. Commited generators
% are designated with "1" and non-commited generators with "0".
% For each hour, program finds the potential feasible states. The potential
% feasible states are the states where demand (and reserve) can be supplied
% by the commited generators. If there are no potential feasible states,
% program displays an error message and terminates.
% For each potential feasible state, program takes all feasible states from
% the previous hour and checks if the transition to the current state (in
% the current hour) is possible. If it is not possible, the corresponding
% transition (start-up) cost as well as production cost are set to Inf.
% However, if the transition is possible, transition costs are calculated.
% Production for the current hour is calculated based on the demand
% taking into account production at previous hour (ramp-up and down
% constraints). Finally, total cost is the sum of the transition cost,
% production cost, (ramp cost) and the total cost at the state in previous
% hour. This procedure is repeated for all the states in previous hour.
% Total costs are then sorted and MN of them are saved (this
% is enhancement comparing to the classical dynamic program where only 1
% previous state is saved). If the transition to a state in current hour is
% not possible from any of the states in previous hour, then current state
% is regarded as infeasible and is not considered anymore. If all the
% states in an hour are infeasible, the program has an option to relax the
% up and down times of generators or display the error message and terminate.
%------------------------------------------------------------------------------------------------
for HOUR = 1:NT
    if HOUR == 1
        PREV_STATES_NUM = INI_STATE_NUM;            % Positions (columns) of feasible states in the list of states, for previous hour
        X_PREV  = GSTATINI;                         % number of hours generators are ON (>0) or OFF (<0).
        PRODUCTION_PREV = zeros(size(X_PREV));      % generator outputs
        TR_PREV = PREV_STATES_NUM;                  % transition path matrix
        FC_PREV = 0;                                % cumulative cost vector
    else
        X_PREV = X_CURR;                            % keep the number of gen. working hours for each state (at previous hour) % NB: X_CURR is defined after HOUR==1
        PRODUCTION_PREV = PRODUCTION_CURR;          % and the gen. outputs for previous hour states                           % NB: PRODUCTION_CURR is defined after HOUR==1
        TR_PREV = TR;                               % rows of matrix TR define the transition path                            % NB: TR is defined after HOUR==1
        FC_PREV = FC;                               % save the cumulative cost vector from the previous hour                  % NB: FC is defined after HOUR==1
        PREV_STATES_NUM = TR_PREV(1:COUNTER,end);   % states in the previous hour are given in the last column of TR
    end
    %
    % FEASIBLE_STATES_NUM = positions (columns) of potentially feasible states in the list of states, for current hour.
    [FEASIBLE_STATES_NUM,SUCCESS,DEMAND,W_on_t,W_of_t,S_t,W_c_t,S_c_t,RES_UP,RES_DN,W_on_c_feas,...
    W_on_c_am_feas,W_of_c_feas,W_of_c_am_feas,S_c_feas,S_c_am_feas,S_c_feas_s,S_c_am_feas_s] = find_feasible_states(GMINlst,GMAXlst,DEMAND,HOUR,RES_UP,RES_DN,...
    RES_UP_CCS,RES_DN_CCS,SEQUENTIAL_RES_CURTAILMENT,R_SP_UP_CAPT_t,R_SP_DN_CAPT_t,NT,W_on_t,W_of_t,...
    W_c_t,S_t,S_u,S_c_t,R_SP_UP_D_t,R_SP_DN_D_t,R_SP_UP_W_t,R_SP_DN_W_t,R_SP_UP_S_t,R_SP_DN_S_t,R_SP_UP_MAX,R_SP_DN_MAX,lambda_up,lambda_dn,S_t_large,S_t_small,W_u,D_t,...
    W_on_c_feas,W_on_c_am_feas,W_of_c_feas,W_of_c_am_feas,S_c_feas,S_c_am_feas,S_c_feas_s,S_c_am_feas_s);
    %
    if SUCCESS == 0                                 % if unable to find any feasible state to match demand plus reserve,
        return                                      % quit the program
    end
    %
    MN = min(length(PREV_STATES_NUM),N_PRED);                      % number of predecessors to be examined
    X_CURR = zeros(NG,MN*length(FEASIBLE_STATES_NUM));             % prepare the number of gen. working hours for each state (at current hour)
    PRODUCTION_CURR = zeros(NG,MN*length(FEASIBLE_STATES_NUM));    % prepare generator outputs for each state at current hour
    TR = zeros(MN*length(FEASIBLE_STATES_NUM),HOUR+1);             % prepare transition path matrix
    FC = zeros(MN*length(FEASIBLE_STATES_NUM),1);                  % prepare cumulative cost vector
    COUNTER = 0;
    % take each feasible (current hour) state and...
    for J = 1: length(FEASIBLE_STATES_NUM)
        GEN_START_SHUT_COST = zeros(NG,1);                         % start-up (shut-down) costs
        TOTAL_COST = zeros(1,length(PREV_STATES_NUM));             % total cost (production cost + start up cost + (ramp cost) + total cost of previous hour)
        % X_TEMP - temporarily stores number of gen. working hours for combination of current state and all previous hour states
        X_TEMP = zeros(NG,length(PREV_STATES_NUM));
        % PRODUCTION_TEMP - temporarily stores gen. outputs for combination of current state and all previous hour states
        PRODUCTION_TEMP = zeros(NG,length(PREV_STATES_NUM));
        % take a state (from all feasible states), one by one; let it be
        % CURRENT_STATE...
        CURRENT_STATE  = LIST_STATES(:,FEASIBLE_STATES_NUM(J));
        %----------------------------------------------------------------------------------
        % ... compare it with each feasible state at previous hour
        for K = 1: length(PREV_STATES_NUM)
            if HOUR == 1
                PREVIOUS_STATE = INI_STATE;
            else
                PREVIOUS_STATE = LIST_STATES(:,PREV_STATES_NUM(K));
            end
            %
            GMINUP   = GMINUP_Original;
            GMINDOWN = GMINDOWN_Original;
            %
            % check if the transition from previous state to the current state is possible regarding min up and down times
            [X,SUCCESS] = check_up_down_time(CURRENT_STATE,PREVIOUS_STATE,X_PREV(:,K),GMINUP,GMINDOWN,NG,RELAXATION);
            %
            GMINUP   = GMINUP_Original;
            GMINDOWN = GMINDOWN_Original;
            %
            if SUCCESS==0                                   % if it is not possible,...
                if RELAXATION==1
                    GEN_START_SHUT_COST(:,K) = 0;
                else
                    GEN_START_SHUT_COST(:,K) = inf;         % ...mark the transition cost and...
                end
                CO2_SU_SD(:,K)           = 0;
                PROD_COST                = ones(NG,1)*Inf;  % ... production cost extremely high
                COST_FUEL                = ones(NG,1)*Inf;  % ... extremely high
                COST_CO2                 = ones(NG,1)*Inf;  % ... extremely high
                COST_GVOM                = ones(NG,1)*Inf;  % ... extremely high
                GEN_PRODUCTION           = ones(NG,1)*nan;
                CO2_PROD                 = ones(NG,1)*Inf;  % ... extremely high
                CO2_CAPT                 = ones(NG,1)*nan;
            else                                            % otherwise, calculate the transition cost
                %
                STATE_DIFF = CURRENT_STATE - PREVIOUS_STATE;
                % STATE_DIFF = 1  means unit is commited
                % STATE_DIFF = -1 means unit is decommited
                %
                if START_UP_COST_METHOD == 1     % start-up costs are constant and equal to cold start costs
                    GEN_START_SHUT_COST(:,K) = (STATE_DIFF > 0) .* (GSC + SU_FUEL_COLD_g.*GFC + SU_CO2_COLD_g.*CoC_g.*(1-Y_CAPT_g_t(HOUR,:)'));
                elseif START_UP_COST_METHOD == 2 % cold/hot start-up costs
                    GEN_START_SHUT_COST(:,K) = ((STATE_DIFF > 0) & (-X_PREV(:,K) >= (GMINDOWN + GCSTIME))) .* (GSC...
                        + SU_FUEL_COLD_g.*GFC + SU_CO2_COLD_g.*CoC_g.*(1-Y_CAPT_g_t(HOUR,:)'));  % cold start-up cost
                    GEN_START_SHUT_COST(:,K) = GEN_START_SHUT_COST(:,K) + ((STATE_DIFF > 0) & (-X_PREV(:,K) <  (GMINDOWN + GCSTIME))) .* (GSH +...
                        +SU_FUEL_COLD_g.*GFC + SU_CO2_COLD_g.*CoC_g.*(1-Y_CAPT_g_t(HOUR,:)'));  % hot start-up cost
                else                             % exponential start-up costs
                    if  EXP_COST_OPTION==0       % all costs inside the brackets
                        GEN_START_SHUT_COST(:,K) = (STATE_DIFF > 0) .* ((GSC + SU_FUEL_COLD_g.*GFC + SU_CO2_COLD_g.*CoC_g.*(1-Y_CAPT_g_t(HOUR,:)')) .* ...
                            (1-exp(X_PREV(:,K) ./ TAU)));
                    elseif EXP_COST_OPTION==1    % start-up fixed cost outside the brackets
                        GEN_START_SHUT_COST(:,K) = (STATE_DIFF > 0) .* (GSC + ((SU_FUEL_COLD_g.*GFC + SU_CO2_COLD_g.*CoC_g.*(1-Y_CAPT_g_t(HOUR,:)'))) .* ...
                            (1-exp(X_PREV(:,K) ./ TAU)));
                    elseif EXP_COST_OPTION==2    % "cold" alpha-beta approach
                        GEN_START_SHUT_COST(:,K) = (STATE_DIFF > 0) .* (ALPHA + (BETA + SU_FUEL_COLD_g.*GFC + SU_CO2_COLD_g.*CoC_g.*(1-Y_CAPT_g_t(HOUR,:)')) .* ...
                            (1-exp(X_PREV(:,K)./ TAU)));
                    else                         % "hot" alpha-beta approach
                        GEN_START_SHUT_COST(:,K) = (STATE_DIFF > 0) .* (ALPHA + SU_FUEL_COLD_g.*GFC + SU_CO2_COLD_g.*CoC_g.*(1-Y_CAPT_g_t(HOUR,:)') + BETA .* ...
                            (1-exp(X_PREV(:,K)./ TAU)));
                    end
                end
                %
                GEN_START_SHUT_COST(:,K) = GEN_START_SHUT_COST(:,K) + (STATE_DIFF  < 0) .* (GSDC + SD_FUEL_g.*GFC + SD_CO2_g.*CoC_g.*(1-Y_CAPT_g_t(HOUR,:)'));   % shut down cost added
                GEN_START_SHUT_COST(isnan(GEN_START_SHUT_COST))=0; % Neglecting this line was the source of mischeduling in previous versions
                %
                % find the generation [MW] and production cost for each unit
                [GEN_PRODUCTION,PROD_COST,~,~] = production(CURRENT_STATE,u_CAPT_g_t,PREVIOUS_STATE,GMIN,GMAX,D_t,W_t,W_on_t,W_of_t,S_t,DEMAND,RES_UP,RES_DN,HOUR,...
                    GNLC,GFC,GINC,m_A_g,m_B_g,m_C_g,c_A_g,c_B_g,c_C_g,NG,GRAMPUP,GRAMPDOWN,PRODUCTION_PREV,GEN_ORDER,COEF_A,COEF_B,COEF_C,DISPATCH_METHOD,ECO2_g,CoC_g,...
                    Y_MAX_CAPT_g,STK_g,VOM_g,Price,CC_FIXED_g,CC_OP_g,P_MIN_CAPT_g,P_MAX_CAPT_g,P_MIN_CAPT_g_t,P_MAX_CAPT_g_t,Y_CAPT_g_t,CC_VOM_g,CC_CSOLV_g,CC_SOLVD_g,...
                    CC_SOLVD_TH_g,CC_TRANS_g,GMIN_0,GMAX_0);
                COST_RAMP = ...
                    (STATE_DIFF == 0) .* ((GEN_PRODUCTION.*CURRENT_STATE - PRODUCTION_PREV(:,K).*PREVIOUS_STATE) > 0) ...% checks that there is a change in status and that it's a ramp upwards
                    .* abs(GEN_PRODUCTION.*CURRENT_STATE - PRODUCTION_PREV(:,K).*PREVIOUS_STATE).*RAMP_UP_g + (STATE_DIFF == 0) ...% calculates cost for ramp upwards
                    .* ((PRODUCTION_PREV(:,K).*PREVIOUS_STATE - GEN_PRODUCTION.*CURRENT_STATE)> 0) ...% checks if there is a change in status and that there is a ramp downwards
                    .* abs(PRODUCTION_PREV(:,K).*PREVIOUS_STATE - GEN_PRODUCTION.*CURRENT_STATE).*RAMP_DN_g; % calculates costs for ramp down
            end
            %
            X_TEMP(:,K) = X; % save the updated gen. work. times when moved from previous state to the current one
            PRODUCTION_TEMP(:,K) = GEN_PRODUCTION; % also save the updated gen. outputs
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
        end  % K
        %
        % among all transitions from each feasible state at previous hour
        % to the current state (at current hour), save up to MN with minimal total cost
        [MM,II] = sort(TOTAL_COST(TOTAL_COST ~= 0));
        for L = 1:MN                                          % L counts from 1 to N_PRED
            if isinf(MM(L)) == 0 % means that min UT and DT from the previous state to the current state are satisfied
                COUNTER = COUNTER +1;% counter counts how many possible transitions from any previous states to any current states are possible
                FC(COUNTER,1) = MM(L);
                TR(COUNTER,1:size(TR_PREV,2)) = TR_PREV(II(L),:);
                TR(COUNTER,end) = FEASIBLE_STATES_NUM(J);
                X_CURR(:,COUNTER) = X_TEMP(:,II(L));
                PRODUCTION_CURR(:,COUNTER) = PRODUCTION_TEMP(:,II(L));
            end % if isinf(MM(L))
        end % L
    end   % J
    %
    if COUNTER == 0
        if RELAXATION == 1
            HOUR_inf(HOUR)=1;
            relax_UT_DT_to_zero; % this reads the script on relaxing UT and DT
        else
            STR = ['NO FEASIBLE STATES FOR HOUR ',num2str(HOUR),'! PROGRAM TERMINATES!'];   % If the rest of the list is empty, then it means there are no feasible states,
            msgbox(STR,'NO FEASIBLE STATES','warn');                                        % and program terminates
            return
        end
    end
end   % HOUR
%============================================
% END OF SEARCHING PROCEDURE
%============================================
%% The search is complete. Now program finds the best solution (the least expensive state) at the last hour of the optimisation horizon.
[~,I]=min(FC(1:COUNTER));
BEST_PATH = TR(I,:).';    % find the best transition path and...
% ...evaluate the solution and generate output files and (print the results)
evaluate_solution(NT,BEST_PATH,LIST_STATES,GMIN,GMAX,D_t,W_t,W_on_t,W_of_t,W_c_t,S_t,S_c_t,W_t_Original,W_on_t_Original,W_of_t_Original,...
    S_t_Original,Tot_c_t,Net_D_t_Original,DEMAND,Price,RES_UP,RES_DN,GEN_ORDER,GNLC,GFC,GINC,m_A_g,m_B_g,m_C_g,c_A_g,c_B_g,c_C_g,GSC,SU_FUEL_COLD_g,SU_CO2_COLD_g,...
    GSDC,SD_FUEL_g,SD_CO2_g,INI_STATE,NG,GMINUP,GMINDOWN,...
    GRAMPUP,GRAMPDOWN,COEF_A,COEF_B,COEF_C,CC_FIXED_g,CC_OP_g,DISPATCH_METHOD,FLEXIBLE_CO2_CAPTURE_OPTION,GSTATINI,...
    GMINUP_Original,GMINDOWN_Original,START_UP_COST_METHOD,GCSTIME,GSH,ALPHA,BETA,TAU,ECO2_g,CoC_g,FLEX_g,Y_MAX_CAPT_g,...
    STK_g,VOM_g,P_MIN_CAPT_g,P_MAX_CAPT_g,Y_CAPT_g_t,u_CAPT_g_t,P_MIN_CAPT_g_t,...
    P_MAX_CAPT_g_t,CC_VOM_g,CC_CSOLV_g,CC_SOLVD_g,CC_SOLVD_TH_g,CC_TRANS_g,RAMP_UP_g,RAMP_DN_g,GMIN_0,GMAX_0,RELAXATION,EXP_COST_OPTION,STACKING_OPTION,RAMP_COSTS_INCLUDED,...
    On_subsidy,Off_subsidy,Uplift_c,Uplift_b,Uplift_VOLL,Uplift_SUM_G,Uplift_k,D_MIN_t,W_on_c_inert,W_of_c_inert,S_c_inert,S_c_inert_s,W_on_c_am_inert,W_of_c_am_inert,S_c_am_inert,...
    S_c_am_inert_s,W_on_c_feas,W_on_c_am_feas,W_of_c_feas,W_of_c_am_feas,S_c_feas,S_c_am_feas,S_c_feas_s,S_c_am_feas_s,HOUR_inf);
warning on
t=toc;
fprintf('\n Elapsed time: %15.4f min.\n\n',t/60);
end
dynamic_flex
%% Flexible Enumeration
function [GMAXlst,GMINlst,LIST_STATES,GEN_ORDER] = dynamic_flex(GNLC,GFC,GMAX,GMIN,GINC,NG,ECO2_g,CoC_g,Y_MAX_CAPT_g,STK_g,VOM_g,...
    CC_VOM_g,CC_TRANS_g,FLEX_g,STACKING_OPTION)
%
if STACKING_OPTION==2 % short-run marginal cost (SRMC)
    MERIT_COST       = 0*GNLC + GFC.*(1./GINC) + VOM_g + ECO2_g.*(1./GINC).*CoC_g.*(1-Y_MAX_CAPT_g) + ECO2_g.*Y_MAX_CAPT_g.*(1./GINC).*(CC_VOM_g + CC_TRANS_g) - STK_g;
else %average full-load cost
    MERIT_COST = (0*GNLC + GFC.*GMAX.*GINC/1000)./GMAX;
end
[M,GEN_ORDER] = sort(MERIT_COST);
%
%NG_FLEX must be less than NG_NON_FLEX or equal
NG_FLEX = nnz(FLEX_g); % number of nonzero matrix elements
NG_NON_FLEX = NG - NG_FLEX;
%
LIST_STATES_Flex = triu(ones(NG_FLEX));
LIST_STATES_Non_Flex = triu(ones(NG_NON_FLEX));
%
LIST_STATES_NON_FLEX_zeros = zeros(NG_FLEX,NG_NON_FLEX);
LIST_STATES_NON_FLEX_zeros = [LIST_STATES_Non_Flex; LIST_STATES_NON_FLEX_zeros];
LIST_STATES_FLEX_zeros_top = zeros(NG_NON_FLEX,NG_NON_FLEX);
LIST_STATES_FLEX_zeros_both = zeros(NG_FLEX,NG_NON_FLEX-NG_FLEX);
LIST_STATES_FLEX_zeros = [LIST_STATES_FLEX_zeros_both LIST_STATES_Flex];
LIST_STATES_FLEX_zeros = [LIST_STATES_FLEX_zeros_top; LIST_STATES_FLEX_zeros];
%
I = length(LIST_STATES_Non_Flex);
J = length(LIST_STATES_Flex);
%
LIST_STATES = (repmat(LIST_STATES_NON_FLEX_zeros(:,1),1,I) + LIST_STATES_FLEX_zeros);
%
for i = 2:I
    dataToAppend = (repmat(LIST_STATES_NON_FLEX_zeros(:,i),1,I) + LIST_STATES_FLEX_zeros);
    LIST_STATES = [LIST_STATES dataToAppend];
end
%
LIST_STATES = unique(LIST_STATES','rows');
LIST_STATES = LIST_STATES';
%
LIST_STATES(GEN_ORDER,:) = LIST_STATES(1:NG,:);
LIST_STATES = logical(LIST_STATES);
GMINlst = LIST_STATES.' * GMIN;
GMAXlst = LIST_STATES.' * GMAX;
%
[GMAXlst,INDEX]=sort(GMAXlst);
GMINlst = GMINlst(INDEX);
LIST_STATES = LIST_STATES(:,INDEX);
%
GMINlst(end:-1:1);
end