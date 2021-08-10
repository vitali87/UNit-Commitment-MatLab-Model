function [FEASIBLE_STATES_NUM,SUCCESS,DEMAND,W_on_t,W_of_t,S_t,W_c_t,S_c_t,RES_UP,RES_DN,W_on_c_feas,...
    W_on_c_am_feas,W_of_c_feas,W_of_c_am_feas,S_c_feas,S_c_am_feas,S_c_feas_s,S_c_am_feas_s] = find_feasible_states(GMINlst,GMAXlst,DEMAND,HOUR,RES_UP,RES_DN,...
    RES_UP_CCS,RES_DN_CCS,SEQUENTIAL_RES_CURTAILMENT,R_SP_UP_CAPT_t,R_SP_DN_CAPT_t,NT,W_on_t,W_of_t,...
    W_c_t,S_t,S_u,S_c_t,R_SP_UP_D_t,R_SP_DN_D_t,R_SP_UP_W_t,R_SP_DN_W_t,R_SP_UP_S_t,R_SP_DN_S_t,R_SP_UP_MAX,R_SP_DN_MAX,lambda_up,lambda_dn,S_t_large,S_t_small,W_u,D_t,...
    W_on_c_feas,W_on_c_am_feas,W_of_c_feas,W_of_c_am_feas,S_c_feas,S_c_am_feas,S_c_feas_s,S_c_am_feas_s)
% --------------------------------------------------------------------------------------------------------------
% Determines all feasible states from the list of possible states. Feasible
% states are the states where demand is between total min and total max
% output of commited generators. If no feasible states found, program
% prepares termination
% AMMENDMENT: If no feasible states are found, there
% is an option to sequentially curtail renewable electricity sources (RES)
% to land on feasible states
% OUTPUT:
% FEASIBLE_STATES_NUM   - vector of positions (columns) of feasible states in the list of states for current hour
% SUCCESS               - indicator: 1 - found at least one feasible states; 0 - no feasible states found
%----------------------------------------------------------------------------------------------------------------
if SEQUENTIAL_RES_CURTAILMENT==1
    Tot_c_feas1          = zeros(NT,1);
    W_c_t1               = zeros(NT,1);
    W_c_t11              = zeros(NT,1);
    S_c_t1               = zeros(NT,1);
    Tot_c_feas2          = zeros(NT,1);
    W_c_t2               = zeros(NT,1);
    W_c_t22              = zeros(NT,1);
    S_c_t2               = zeros(NT,1);
    Tot_c_feas3          = zeros(NT,1);
    W_c_t3               = zeros(NT,1);
    W_c_t33              = zeros(NT,1);
    S_c_t3               = zeros(NT,1);
    Tot_c_feas4          = zeros(NT,1);
    W_c_t4               = zeros(NT,1);
    W_c_t44              = zeros(NT,1);
    S_c_t4               = zeros(NT,1);
    Tot_c_feas5          = zeros(NT,1);
    W_c_t5               = zeros(NT,1);
    W_c_t55              = zeros(NT,1);
    S_c_t5               = zeros(NT,1);
    W_on_c_feas1         = zeros(NT,1);
    W_on_c_am_feas1      = zeros(NT,1);
    W_of_c_feas1         = zeros(NT,1);
    W_of_c_am_feas1      = zeros(NT,1);
    S_c_feas1            = zeros(NT,1);
    S_c_feas1_s          = zeros(NT,1);
    S_c_am_feas1         = zeros(NT,1);
    S_c_am_feas1_s       = zeros(NT,1);
    W_on_c_feas2         = zeros(NT,1);
    W_on_c_am_feas2      = zeros(NT,1);
    W_of_c_feas2         = zeros(NT,1);
    W_of_c_am_feas2      = zeros(NT,1);
    S_c_feas2            = zeros(NT,1);
    S_c_feas2_s          = zeros(NT,1);
    S_c_am_feas2         = zeros(NT,1);
    S_c_am_feas2_s       = zeros(NT,1);
    W_on_c_feas3         = zeros(NT,1);
    W_on_c_am_feas3      = zeros(NT,1);
    W_of_c_feas3         = zeros(NT,1);
    W_of_c_am_feas3      = zeros(NT,1);
    S_c_feas3            = zeros(NT,1);
    S_c_feas3_s          = zeros(NT,1);
    S_c_am_feas3         = zeros(NT,1);
    S_c_am_feas3_s       = zeros(NT,1);
    W_on_c_feas4         = zeros(NT,1);
    W_on_c_am_feas4      = zeros(NT,1);
    W_of_c_feas4         = zeros(NT,1);
    W_of_c_am_feas4      = zeros(NT,1);
    S_c_feas4            = zeros(NT,1);
    S_c_feas4_s          = zeros(NT,1);
    S_c_am_feas4         = zeros(NT,1);
    S_c_am_feas4_s       = zeros(NT,1);
    W_on_c_feas5         = zeros(NT,1);
    W_on_c_am_feas5      = zeros(NT,1);
    W_of_c_feas5         = zeros(NT,1);
    W_of_c_am_feas5      = zeros(NT,1);
    S_c_feas5            = zeros(NT,1);
    S_c_feas5_s          = zeros(NT,1);
    S_c_am_feas5         = zeros(NT,1);
    S_c_am_feas5_s       = zeros(NT,1);
    W_on_t_after_inert   = W_on_t;
    W_of_t_after_inert   = W_of_t;
    S_t_after_inert      = S_t_large;
    S_t_after_inert_s    = S_t_small;
    %
    FEASIBLE_STATES_NUM = find((GMINlst <= DEMAND(HOUR) - RES_DN(HOUR) - RES_DN_CCS(HOUR)) & (DEMAND(HOUR) + RES_UP(HOUR) + RES_UP_CCS(HOUR) <= GMAXlst),1);
    if isempty(FEASIBLE_STATES_NUM)
%         R_SP_UP_W_t(HOUR) = R_SP_UP_W_t(HOUR)*0.5;
%         R_SP_DN_W_t(HOUR) = R_SP_DN_W_t(HOUR)*0.5;
        R_SP_UP_W_t(HOUR) = R_SP_UP_W_t(HOUR); % Consider all the upward reserve needed for wind
        R_SP_DN_W_t(HOUR) = R_SP_DN_W_t(HOUR); % Consider all the downward reserve needed for wind
        R_SP_UP_W_t(R_SP_UP_W_t<0) = 0;
%         R_SP_UP_S_t(HOUR) = R_SP_UP_S_t(HOUR)*0.5; % Solar reserve requirements are reduced to half (LV)
%         R_SP_DN_S_t(HOUR) = R_SP_DN_S_t(HOUR)*0.5; % (LV)
        R_SP_UP_S_t(HOUR) = R_SP_UP_S_t(HOUR); % Consider all the upward reserve needed for solar
        R_SP_DN_S_t(HOUR) = R_SP_DN_S_t(HOUR); % Consider all the downward reserve needed for solar
        R_SP_UP_S_t(R_SP_UP_S_t<0) = 0; % In case solar reserve becomes negative (LV)
        R_SP_UP_D_W_S_t(HOUR) = ((R_SP_UP_D_t(HOUR)).^2 + (R_SP_UP_W_t(HOUR)).^2 + (R_SP_UP_S_t(HOUR)).^2).^0.5; % Solar added LV
        R_SP_DN_D_W_S_t(HOUR) = ((R_SP_DN_D_t(HOUR)).^2 + (R_SP_DN_W_t(HOUR)).^2 + (R_SP_DN_S_t(HOUR)).^2).^0.5; % Solar added LV
        RES_UP(HOUR) = R_SP_UP_MAX(HOUR) + lambda_up*R_SP_UP_D_W_S_t(HOUR) - R_SP_UP_CAPT_t(HOUR);
        RES_UP(RES_UP<0) = 0;
        RES_DN(HOUR) = R_SP_DN_MAX(HOUR) + lambda_dn*R_SP_DN_D_W_S_t(HOUR) - R_SP_DN_CAPT_t(HOUR);
        RES_DN(RES_DN<0) = 0;
    else
        SUCCESS = 1;
    end
    %
    FEASIBLE_STATES_NUM = find((GMINlst <= DEMAND(HOUR) - RES_DN(HOUR) - RES_DN_CCS(HOUR)) & (DEMAND(HOUR) + RES_UP(HOUR) + RES_UP_CCS(HOUR) <= GMAXlst), 1);
    if isempty(FEASIBLE_STATES_NUM) % Large solar PV is curtailed first if infeasible (LV)
        Tot_c_feas1(HOUR) = 1000;
        S_t_large(HOUR) = S_t_large(HOUR) - min(S_t_large(HOUR),Tot_c_feas1(HOUR)); % (LV)
        S_c_feas1(HOUR)=(S_t_after_inert(HOUR)>S_t_large(HOUR)); % count how many times large solar is curtailed for feasibility requirement 1
        S_c_am_feas1(HOUR) = (S_t_after_inert(HOUR)>S_t_large(HOUR)).*((S_t_after_inert(HOUR)-S_t_large(HOUR)));% how much large solar is curtailed for feasibility requirement 1
        if S_t_large (HOUR)==0
            W_c_t1(HOUR) = Tot_c_feas1(HOUR) - S_t_after_inert(HOUR);
            W_on_t(HOUR) = W_on_t(HOUR) - min(W_c_t1(HOUR),W_on_t(HOUR));
            W_on_c_feas1(HOUR)=(W_on_t_after_inert(HOUR)>W_on_t(HOUR)); % count how many times offshore wind is curtailed for feasibility requirement 1
            W_on_c_am_feas1(HOUR) = (W_on_t_after_inert(HOUR)>W_on_t(HOUR)).*((W_on_t_after_inert(HOUR)-W_on_t(HOUR)));% how much onshore wind is curtailed for feasibility requirement 1
            if W_on_t(HOUR) == 0
                W_c_t11(HOUR) = W_c_t1(HOUR) - W_on_t_after_inert(HOUR);
                W_of_t(HOUR) = W_of_t(HOUR) - min(W_of_t(HOUR),W_c_t11(HOUR));
                W_of_c_feas1(HOUR)=(W_of_t_after_inert(HOUR)>W_of_t(HOUR)); % count how many times offshore wind is curtailed for feasibility requirement 1
                W_of_c_am_feas1(HOUR) = (W_of_t_after_inert(HOUR)>W_of_t(HOUR)).*((W_of_t_after_inert(HOUR)-W_of_t(HOUR)));% how much onshore wind is curtailed for feasibility requirement 1
                if W_of_t(HOUR) == 0 % Small scale solar curtailment in case large solar and wind are not enough to curtail the 1 GW
                    S_c_t1(HOUR) = Tot_c_feas1(HOUR) - S_t_after_inert(HOUR) - W_on_t_after_inert(HOUR) - W_of_t_after_inert(HOUR);
                    S_t_small(HOUR) = S_t_small(HOUR) - min(S_t_small(HOUR),S_c_t1(HOUR));
                    S_c_feas1_s(HOUR)=(S_t_after_inert_s(HOUR)>S_t_small(HOUR)); % count how many times solar is curtailed for feasibility requirement 1
                    S_c_am_feas1_s(HOUR) = (S_t_after_inert_s(HOUR)>S_t_small(HOUR)).*((S_t_after_inert_s(HOUR)-S_t_small(HOUR)));% how much solar is curtailed for feasibility requirement 1
                end
            end
        end
        W_c_t(HOUR) = W_c_t(HOUR) + W_c_t1(HOUR) - S_c_t1(HOUR); % Wind curtailed initially (min load) plus wind curtailed for feasibility
        W_t(HOUR) = W_on_t(HOUR) + W_of_t(HOUR);
        S_t(HOUR) = S_t_large(HOUR) + S_t_small(HOUR);
        S_c_t(HOUR) = S_c_t(HOUR) + Tot_c_feas1(HOUR) - W_c_t1(HOUR) + S_c_t1(HOUR); % Solar curtailed initially (min load) plus solar curtailed for feasibility
        R_SP_UP_S_t(HOUR) = S_u*(S_t(HOUR)); % S_t changes so we need to define the reserve requirements again (LV) Why is this not consistent with wind? (VA)
%         R_SP_UP_S_t(HOUR) = min(S_u*S_t(HOUR),max(0,S_t(HOUR) - S_c_t(HOUR))); % To define solar reserves the same way as wind (it does not make a difference)
        R_SP_DN_S_t(HOUR) = S_u*S_t(HOUR).*0.5; % (LV) I did not change this back to the entire reserve becasue the difference is negligible
        R_SP_UP_S_t(R_SP_UP_S_t<0) = 0; % (LV)
        DEMAND(HOUR) = D_t(HOUR) - W_t(HOUR) - S_t(HOUR); % Modified (LV)
        R_SP_UP_W_t(HOUR) = min(W_u*(W_t(HOUR)),W_t(HOUR) - W_c_t(HOUR));
%         R_SP_DN_W_t(HOUR) = W_u*(W_t(HOUR)).*0.5; % 
        R_SP_DN_W_t(HOUR) = W_u*(W_t(HOUR)); %  Consider all the upward reserve needed for wind
        R_SP_UP_W_t(R_SP_UP_W_t<0) = 0;
        R_SP_UP_D_W_S_t(HOUR) = ((R_SP_UP_D_t(HOUR)).^2 + (R_SP_UP_W_t(HOUR)).^2 + (R_SP_UP_S_t(HOUR)).^2).^0.5; % Solar added LV
        R_SP_DN_D_W_S_t(HOUR) = ((R_SP_DN_D_t(HOUR)).^2 + (R_SP_DN_W_t(HOUR)).^2 + (R_SP_DN_S_t(HOUR)).^2).^0.5; % Solar added LV
        RES_UP(HOUR) = R_SP_UP_MAX(HOUR) + lambda_up*R_SP_UP_D_W_S_t(HOUR) - R_SP_UP_CAPT_t(HOUR);
        RES_UP(RES_UP<0) = 0;
        RES_DN(HOUR) = R_SP_DN_MAX(HOUR) + lambda_dn*R_SP_DN_D_W_S_t(HOUR) - R_SP_DN_CAPT_t(HOUR);
        RES_DN(RES_DN<0) = 0;
        W_on_t_after_feas1=W_on_t;
        W_of_t_after_feas1=W_of_t;
        S_t_after_feas1=S_t_large;
        S_t_after_feas1_s=S_t_small;
    else
        SUCCESS = 1;
    end
    %
    FEASIBLE_STATES_NUM = find((GMINlst <= DEMAND(HOUR) - RES_DN(HOUR) - RES_DN_CCS(HOUR)) & (DEMAND(HOUR) + RES_UP(HOUR) + RES_UP_CCS(HOUR) <= GMAXlst), 1);
    if isempty(FEASIBLE_STATES_NUM)
        Tot_c_feas2(HOUR) = 1000;
        S_t_large(HOUR)= S_t_large(HOUR) - min(S_t_large(HOUR),Tot_c_feas2(HOUR)); % (LV)
        S_c_feas2(HOUR)=(S_t_after_feas1(HOUR)>S_t_large(HOUR)); % count how many times large solar is curtailed for feasibility requirement 2
        S_c_am_feas2(HOUR) = (S_t_after_feas1(HOUR)>S_t_large(HOUR)).*((S_t_after_feas1(HOUR)-S_t_large(HOUR)));% how much large solar is curtailed for feasibility requirement 2
        if S_t_large (HOUR)== 0
            W_c_t2(HOUR) = Tot_c_feas2(HOUR) - S_t_after_feas1(HOUR);
            W_on_t(HOUR) = W_on_t(HOUR) - min(W_c_t2(HOUR),W_on_t(HOUR));
            W_on_c_feas2(HOUR)=(W_on_t_after_feas1(HOUR)>W_on_t(HOUR)); % count how many times offshore wind is curtailed for feasibility requirement 2
            W_on_c_am_feas2(HOUR) = (W_on_t_after_feas1(HOUR)>W_on_t(HOUR)).*((W_on_t_after_feas1(HOUR)-W_on_t(HOUR)));% how much onshore wind is curtailed for feasibility requirement 2
            if W_on_t(HOUR) == 0
                W_c_t22(HOUR) = W_c_t2(HOUR) - W_on_t_after_feas1(HOUR);
                W_of_t(HOUR) = W_of_t(HOUR) - min(W_of_t(HOUR),W_c_t22(HOUR));
                W_of_c_feas2(HOUR)=(W_of_t_after_feas1(HOUR)>W_of_t(HOUR)); % count how many times offshore wind is curtailed for feasibility requirement 2
                W_of_c_am_feas2(HOUR) = (W_of_t_after_feas1(HOUR)>W_of_t(HOUR)).*((W_of_t_after_feas1(HOUR)-W_of_t(HOUR)));% how much onshore wind is curtailed for feasibility requirement 2
                if W_of_t(HOUR) == 0 % Small scale solar curtailment in case large solar and wind are not enough to curtail the 1 GW
                    S_c_t2(HOUR) = Tot_c_feas2(HOUR) - S_t_after_feas1(HOUR) - W_on_t_after_feas1(HOUR) - W_of_t_after_feas1(HOUR);
                    S_t_small(HOUR) = S_t_small(HOUR) - min(S_t_small(HOUR),S_c_t2(HOUR));
                    S_c_feas2_s(HOUR)=(S_t_after_feas1_s(HOUR)>S_t_small(HOUR)); % count how many times solar is curtailed for feasibility requirement 2
                    S_c_am_feas2_s(HOUR) = (S_t_after_feas1_s(HOUR)>S_t_small(HOUR)).*((S_t_after_feas1_s(HOUR)-S_t_small(HOUR)));% how much solar is curtailed for feasibility requirement 2
                end
            end
        end
        W_c_t(HOUR) = W_c_t(HOUR) + W_c_t2(HOUR) - S_c_t2(HOUR); % Wind curtailed initially (min load) plus wind curtailed for feasibility
        W_t(HOUR) = W_on_t(HOUR) + W_of_t(HOUR);
        S_t(HOUR) = S_t_large(HOUR) + S_t_small(HOUR);
        S_c_t(HOUR) = S_c_t(HOUR) + Tot_c_feas2(HOUR) - W_c_t2(HOUR) + S_c_t2(HOUR); % Solar curtailed initially (min load) plus solar curtailed for feasibility
        R_SP_UP_S_t(HOUR) = S_u*(S_t(HOUR));
%         R_SP_UP_S_t(HOUR) = min(S_u*S_t(HOUR),max(0,S_t(HOUR) - S_c_t(HOUR))); 
        R_SP_DN_S_t(HOUR) = S_u*S_t(HOUR).*0.5; % (LV)
        R_SP_UP_S_t(R_SP_UP_S_t<0) = 0; % (LV)
        DEMAND(HOUR) = D_t(HOUR) - W_t(HOUR)- S_t(HOUR); % Modified (LV)
        R_SP_UP_W_t(HOUR) = min(W_u*(W_t(HOUR)),W_t(HOUR) - W_c_t(HOUR));
%         R_SP_DN_W_t(HOUR) = W_u*(W_t(HOUR)).*0.5;        
        R_SP_DN_W_t(HOUR) = W_u*(W_t(HOUR));%  Consider all the upward reserve needed for wind
        R_SP_UP_W_t(R_SP_UP_W_t<0) = 0;
        R_SP_UP_D_W_S_t(HOUR) = ((R_SP_UP_D_t(HOUR)).^2 + (R_SP_UP_W_t(HOUR)).^2 + (R_SP_UP_S_t(HOUR)).^2).^0.5; % Solar added LV
        R_SP_DN_D_W_S_t(HOUR) = ((R_SP_DN_D_t(HOUR)).^2 + (R_SP_DN_W_t(HOUR)).^2 + (R_SP_DN_S_t(HOUR)).^2).^0.5; % Solar added LV
        RES_UP(HOUR) = R_SP_UP_MAX(HOUR) + lambda_up*R_SP_UP_D_W_S_t(HOUR) - R_SP_UP_CAPT_t(HOUR);
        RES_UP(RES_UP<0) = 0;
        RES_DN(HOUR) = R_SP_DN_MAX(HOUR) + lambda_dn*R_SP_DN_D_W_S_t(HOUR) - R_SP_DN_CAPT_t(HOUR);
        RES_DN(RES_DN<0) = 0;
        W_on_t_after_feas2=W_on_t;
        W_of_t_after_feas2=W_of_t;
        S_t_after_feas2=S_t_large;
        S_t_after_feas2_s=S_t_small;
    else
        SUCCESS = 1;
    end
    %
    FEASIBLE_STATES_NUM = find((GMINlst <= DEMAND(HOUR) - RES_DN(HOUR) - RES_DN_CCS(HOUR)) & (DEMAND(HOUR) + RES_UP(HOUR) + RES_UP_CCS(HOUR) <= GMAXlst), 1);
    if isempty(FEASIBLE_STATES_NUM)
        Tot_c_feas3(HOUR) = 1000;
        S_t_large(HOUR)= S_t_large(HOUR) - min(S_t_large(HOUR),Tot_c_feas3(HOUR)); % (LV)
        S_c_feas3(HOUR)=(S_t_after_feas2(HOUR)>S_t_large(HOUR)); % count how many times large solar is curtailed for feasibility requirement 2
        S_c_am_feas3(HOUR) = (S_t_after_feas2(HOUR)>S_t_large(HOUR)).*((S_t_after_feas2(HOUR)-S_t_large(HOUR)));% how much large solar is curtailed for feasibility requirement 2
        if S_t_large (HOUR)== 0
            W_c_t3(HOUR) = Tot_c_feas3(HOUR) - S_t_after_feas2(HOUR);
            W_on_t(HOUR) = W_on_t(HOUR) - min(W_c_t3(HOUR),W_on_t(HOUR));
            W_on_c_feas3(HOUR)=(W_on_t_after_feas2(HOUR)>W_on_t(HOUR)); % count how many times offshore wind is curtailed for feasibility requirement 2
            W_on_c_am_feas3(HOUR) = (W_on_t_after_feas2(HOUR)>W_on_t(HOUR)).*((W_on_t_after_feas2(HOUR)-W_on_t(HOUR)));% how much onshore wind is curtailed for feasibility requirement 2
            if W_on_t(HOUR) == 0
                W_c_t33(HOUR) = W_c_t3(HOUR) - W_on_t_after_feas2(HOUR);
                W_of_t(HOUR) = W_of_t(HOUR) - min(W_of_t(HOUR),W_c_t33(HOUR));
                W_of_c_feas3(HOUR)=(W_of_t_after_feas2(HOUR)>W_of_t(HOUR)); % count how many times offshore wind is curtailed for feasibility requirement 2
                W_of_c_am_feas3(HOUR) = (W_of_t_after_feas2(HOUR)>W_of_t(HOUR)).*((W_of_t_after_feas2(HOUR)-W_of_t(HOUR)));% how much onshore wind is curtailed for feasibility requirement 2
                if W_of_t(HOUR) == 0 % Small scale solar curtailment in case large solar and wind are not enough to curtail the 1 GW
                    S_c_t3(HOUR) = Tot_c_feas3(HOUR) - S_t_after_feas2(HOUR) - W_on_t_after_feas2(HOUR) - W_of_t_after_feas2(HOUR);
                    S_t_small(HOUR) = S_t_small(HOUR) - min(S_t_small(HOUR),S_c_t3(HOUR));
                    S_c_feas3_s(HOUR)=(S_t_after_feas2_s(HOUR)>S_t_small(HOUR)); % count how many times solar is curtailed for feasibility requirement 2
                    S_c_am_feas3_s(HOUR) = (S_t_after_feas2_s(HOUR)>S_t_small(HOUR)).*((S_t_after_feas2_s(HOUR)-S_t_small(HOUR)));% how much solar is curtailed for feasibility requirement 2
                end
            end
        end
        W_c_t(HOUR) = W_c_t(HOUR) + W_c_t3(HOUR) - S_c_t3(HOUR); % Wind curtailed initially (min load) plus wind curtailed for feasibility
        W_t(HOUR) = W_on_t(HOUR) + W_of_t(HOUR);
        S_t(HOUR) = S_t_large(HOUR) + S_t_small(HOUR);
        S_c_t(HOUR) = S_c_t(HOUR) + Tot_c_feas3(HOUR) - W_c_t3(HOUR) + S_c_t3(HOUR); % Solar curtailed initially (min load) plus solar curtailed for feasibility
        R_SP_UP_S_t(HOUR) = S_u*(S_t(HOUR));
%         R_SP_UP_S_t(HOUR) = min(S_u*S_t(HOUR),max(0,S_t(HOUR) - S_c_t(HOUR))); 
        R_SP_DN_S_t(HOUR) = S_u*S_t(HOUR).*0.5; % (LV)
        R_SP_UP_S_t(R_SP_UP_S_t<0) = 0; % (LV)
        DEMAND(HOUR) = D_t(HOUR) - W_t(HOUR)- S_t(HOUR); % Modified (LV)
        R_SP_UP_W_t(HOUR) = min(W_u*(W_t(HOUR)),W_t(HOUR) - W_c_t(HOUR));
%         R_SP_DN_W_t(HOUR) = W_u*(W_t(HOUR)).*0.5;
        R_SP_DN_W_t(HOUR) = W_u*(W_t(HOUR));%  Consider all the upward reserve needed for wind
        R_SP_UP_W_t(R_SP_UP_W_t<0) = 0;
        R_SP_UP_D_W_S_t(HOUR) = ((R_SP_UP_D_t(HOUR)).^2 + (R_SP_UP_W_t(HOUR)).^2 + (R_SP_UP_S_t(HOUR)).^2).^0.5; % Solar added LV
        R_SP_DN_D_W_S_t(HOUR) = ((R_SP_DN_D_t(HOUR)).^2 + (R_SP_DN_W_t(HOUR)).^2 + (R_SP_DN_S_t(HOUR)).^2).^0.5; % Solar added LV
        RES_UP(HOUR) = R_SP_UP_MAX(HOUR) + lambda_up*R_SP_UP_D_W_S_t(HOUR) - R_SP_UP_CAPT_t(HOUR);
        RES_UP(RES_UP<0) = 0;
        RES_DN(HOUR) = R_SP_DN_MAX(HOUR) + lambda_dn*R_SP_DN_D_W_S_t(HOUR) - R_SP_DN_CAPT_t(HOUR);
        RES_DN(RES_DN<0) = 0;
        W_on_t_after_feas3=W_on_t;
        W_of_t_after_feas3=W_of_t;
        S_t_after_feas3=S_t_large;
        S_t_after_feas3_s=S_t_small;
    else
        SUCCESS = 1;
    end
    %
    FEASIBLE_STATES_NUM = find((GMINlst <= DEMAND(HOUR) - RES_DN(HOUR) - RES_DN_CCS(HOUR)) & (DEMAND(HOUR) + RES_UP(HOUR) + RES_UP_CCS(HOUR) <= GMAXlst), 1);
    if isempty(FEASIBLE_STATES_NUM)
        Tot_c_feas4(HOUR) = 2000;
        S_t_large(HOUR)= S_t_large(HOUR) - min(S_t_large(HOUR),Tot_c_feas4(HOUR)); % (LV)
        S_c_feas4(HOUR)=(S_t_after_feas3(HOUR)>S_t_large(HOUR)); % count how many times large solar is curtailed for feasibility requirement 2
        S_c_am_feas4(HOUR) = (S_t_after_feas3(HOUR)>S_t_large(HOUR)).*((S_t_after_feas3(HOUR)-S_t_large(HOUR)));% how much large solar is curtailed for feasibility requirement 2
        if S_t_large (HOUR)== 0
            W_c_t4(HOUR) = Tot_c_feas4(HOUR) - S_t_after_feas3(HOUR);
            W_on_t(HOUR) = W_on_t(HOUR) - min(W_c_t4(HOUR),W_on_t(HOUR));
            W_on_c_feas4(HOUR)=(W_on_t_after_feas3(HOUR)>W_on_t(HOUR)); % count how many times offshore wind is curtailed for feasibility requirement 2
            W_on_c_am_feas4(HOUR) = (W_on_t_after_feas3(HOUR)>W_on_t(HOUR)).*((W_on_t_after_feas3(HOUR)-W_on_t(HOUR)));% how much onshore wind is curtailed for feasibility requirement 2
            if W_on_t(HOUR) == 0
                W_c_t44(HOUR) = W_c_t4(HOUR) - W_on_t_after_feas3(HOUR);
                W_of_t(HOUR) = W_of_t(HOUR) - min(W_of_t(HOUR),W_c_t44(HOUR));
                W_of_c_feas4(HOUR)=(W_of_t_after_feas3(HOUR)>W_of_t(HOUR)); % count how many times offshore wind is curtailed for feasibility requirement 2
                W_of_c_am_feas4(HOUR) = (W_of_t_after_feas3(HOUR)>W_of_t(HOUR)).*((W_of_t_after_feas3(HOUR)-W_of_t(HOUR)));% how much onshore wind is curtailed for feasibility requirement 2
                if W_of_t(HOUR) == 0 % Small scale solar curtailment in case large solar and wind are not enough to curtail the 1 GW
                    S_c_t4(HOUR) = Tot_c_feas4(HOUR) - S_t_after_feas3(HOUR) - W_on_t_after_feas3(HOUR) - W_of_t_after_feas3(HOUR);
                    S_t_small(HOUR) = S_t_small(HOUR) - min(S_t_small(HOUR),S_c_t4(HOUR));
                    S_c_feas4_s(HOUR)=(S_t_after_feas3_s(HOUR)>S_t_small(HOUR)); % count how many times solar is curtailed for feasibility requirement 2
                    S_c_am_feas4_s(HOUR) = (S_t_after_feas3_s(HOUR)>S_t_small(HOUR)).*((S_t_after_feas3_s(HOUR)-S_t_small(HOUR)));% how much solar is curtailed for feasibility requirement 2
                end
            end
        end
        W_c_t(HOUR) = W_c_t(HOUR) + W_c_t4(HOUR) - S_c_t4(HOUR); % Wind curtailed initially (min load) plus wind curtailed for feasibility
        W_t(HOUR) = W_on_t(HOUR) + W_of_t(HOUR);
        S_t(HOUR) = S_t_large(HOUR) + S_t_small(HOUR);
        S_c_t(HOUR) = S_c_t(HOUR) + Tot_c_feas4(HOUR) - W_c_t4(HOUR) + S_c_t4(HOUR); % Solar curtailed initially (min load) plus solar curtailed for feasibility
        R_SP_UP_S_t(HOUR) = S_u*(S_t(HOUR));
%         R_SP_UP_S_t(HOUR) = min(S_u*S_t(HOUR),max(0,S_t(HOUR) - S_c_t(HOUR))); 
        R_SP_DN_S_t(HOUR) = S_u*S_t(HOUR).*0.5; % (LV)
        R_SP_UP_S_t(R_SP_UP_S_t<0) = 0; % (LV)
        DEMAND(HOUR) = D_t(HOUR) - W_t(HOUR)- S_t(HOUR); % Modified (LV)
        R_SP_UP_W_t(HOUR) = min(W_u*(W_t(HOUR)),W_t(HOUR) - W_c_t(HOUR));
%         R_SP_DN_W_t(HOUR) = W_u*(W_t(HOUR)).*0.5;          
        R_SP_DN_W_t(HOUR) = W_u*(W_t(HOUR));%  Consider all the upward reserve needed for wind
        R_SP_UP_W_t(R_SP_UP_W_t<0) = 0;
        R_SP_UP_D_W_S_t(HOUR) = ((R_SP_UP_D_t(HOUR)).^2 + (R_SP_UP_W_t(HOUR)).^2 + (R_SP_UP_S_t(HOUR)).^2).^0.5; % Solar added LV
        R_SP_DN_D_W_S_t(HOUR) = ((R_SP_DN_D_t(HOUR)).^2 + (R_SP_DN_W_t(HOUR)).^2 + (R_SP_DN_S_t(HOUR)).^2).^0.5; % Solar added LV
        RES_UP(HOUR) = R_SP_UP_MAX(HOUR) + lambda_up*R_SP_UP_D_W_S_t(HOUR) - R_SP_UP_CAPT_t(HOUR);
        RES_UP(RES_UP<0) = 0;
        RES_DN(HOUR) = R_SP_DN_MAX(HOUR) + lambda_dn*R_SP_DN_D_W_S_t(HOUR) - R_SP_DN_CAPT_t(HOUR);
        RES_DN(RES_DN<0) = 0;
        W_on_t_after_feas4=W_on_t;
        W_of_t_after_feas4=W_of_t;
        S_t_after_feas4=S_t_large;
        S_t_after_feas4_s=S_t_small;
    else
        SUCCESS = 1;
    end
    %
    FEASIBLE_STATES_NUM = find((GMINlst <= DEMAND(HOUR) - RES_DN(HOUR) - RES_DN_CCS(HOUR)) & (DEMAND(HOUR) + RES_UP(HOUR) + RES_UP_CCS(HOUR) <= GMAXlst), 1);
    if isempty(FEASIBLE_STATES_NUM)
        Tot_c_feas5(HOUR) = 2000;
        S_t_large(HOUR)= S_t_large(HOUR) - min(S_t_large(HOUR),Tot_c_feas5(HOUR)); % (LV)
        S_c_feas5(HOUR)=(S_t_after_feas4(HOUR)>S_t_large(HOUR)); % count how many times large solar is curtailed for feasibility requirement 2
        S_c_am_feas5(HOUR) = (S_t_after_feas4(HOUR)>S_t_large(HOUR)).*((S_t_after_feas4(HOUR)-S_t_large(HOUR)));% how much large solar is curtailed for feasibility requirement 2
        if S_t_large (HOUR)== 0
            W_c_t5(HOUR) = Tot_c_feas5(HOUR) - S_t_after_feas4(HOUR);
            W_on_t(HOUR) = W_on_t(HOUR) - min(W_c_t5(HOUR),W_on_t(HOUR));
            W_on_c_feas5(HOUR)=(W_on_t_after_feas4(HOUR)>W_on_t(HOUR)); % count how many times offshore wind is curtailed for feasibility requirement 2
            W_on_c_am_feas5(HOUR) = (W_on_t_after_feas4(HOUR)>W_on_t(HOUR)).*((W_on_t_after_feas4(HOUR)-W_on_t(HOUR)));% how much onshore wind is curtailed for feasibility requirement 2
            if W_on_t(HOUR) == 0
                W_c_t55(HOUR) = W_c_t5(HOUR) - W_on_t_after_feas4(HOUR);
                W_of_t(HOUR) = W_of_t(HOUR) - min(W_of_t(HOUR),W_c_t55(HOUR));
                W_of_c_feas5(HOUR)=(W_of_t_after_feas4(HOUR)>W_of_t(HOUR)); % count how many times offshore wind is curtailed for feasibility requirement 2
                W_of_c_am_feas5(HOUR) = (W_of_t_after_feas4(HOUR)>W_of_t(HOUR)).*((W_of_t_after_feas4(HOUR)-W_of_t(HOUR)));% how much onshore wind is curtailed for feasibility requirement 2
                if W_of_t(HOUR) == 0 % Small scale solar curtailment in case large solar and wind are not enough to curtail the 1 GW
                    S_c_t5(HOUR) = Tot_c_feas5(HOUR) - S_t_after_feas4(HOUR) - W_on_t_after_feas4(HOUR) - W_of_t_after_feas4(HOUR);
                    S_t_small(HOUR) = S_t_small(HOUR) - min(S_t_small(HOUR),S_c_t5(HOUR));
                    S_c_feas5_s(HOUR)=(S_t_after_feas4_s(HOUR)>S_t_small(HOUR)); % count how many times solar is curtailed for feasibility requirement 2
                    S_c_am_feas5_s(HOUR) = (S_t_after_feas4_s(HOUR)>S_t_small(HOUR)).*((S_t_after_feas4_s(HOUR)-S_t_small(HOUR)));% how much solar is curtailed for feasibility requirement 2
                end
            end
        end
        W_c_t(HOUR) = W_c_t(HOUR) + W_c_t5(HOUR) - S_c_t5(HOUR); % Wind curtailed initially (min load) plus wind curtailed for feasibility
        W_t(HOUR) = W_on_t(HOUR) + W_of_t(HOUR);
        S_t(HOUR) = S_t_large(HOUR) + S_t_small(HOUR);
        S_c_t(HOUR) = S_c_t(HOUR) + Tot_c_feas5(HOUR) - W_c_t5(HOUR) + S_c_t5(HOUR); % Solar curtailed initially (min load) plus solar curtailed for feasibility
        R_SP_UP_S_t(HOUR) = S_u*(S_t(HOUR));
%         R_SP_UP_S_t(HOUR) = min(S_u*S_t(HOUR),max(0,S_t(HOUR) - S_c_t(HOUR))); 
        R_SP_DN_S_t(HOUR) = S_u*S_t(HOUR).*0.5; % (LV)
        R_SP_UP_S_t(R_SP_UP_S_t<0) = 0; % (LV)
        DEMAND(HOUR) = D_t(HOUR) - W_t(HOUR)- S_t(HOUR); % Modified (LV)
        R_SP_UP_W_t(HOUR) = min(W_u*(W_t(HOUR)),W_t(HOUR) - W_c_t(HOUR));
%         R_SP_DN_W_t(HOUR) = W_u*(W_t(HOUR)).*0.5; 
        R_SP_DN_W_t(HOUR) = W_u*(W_t(HOUR));%  Consider all the upward reserve needed for wind
        R_SP_UP_W_t(R_SP_UP_W_t<0) = 0;
        R_SP_UP_D_W_S_t(HOUR) = ((R_SP_UP_D_t(HOUR)).^2 + (R_SP_UP_W_t(HOUR)).^2 + (R_SP_UP_S_t(HOUR)).^2).^0.5; % Solar added LV
        R_SP_DN_D_W_S_t(HOUR) = ((R_SP_DN_D_t(HOUR)).^2 + (R_SP_DN_W_t(HOUR)).^2 + (R_SP_DN_S_t(HOUR)).^2).^0.5; % Solar added LV
        RES_UP(HOUR) = R_SP_UP_MAX(HOUR) + lambda_up*R_SP_UP_D_W_S_t(HOUR) - R_SP_UP_CAPT_t(HOUR);
        RES_UP(RES_UP<0) = 0;
        RES_DN(HOUR) = R_SP_DN_MAX(HOUR) + lambda_dn*R_SP_DN_D_W_S_t(HOUR) - R_SP_DN_CAPT_t(HOUR);
        RES_DN(RES_DN<0) = 0;
    else
        SUCCESS = 1;
    end
    %
    FEASIBLE_STATES_NUM = find((GMINlst <= DEMAND(HOUR) - RES_DN(HOUR) - RES_DN_CCS(HOUR)) & (DEMAND(HOUR) + RES_UP(HOUR) + RES_UP_CCS(HOUR) <= GMAXlst));
    if isempty(FEASIBLE_STATES_NUM)
        SUCCESS = 0;
        STR = ['NO FEASIBLE STATES FOR HOUR',num2str(HOUR), ' PROGRAM TERMINATES!'];
        msgbox(STR,'NO FEASIBLE STATES','warn');
        return
    else
        SUCCESS = 1;
        W_on_c_feas(HOUR) = max([W_on_c_feas1(HOUR) W_on_c_feas2(HOUR) W_on_c_feas3(HOUR) W_on_c_feas4(HOUR) W_on_c_feas5(HOUR)]);
        W_on_c_am_feas(HOUR) = sum([W_on_c_am_feas1(HOUR) W_on_c_am_feas2(HOUR) W_on_c_am_feas3(HOUR) W_on_c_am_feas4(HOUR) W_on_c_am_feas5(HOUR)]);
        W_of_c_feas(HOUR) = max([W_of_c_feas1(HOUR) W_of_c_feas2(HOUR) W_of_c_feas3(HOUR) W_of_c_feas4(HOUR) W_of_c_feas5(HOUR)]);
        W_of_c_am_feas(HOUR) = sum([W_of_c_am_feas1(HOUR) W_of_c_am_feas2(HOUR) W_of_c_am_feas3(HOUR) W_of_c_am_feas4(HOUR) W_of_c_am_feas5(HOUR)]);
        S_c_feas(HOUR) = max([S_c_feas1(HOUR) S_c_feas2(HOUR) S_c_feas3(HOUR) S_c_feas4(HOUR) S_c_feas5(HOUR)]);
        S_c_am_feas(HOUR) = sum([S_c_am_feas1(HOUR) S_c_am_feas2(HOUR) S_c_am_feas3(HOUR) S_c_am_feas4(HOUR) S_c_am_feas5(HOUR)]);
        S_c_feas_s(HOUR) = max([S_c_feas1_s(HOUR) S_c_feas2_s(HOUR) S_c_feas3_s(HOUR) S_c_feas4_s(HOUR) S_c_feas5_s(HOUR)]);
        S_c_am_feas_s(HOUR) = sum([S_c_am_feas1_s(HOUR) S_c_am_feas2_s(HOUR) S_c_am_feas3_s(HOUR) S_c_am_feas4_s(HOUR) S_c_am_feas5_s(HOUR)]);
    end
    
    if SUCCESS == 0
        return
    end
else
    FEASIBLE_STATES_NUM = find((GMINlst <= DEMAND(HOUR) - RES_DN(HOUR) - RES_DN_CCS(HOUR)) & (DEMAND(HOUR) + RES_UP(HOUR) + RES_UP_CCS(HOUR) <= GMAXlst));
    if isempty(FEASIBLE_STATES_NUM)         % if there are no feasible states found
        SUCCESS = 0;                        % prepare for program termination
        STR = ['NO FEASIBLE STATES FOR HOUR ',num2str(HOUR),'! PROGRAM TERMINATES!'];
        msgbox(STR,'NO FEASIBLE STATES','warn');
        return
    else
        SUCCESS = 1;
    end
end