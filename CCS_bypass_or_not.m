%CCS_bypass_or_not script
Y_CAPT_g_t          = zeros(NT,NG);
P_MIN_CAPT_g_t      = zeros(NT,NG);
P_MAX_CAPT_g_t      = zeros(NT,NG);
u_CAPT_g_t          = zeros(NT,NG);
R_SP_UP_CAPT_t      = zeros(NT,1);
R_SP_DN_CAPT_t      = zeros(NT,1);
Y_CAPT_g            = zeros(1,NG);
if FLEXIBLE_CO2_CAPTURE_OPTION == 1  % comparison bypass or not one by one for ccs plants
    for HOUR = 1:NT
        Y_MIN_CAPT_g = Y_MAX_CAPT_g.*K_CCS_LB_CAPTURE_RATE;% minimum capture rate is a proportion of the max capture rate 
        Profit_Y_MAX_CAPT_g = (GMAX - (CC_FIXED_g + CC_OP_g.*(((GMAX.*m_C_g)+c_C_g).*ECO2_g).*Y_MAX_CAPT_g)).*Price(HOUR) - ((GMAX.*m_C_g)+...
            c_C_g).*GFC -((GMAX.*m_C_g)+c_C_g).*ECO2_g.*CoC_g.*(1-Y_MAX_CAPT_g) - GMAX.*VOM_g - ...
            ((GMAX.*m_C_g)+c_C_g).*ECO2_g.*Y_MAX_CAPT_g.*CC_VOM_g -((GMAX.*m_C_g)+c_C_g).*ECO2_g.*Y_MAX_CAPT_g.*CC_TRANS_g - ...
            ((GMAX.*m_C_g)+c_C_g).*ECO2_g.*Y_MAX_CAPT_g.*CC_CSOLV_g.*(((CC_SOLVD_g - CC_SOLVD_TH_g).*Y_MAX_CAPT_g) + CC_SOLVD_TH_g);
        Profit_Y_MIN_CAPT_g = (GMAX - (CC_FIXED_g + CC_OP_g.*(((GMAX.*m_C_g)+c_C_g).*ECO2_g).*Y_MIN_CAPT_g)).*Price(HOUR) - ((GMAX.*m_C_g)+c_C_g).*GFC ...
            - ((GMAX.*m_C_g)+c_C_g).*ECO2_g.*CoC_g.*(1-Y_MIN_CAPT_g) - GMAX.*VOM_g - ((GMAX.*m_C_g)+c_C_g).*ECO2_g.*Y_MIN_CAPT_g.*CC_VOM_g - ...
            ((GMAX.*m_C_g)+c_C_g).*ECO2_g.*Y_MIN_CAPT_g.*CC_TRANS_g - ((GMAX.*m_C_g)+c_C_g).*ECO2_g.*Y_MAX_CAPT_g.*CC_CSOLV_g.*(((CC_SOLVD_g - ...
            CC_SOLVD_TH_g).*Y_MIN_CAPT_g) + CC_SOLVD_TH_g);
        % Capture rate is either 90% or 0%
        for i = 1:NG
            if Profit_Y_MAX_CAPT_g(i) > Profit_Y_MIN_CAPT_g(i)
                Y_CAPT_g(i)     = Y_MAX_CAPT_g(i);
                P_MAX_CAPT_g(i) = round(CC_FIXED_g(i) + (CC_OP_g(i).*(((GMAX(i).*m_C_g(i))+c_C_g(i)).*ECO2_g(i)).*Y_MAX_CAPT_g(i)));
                P_MIN_CAPT_g(i) = round(CC_FIXED_g(i) + (CC_OP_g(i).*(((GMIN(i).*m_A_g(i))+c_A_g(i)).*ECO2_g(i)).*Y_MAX_CAPT_g(i)));
            else
                Y_CAPT_g(i)     = 0;
                P_MAX_CAPT_g(i) = CC_FIXED_g(i);
                P_MIN_CAPT_g(i) = round(CC_FIXED_g(i) + (CC_OP_g(i).*(((GMIN(i).*m_A_g(i))+c_A_g(i)).*ECO2_g(i)).*Y_MAX_CAPT_g(i)));
            end
        end
        Y_CAPT_g_t(HOUR,:) = Y_CAPT_g';
        P_MAX_CAPT_g_t(HOUR,:) = P_MAX_CAPT_g';
        P_MIN_CAPT_g_t(HOUR,:) = P_MIN_CAPT_g';
        %Assumes base power plant is at full output
        R_SP_UP_CAPT_t(HOUR) = sum(CC_OP_g.*(((GMAX.*m_C_g)+c_C_g).*ECO2_g.*Y_CAPT_g'));
        R_SP_DN_CAPT_t(HOUR) = sum(CC_OP_g.*(((GMAX.*m_C_g)+c_C_g).*ECO2_g.*Y_MAX_CAPT_g)) - sum(CC_OP_g.*(((GMAX.*m_C_g)+c_C_g).*ECO2_g.*Y_CAPT_g'));
        RES_UP_CCS = zeros(NT,1); % Adjustment/contribution of CCS units to upwards spinning reserve in feasible states
        RES_DN_CCS = zeros(NT,1); %Adjustment/contribution of CCS units to downwards spinning reserve in feasible states
    end
else
    Y_CAPT_g_t = ones(NT,1)*Y_MAX_CAPT_g';
    P_MAX_CAPT_g_t = ones(NT,1)*round(CC_FIXED_g + CC_OP_g.*(((GMAX.*m_C_g)+c_C_g).*ECO2_g.*Y_MAX_CAPT_g))';
    P_MIN_CAPT_g_t = ones(NT,1)*round(CC_FIXED_g + CC_OP_g.*(((GMIN.*m_A_g)+c_A_g).*ECO2_g.*Y_MAX_CAPT_g))';
    %Adjustment/contribution of CCS units to upwards spinning reserve in feasible states
    RES_UP_CCS = ones(NT,1)*sum(CC_OP_g.*(((GMAX.*m_C_g)+c_C_g).*ECO2_g.*Y_MAX_CAPT_g));
    %Adjustment/contribution of CCS units to downwards spinning reserve in feasible states
    RES_DN_CCS = zeros(NT,1);
end
u_CAPT_g_t = round(Y_CAPT_g_t);
%
if MIN_UP_DOWN_TIME_FLAG == 0
    GMINUP(:) = 0;
    GMINDOWN(:) = 0;
end
%
if RAMP_UP_DOWN_FLAG == 0
    GRAMPUP(:) = Inf;
    GRAMPDOWN(:) = Inf;
end
%
GMINUP_Original  = GMINUP;
GMINDOWN_Original  = GMINDOWN;
%
if COMPLETE_ENUMERATION_FLAG == 0 
    [GMAXlst,GMINlst,LIST_STATES,GEN_ORDER] = priority_list(GNLC,GFC,GMAX,GMIN,GINC,NG,ECO2_g,CoC_g,Y_MAX_CAPT_g,STK_g,VOM_g,CC_VOM_g,CC_TRANS_g,STACKING_OPTION);
elseif COMPLETE_ENUMERATION_FLAG == 1
    [GMAXlst,GMINlst,LIST_STATES,GEN_ORDER] = complete_enumeration(GNLC,GFC,GMAX,GMIN,GINC,NG,ECO2_g,CoC_g,Y_MAX_CAPT_g,STK_g,VOM_g,CC_VOM_g,...
        CC_TRANS_g,STACKING_OPTION);
else 
    [GMAXlst,GMINlst,LIST_STATES,GEN_ORDER] = dynamic_flex(GNLC,GFC,GMAX,GMIN,GINC,NG,ECO2_g,CoC_g,Y_MAX_CAPT_g,STK_g,VOM_g,CC_VOM_g,...
        CC_TRANS_g,FLEX_g,STACKING_OPTION); 
end
%
R_SP_UP_D_t = zeros(NT,1);
R_SP_DN_D_t = zeros(NT,1);
R_SP_UP_W_t = zeros(NT,1);
R_SP_DN_W_t = zeros(NT,1);
R_SP_UP_S_t = zeros(NT,1);
R_SP_DN_S_t = zeros(NT,1);
%
x = min_period:max_period; 
R_SP_UP_D_t(x) = D_u*(D_t(x)); % Demand forecast uncertainty (std) is D_u*100% of current demand
R_SP_DN_D_t(x) = D_u*(D_t(x)); % Demand forecast uncertainty (std) is D_u*100% of current demand
R_SP_UP_W_t(x) = min(lambda_up*W_u*(W_t(x)),W_t(x) - W_c_t(x));
R_SP_DN_W_t(x) = W_u*(W_t(x));  % Wind uncertainty is W_u*100% of forecast wind
R_SP_DN_S_t(x) = S_u*(S_t(x));  % Solar forecast uncertainty is S_u*100% of forecast solar
% R_SP_UP_S_t(x) = min(S_u*S_t(x),max(0,S_t(x) - S_c_t(x))); 
R_SP_DN_S_t(x) = S_u*S_t(x); % Solar forecast uncertainty is S_u*100% of forecast solar
%
R_SP_UP_W_t(R_SP_UP_W_t<0) = 0;
R_SP_UP_S_t(R_SP_UP_S_t<0) = 0; % In case solar reserves are negative
%
R_SP_UP_D_W_S_t = ((R_SP_UP_D_t).^2 + (R_SP_UP_W_t).^2 + (R_SP_UP_S_t).^2).^0.5;
R_SP_DN_D_W_S_t = ((R_SP_DN_D_t).^2 + (R_SP_DN_W_t).^2 + (R_SP_DN_S_t).^2).^0.5;
%
RES_UP = R_SP_UP_MAX + lambda_up*R_SP_UP_D_W_S_t;
RES_UP(RES_UP<0) = 0;
RES_DN = R_SP_DN_MAX + lambda_dn*R_SP_DN_D_W_S_t;
RES_DN(RES_DN<0) = 0;