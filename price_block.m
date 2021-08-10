%price block script
if STACKING_OPTION==2 % short-run marginal cost (SRMC)
    MERIT_COST       = 0*GNLC + GFC.*(1./GINC) + VOM_g + ECO2_g.*(1./GINC).*CoC_g.*(1-Y_MAX_CAPT_g) + ECO2_g.*Y_MAX_CAPT_g.*(1./GINC).*(CC_VOM_g + CC_TRANS_g) - STK_g;
else %average full-load cost
    MERIT_COST       = (0*GNLC + GFC.*GMAX.*GINC/1000)./GMAX;
end
Merit_Order  = cumsum(GMAX);
%
Price     = zeros(NT,1);
Status    = ones(NT,1)*Merit_Order';
Start_up  = zeros(NT,NG);
Shut_dn   = zeros(NT,NG);
%
for i = 1:NT %which units are on or off
    for j = 1:NG
        if Status(i,j) > DEMAND(i) + RES_UP(i)/2
            Status(i,j) = 0;
        else
            Status(i,j) = 1;
        end
    end
end
%
zz=zeros(length(Status(:,1)),1);
up_down=[zz diff(Status,1,2)];
Start_up(up_down>0)=1;
Start_up(up_down<=0)=0;
Shut_dn(up_down<0)=1;
Shut_dn(up_down>=0)=0;
%
Start_up = [Start_up;zeros(1,NG)];% to make it 168 hours if run for 7 days 
Shut_dn = [Shut_dn;zeros(1,NG)];  % to make it 168 hours if run for 7 days
%
Op_hours_Num_g = sum(Status,1)';  % number of hours in operation 
Start_up_Num_g = sum(Start_up,1)';% number of start-ups 
Op_hours_Num_g(Op_hours_Num_g==0)=1;
Start_up_Num_g(Start_up_Num_g==0)=1;
%
%start-up and shut-down costs
GEN_START_SHUT_COST = (((GSC + (SU_FUEL_COLD_g.*GFC + SU_CO2_COLD_g.*CoC_g.*(1-Y_MAX_CAPT_g)))./2 + (GSDC + SD_FUEL_g.*GFC + ...
    SD_CO2_g.*CoC_g.*(1-Y_MAX_CAPT_g))).*Start_up_Num_g)./(Op_hours_Num_g.*GMAX_0);
%
Uplift_g = zeros(NG,1);
%
MC = zeros(NT,1);
SC = zeros(NT,1);
UL = zeros(NT,1);
UL_e = zeros(NT,1);
UL_max = zeros(NT,1);
CR = zeros(NT,1);
%
x_g = 1:NG;
Uplift_g(x_g) = ((Uplift_k.*Merit_Order(x_g))./(Uplift_SUM_G-(Merit_Order(x_g))));
%
%price calcualtion block
for HOUR = 1:NT
    Merit_Index = Merit_Order(DEMAND(HOUR)+RES_UP(HOUR) > Merit_Order);
    MG = min(length(Merit_Index)+1,length(GMIN));
    P_MARGINAL = GMIN(MG);                                                                % marginal power array
    MC(HOUR) = MERIT_COST(MG);                                                            % marginal cost at each hour
    SC(HOUR) = GEN_START_SHUT_COST(MG);                                                   % start-up and shut-down cost at each hour
    UL(HOUR) = ((Uplift_k*DEMAND(HOUR))/(Uplift_SUM_G - (DEMAND(HOUR))));                 % hyperbolic price mark-up uplift
    UL_e(HOUR) = Uplift_VOLL*exp(-Uplift_b*((Uplift_SUM_G - DEMAND(HOUR))/Uplift_SUM_G)); % exponential price mark-up uplift
    UL_max(HOUR) = max(UL(HOUR),UL_e(HOUR));
    %
    % price mark-down is equal to -1 ROC = - On_subsidy = - 50, when the amount of...
    %...available wind generation (W_t_Original) is equal to demand (D_t)
    CR(HOUR) = -On_subsidy*exp(Uplift_c*((W_t_Original(HOUR)-D_t(HOUR)))/D_t(HOUR));
    Price(HOUR) =  MC(HOUR) + SC(HOUR) + UL_max(HOUR) + CR(HOUR);
end
Price(Price<-Off_subsidy) = -Off_subsidy;   % electricity prices are capped from falling below 100 £/MWh
Price_Original = Price;                     % finishes simulating the price
