%% Input Data for Dynamic Programming based Unit Commitment (v_9)
% Author: Vitali Avagyan, June 2017
%% Initialisation and flagging
Day_Start                   = 0;            % which day to start
Time_Start                  = Day_Start*24; % t_0
Days                        = 365;          % days to schedule
Time_Length                 = Days*24;
%
MIN_LOAD                    = 8500;        % min load requirement for the system inertia (MW)
K_WIND_INERTIA              = -0.1;            % Proportion of wind contribution towards synthetic inertia. This gets multiplied by the...
                                            % ...installed wind capacity and added to the 
%                                           % minimum load, so it should bear a minus sign.
P_MAX                       = 1800;         % largest online generator (credible loss) (MW)
%
S_u                         = 0.08;         % Solar forecast uncertainty coefficient(IEA, 2011)
W_u                         = 0.1;          % Wind forecast uncertainty coefficient
D_u                         = 0.01;         % Demand forecast uncertainty coefficient
%
S_large                     = 0.48;         % Proportion of large scale solar PV (Solar PV deployment from BEIS, 2017)
%
A_w_on                      = 0.98;         % Availability factor of onshore wind
A_w_of                      = 0.90;         % Availability factor of offshore wind
%
D_index                     = 1;            % Demand database index: 1 - INDO, 2 - IO14_DEM
S_index                     = 20;           % Solar database index: 19 - MERRA1, 20 - MERRA2, 21 - Sarah, 22 - Mean_MERRA2
year_index                  = 10;           % year data: 2 - 2002, 3 - 2003, 4 - 2004, 5 - 2005, 6 - 2006, 7 - 2007, 8 - 2008, 9 - 2009, 10 - 2010
%
K_CCS_LB_CAPTURE_RATE       = 0;            % min capture rate of CCS as a proportion of its max capture rate
%
ins_solar_MW                = 30000;        % installed capacity of solar (MW)
ins_wind_on_of_GW           = [11 12];      % [3 4] - 10GW, [5 6] - 15GW, [7 8] - 20GW, [9 10] - 25GW, [11 12] - 30GW, [13 14] - 35GW, [15 16] - 40GW, [17 18] - 45GW % installed capacity of wind
%
Uplift_k                    = 5;            % Scale factor for hyperbolic price mark-up for simulating electricity price
Uplift_VOLL                 = 4000;         % National Grid [2013] assumes VOLL is equal to 4000 ï¿½/MWh which is approximately the GDP per unit of electricity consumed
Uplift_b                    = 20;           % value of 20 gives results that agree well with Eager (2011) and Grunewald (2012)
Uplift_c                    = 10;           % value of 10 adjusts the shape of the price mark-down function so that negative prices only occur as frequently as reported in Baringa (2015)
%
Off_subsidy                 = 100;          % subsidy for electricity generation from offshore wind which is equal to 2*ROC (ï¿½/MWh)
On_subsidy                  = 50;           % subsidy for electricity generation from onshore which is equal to 1*ROC (ï¿½/MWh)
% flags
MIN_UP_DOWN_TIME_FLAG       = 1;            % take minimum up and down times into account (1) or not (0)
RELAXATION                  = 1;            % 0 - terminate if infeasible, 1 - relax up/down time if infeasible
RAMP_UP_DOWN_FLAG           = 1;            % take ramp up and down rates into account (1) or not (0)
RAMP_COSTS_INCLUDED         = 1;            % 1 - ramp costs included in the system costs, 0 - ramp costs excluded from the system costs
N_PRED                      = 1;            % number of predecesors to be searched (N_PRED >= 1)
COMPLETE_ENUMERATION_FLAG   = 0;            % 2 - flexible enumeration, 1 - complete enumeration, 0 - priority list
STACKING_OPTION             = 2;            % 1 - average full-load cost, 2 - short-run marginal cost (SRMC)
DETAIL_PRINT_FLAG           = 0;            % detailed results printing: 0 - no, 1 - yes 
DISPATCH_METHOD             = 0;            % 0 - piece-wise quadratic approximation, 1 - quadprog, 2 - linprog/cplex, 3 - quick dispatch
LINPROG_OR_CPLEX            = 0;            % 1 - linprog,0 - cplex
RESERVE_FLAG                = 1;            % take spinning reserve in calculation (1) or not (0)
START_UP_COST_METHOD        = 3;            % 1 - cold start-up (const), 2 - cold/hot start-up, 3 - exponential start-up
EXP_COST_OPTION             = 1;            % 0 - all costs inside the brackets, 1 - start-up fixed cost outside the brackets, 2 - "cold" alpha-beta approach, 3 - "hot" alpha-beta approach, 
FLEXIBLE_CO2_CAPTURE_OPTION = 0;            % 0 - constant capture, 1 - bypass allowed
SEQUENTIAL_RES_CURTAILMENT  = 1;            % 0 - no RES curtailment if infeasible, 1 - RES curtailment if infeasible
%% Data on demand, wind, solar and generators
y_start          = 2 + 8760*(year_index-2) + floor((year_index-1)/4)*24;
y_end            = 1 + 8760*(year_index-1) + floor(year_index/4)*24;
D_2010           = csvread('DP_DEMAND_WIND_2002_2010_Solar.csv',y_start,D_index,[y_start D_index y_end D_index]);
W_e_on_2010      = csvread('DP_DEMAND_WIND_2002_2010_Solar.csv',y_start,ins_wind_on_of_GW(1),[y_start ins_wind_on_of_GW(1) y_end ins_wind_on_of_GW(1)])*A_w_on; 
W_e_of_2010      = csvread('DP_DEMAND_WIND_2002_2010_Solar.csv',y_start,ins_wind_on_of_GW(2),[y_start ins_wind_on_of_GW(2) y_end ins_wind_on_of_GW(2)])*A_w_of; 
S_2010           = csvread('DP_DEMAND_WIND_2002_2010_Solar.csv',y_start,S_index,[y_start S_index y_end S_index])*ins_solar_MW;
gen_data         = csvread('UCED_Parameters_v9.csv',3,0,'A4..BV71');
%% Generator column data
GMIN            = gen_data(:,2);            % generator min power                                (MW)
GMAX            = gen_data(:,3);            % generator max. power                               (MW)
P_MIN_CAPT_g    = gen_data(:,4);            % generator minimum power output without capture     (MW_e)
P_MAX_CAPT_g    = gen_data(:,5);            % generator max imum power output without capture    (MW_e)
GSTATINI        = gen_data(:,6);            % generator initial status (time on/off)             (h)
GMINUP          = gen_data(:,7);            % generator min. up time                             (h)
GMINDOWN        = gen_data(:,8);            % generator min. down time                           (h)
GMINOFF         = gen_data(:,52);           % generator min. offline time                        (h)
GRAMPUP         = gen_data(:,9);            % generator ramp up rate                             (MW/h)
GRAMPDOWN       = gen_data(:,10);           % generator ramp down rate                           (MW/h)
GNLC            = gen_data(:,11);           % generator no load cost                             (ï¿½/h)
GFC             = gen_data(:,12);           % generator fuel cost                                (ï¿½/MBTU)
GFC_TH          = gen_data(:,13);           % generator fuel cost                                (ï¿½/MWhth)
VOM_g           = gen_data(:,14);           % generator operating and maintenance cost           (ï¿½/MWh_e)
GINC            = gen_data(:,15);           % incremental heat rate                              (BTU/kWh)
GINC_E_TO_TH    = gen_data(:,16);           % incremental heat rate                              (MW_e/MW_th)
GCSTIME         = gen_data(:,17);           % generator cold start time                          (h)
GSDTIME         = gen_data(:,22);           % generator shut down time                           (h)
GSC             = gen_data(:,18);           % generator start up cost (cold), also BETA          (ï¿½)
SU_FUEL_COLD_g  = gen_data(:,19);           % generator fuel consumption during start-up (cold)  (MWh_th) 
SU_CO2_COLD_g   = gen_data(:,21);           % generator CO2 emissions during start-up cold       (tCO2)
GSDC            = gen_data(:,23);           % generator shut down cost                           (ï¿½)
SD_FUEL_g       = gen_data(:,24);           % generator shut down fuel consumption               (MWh_th)
GSH             = gen_data(:,25);           % generator start up cost (hot), also ALPHA          (ï¿½)
SD_CO2_g        = gen_data(:,26);           % generator CO2 emissions during shut-down           (tCO2)
TAU             = gen_data(:,27);           % generator start up cost exp. coef.                 (ï¿½) modified by (VA)
ECO2_g          = gen_data(:,28);           % Specific ECO2_g emissions                          (tCO2/MWh_th)
CoC_g           = gen_data(:,29);           % Cost of Carbon                                     (ï¿½/tCO2)
%
Y_MAX_CAPT_g    = gen_data(:,30);           % Maximum capture rate                               (%)
CC_FIXED_g      = gen_data(:,31);           % Capture penalty fixed (Fixed CO2 Capture Costs)    (MW_e)
CC_OP_g         = gen_data(:,32);           % Capture penalty operating                          (MW_e/tCO2)
CC_VOM_g        = gen_data(:,33);           % Capture variable operating and maintenance cost    (ï¿½/tCO2)
CC_CSOLV_g      = gen_data(:,34);           % Capture solvent cost                               (ï¿½/kg)
CC_SOLVD_g      = gen_data(:,35);           % Capture solvent degradation rate                   (kg/tCO2)
CC_SOLVD_TH_g   = gen_data(:,36);           % Capture solvent thermal degradation rate           (kg/tCO2)
CC_TRANS_g      = gen_data(:,37);           % Capture transport and storage cost                 (ï¿½/tCO2)
STK_g           = gen_data(:,38);           % generator strike price                             (ï¿½/MWh_e)
%
RAMP_UP_g       = gen_data(:,39);           % generator ramping up cost                          (ï¿½/MW_e)
RAMP_DN_g       = gen_data(:,40);           % generator ramping dn cost                          (ï¿½/MW_e)
m_A_g           = gen_data(:,41);           % Fuel gradient coefficient A                        (MWh_th/MWh_e)
m_B_g           = gen_data(:,42);           % Fuel gradient coefficient B                        (MWh_th/MWh_e)
m_C_g           = gen_data(:,43);           % Fuel gradient coefficient C                        (MWh_th/MWh_e)
c_A_g           = gen_data(:,44);           % Fuel intercept coefficient A                       (MWh_th)
c_B_g           = gen_data(:,45);           % Fuel intercept coefficient B                       (MWh_th)
c_C_g           = gen_data(:,46);           % Fuel intercept coefficient C                       (MWh_th)
Inertia_g       = gen_data(:,47);           % Inertia constant                                   (MJ/MW or s)
%
COEF_A          = gen_data(:,48);           % free term in quadratic-cost function               (ï¿½)
COEF_B          = gen_data(:,49);           % linear term in quadratic-cost function             (ï¿½/MWh]
COEF_C          = gen_data(:,50);           % 2nd order term in quadratic-cost function          (ï¿½/(MW^2)h]
%
FLEX_g          = gen_data(:,51);           % Flexible plant binary indicator
%
Uplift_SUM_G    = sum(GMAX);                % Installed capacity                                 (MW)
%% Unifying block 
if STACKING_OPTION==2
    GFC=GFC_TH;
    GINC=GINC_E_TO_TH;
end
for x=1:length(GCSTIME)
    if GCSTIME(x)~=0
        GMINDOWN(x)=GSDTIME(x)+GMINOFF(x); % this is to facilitate the modelling process as some of the notations are based on GMINDOWN
    end
end