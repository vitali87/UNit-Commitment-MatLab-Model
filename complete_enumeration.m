function [GMAXlst,GMINlst,LIST_STATES,GEN_ORDER] = complete_enumeration(GNLC,GFC,GMAX,GMIN,GINC,NG,ECO2_g,CoC_g,Y_MAX_CAPT_g,STK_g,...
    VOM_g,CC_VOM_g,CC_TRANS_g,STACKING_OPTION)
% --------------------------------------------------------------------------------------------
% Creates the complete list of states (for NG generators, there are totally 2^NG states)
% The list is column-based (each column represents a state).
% OUTPUT:
% LIST_STATES [NG x 2^NG] - matrix of states; each column represents one state
% GEN_ORDER [NG x 1]      - order of generators to be commited (least expensive gen. the first)
% GMINlst [2^NG x 1]      - total min. generator output for each state
% GMAXlst [2^NG x 1]      - total max. generator output for each state
% Example: if there are 3 generators, then there are totally 8 states:
% LIST_STATES = [0 0 0 0 1 1 1 1
%                0 0 1 1 0 0 1 1
%                0 1 0 1 0 1 0 1]
%--------------------------------------------------------------------------------------------
if STACKING_OPTION==1
    MERIT_COST = (0*GNLC + GFC.*GMAX.*GINC/1000)./GMAX; % Calculate full load average cost for each unit
else
    % Calculate short-run marginal cost for each unit
    MERIT_COST=0*GNLC + GFC.*(1./GINC) + VOM_g + ECO2_g.*(1./GINC).*CoC_g.*(1-Y_MAX_CAPT_g) + ECO2_g.*Y_MAX_CAPT_g.*(1./GINC).*(CC_VOM_g + CC_TRANS_g) - STK_g;
end
[M,GEN_ORDER] = sort(MERIT_COST);                   % sort them (make a priority list of gen. commitment)
LIST_STATES=dec2bin(0:2^NG-1)';                       % all possible combinations of NG generators
LIST_STATES = logical(sscanf(LIST_STATES,'%1d',size(LIST_STATES)));

GMINlst = LIST_STATES.' * GMIN; % for each state (combination of generators) in the list,
GMAXlst = LIST_STATES.' * GMAX; % find the min. and max. power of the combination
% next 3 lines are not neccessary, but it is nice to see max. possible output of generators in an increasing order
[GMAXlst,INDEX]=sort(GMAXlst);      % sort the states according to increasing total max. power
GMINlst = GMINlst(INDEX);           % re-order min. power accordingly
LIST_STATES = LIST_STATES(:,INDEX); % and re-order the list of states as well
%
%prints_states(NG,GMINlst,GMAXlst,LIST_STATES) % optional
end