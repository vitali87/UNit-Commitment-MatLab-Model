%priority_list
function [GMAXcum,GMINcum,LIST_STATES,LIST_INDEX] = priority_list(GNLC,GFC,GMAX,GMIN,GINC,NG,ECO2_g,CoC_g,Y_MAX_CAPT_g,STK_g,VOM_g,CC_VOM_g,CC_TRANS_g,STACKING_OPTION)
% --------------------------------------------------------------------------------------------
% Creates the list of states, where each state has one unit commited more than
% the previous state. Generators are ordered according to their full load 
% average cost or short-run marginal cost bases (VA),the least expensive coming first.
% The list is column-based: the 1st column contains only the cheapest gen.;
% the 2nd column contains two cheapest etc.; the last one containts all NG gen.
% OUTPUT:
% LIST_STATES [NG x NG] - matrix of states; each column represents one state
% LIST_INDEX [NG x 1]   - order of generator indices (least expensive gen. is the first)
% GMINcum [NG x 1]      - total min. generator output for each state
% GMAXcum [NG x 1]      - total max. generator output for each state
% Example: if the order of generators average costs gives G3,G1,G2, then:
% LIST_INDEX = 3,1,2
% LIST_STATES = [0 1 1       - 1st column has only          G3 commited
%                0 0 1       - 2nd column has commited      G3+G1
%                1 1 1]      - 3rd column has commited      G3+G1+G2
%--------------------------------------------------------------------------------------------
if STACKING_OPTION==2 % short-run marginal cost (SRMC)
    MERIT_COST       = 0*GNLC + GFC.*(1./GINC) + VOM_g + ECO2_g.*(1./GINC).*CoC_g.*(1-Y_MAX_CAPT_g) + ECO2_g.*Y_MAX_CAPT_g.*(1./GINC).*(CC_VOM_g + CC_TRANS_g) - STK_g;
else %average full-load cost
    MERIT_COST = (0*GNLC + GFC.*GMAX.*GINC/1000)./GMAX;
end
[M,LIST_INDEX] = sort(MERIT_COST);                    % sort them (make a priority list)
LIST_STATES = triu(ones(NG));                           % prepare matrix of the states [NGxNG]
LIST_STATES(LIST_INDEX,:) = LIST_STATES(1:NG,:);        % create the list of states
LIST_STATES = logical(LIST_STATES);                     % not neccessarily, but it is a good manner
GMAXcum = cumsum(GMAX(LIST_INDEX));                     % create total max. generator output for each state
GMINcum = cumsum(GMIN(LIST_INDEX));                     % create total min. generator output for each state
%prints_states(NG,GMINcum,GMAXcum,LIST_STATES)% can be activated if you
%want to display the results on matlab command window
end