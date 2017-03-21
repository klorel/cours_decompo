reset;
model study.mod;
data study.dat;

problem simulation:
	simulation_cost,
	p, lol, p_new, p_interco,
	ctr_eod, ctr_pmax, ctr_pmax_new;
	
problem coordination: p_max_new;


param LB default -1e20;
param UB default +1e20;
param BEST_UB default +1e20;
param EPS := 1e-3;

param CURRENT_ITE default 0;
param STOP default 0;

var alpha >= -1e10;

set ITERATIONS default {};
param CUT_RHS{ITERATIONS}				default 0;
param CUT_COEFF{ITERATIONS, AREA, NEW_PROD}	default 0;
param CUT_POINT{ITERATIONS, AREA, NEW_PROD}	default 0;


set OTHER_ITERATIONS default {};
param OTHER_RHS{OTHER_ITERATIONS}				default 0;
param OTHER_COEFF{OTHER_ITERATIONS, AREA, NEW_PROD}	default 0;
param OTHER_POINT{OTHER_ITERATIONS, AREA, NEW_PROD}	default 0;

param OTHER_CENTER{AREA, NEW_PROD}	default 0;
param OTHER_LB{AREA, NEW_PROD} default 0;
param OTHER_UB{AREA, NEW_PROD} default 0;

subject to p_max_max{area in AREA, prod in NEW_PROD}: 0 <= p_max_new[area, prod]<= 60e3;

subject to ctr_cut{ite in ITERATIONS}:
	alpha 
	>=
	+CUT_RHS[ite]	
	+sum{area in AREA,prod in NEW_PROD}CUT_COEFF[ite, area, prod]*(p_max_new[area, prod]-CUT_POINT[ite, area, prod])
		
	+sum{area in AREA,prod in NEW_PROD}POWER_COST[prod]*p_max_new[area, prod]	
	;
	
subject to ctr_other_cut{ite in OTHER_ITERATIONS:1=0}:
	alpha 
	>=
	+OTHER_RHS[ite]	
	+sum{area in AREA,prod in NEW_PROD}OTHER_COEFF[ite, area, prod]*(p_max_new[area, prod]-OTHER_POINT[ite, area, prod])
	;
	
problem other_simulation:
	overall_cost, 
	p, lol, p_new, p_interco,
	ctr_eod, ctr_pmax, ctr_pmax_new, p_max_new;
	
subject to lower_p_max_new{area in AREA, prod in NEW_PROD}:  p_max_new[area, prod] >= max(OTHER_CENTER[area, prod]-OTHER_LB[area, prod], 0);
subject to upper_p_max_new{area in AREA, prod in NEW_PROD}:  p_max_new[area, prod] <= min(60e3, OTHER_CENTER[area, prod]+OTHER_UB[area, prod]);

problem coordination;