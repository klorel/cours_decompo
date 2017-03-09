reset;
option solver amplxpress;
option presolve 0;

option solver_msg 0;

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

subject to p_max_max{area in AREA, prod in NEW_PROD}: 0 <= p_max_new[area, prod]<= 60e3;

subject to ctr_cut{ite in ITERATIONS}:
	alpha 
	>=
	+CUT_RHS[ite]	
	+sum{area in AREA,prod in NEW_PROD}CUT_COEFF[ite, area, prod]*(p_max_new[area, prod]-CUT_POINT[ite, area, prod])
		
	+sum{area in AREA,prod in NEW_PROD}POWER_COST[prod]*p_max_new[area, prod]	
	;
	