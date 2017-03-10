set AREA_TS  dimen 2;

set TS   := union{(area, ts) in AREA_TS}{ts};
set AREA := union{(area, ts) in AREA_TS}{area};

param MIN_TS := min{ts in TS}ts;
param MAX_TS := max{ts in TS}ts;


set WEEK := {week in TS: (week+1)*168<=MAX_TS and week<2};

param N_WEEK := card(WEEK);
param N_PDT := 1;
set WEEK_START_END := setof{week in WEEK}(week, week*168, week*168+N_PDT);

set USED_TS := union{(week, start, end) in WEEK_START_END}{start..end-1};

param CONSO {AREA_TS};
param EOLIEN{AREA_TS};


set PROD;

param OLD_PMAX{AREA, {"THERMIQUE", "NUCLEAIRE"} };

param ENERGY_COST{PROD};
param POWER_COST {PROD};

param VALUE_OF_LOSS_LOAD := 20e3;
 
set NEW_PROD within PROD;

param PMAX_AREA{area in AREA, prod in PROD} := 
	if prod == "NUCLEAIRE" then 
		OLD_PMAX[area, "NUCLEAIRE"]
	else if prod == "CHARBON" then
		0.3*OLD_PMAX[area, "THERMIQUE"]
	else if prod == "CCGT" then
		0.5*OLD_PMAX[area, "THERMIQUE"]
	else if prod == "TAC" then
		0.2*OLD_PMAX[area, "THERMIQUE"]
	else 1/0;

set INTERCO dimen 2;

set NEIGHBOR{area in AREA} := 
	union{(area, neighbor) in INTERCO}{neighbor}
	union
	union{(neighbor, area) in INTERCO}{neighbor}
	;

param NB_NEIGHBOR{area in AREA} :=card(NEIGHBOR[area]);

###
# full model
###
var p{AREA, PROD, USED_TS} >= 0;

var lol{AREA, USED_TS} >= 0;

var p_new{AREA, NEW_PROD, USED_TS} >= 0;

var p_max_new{AREA, NEW_PROD};

var p_interco{INTERCO, USED_TS} >= 0, <=1e3;

minimize overall_cost:
	(
		0
		+sum{area in AREA, prod in PROD, (week, start, end) in WEEK_START_END, ts in start .. end-1}p[area, prod, ts]*ENERGY_COST[prod]
		+sum{area in AREA, prod in NEW_PROD, (week, start, end) in WEEK_START_END, ts in start .. end-1}p_new[area, prod, ts]*ENERGY_COST[prod]
		+sum{area in AREA, (week, start, end) in WEEK_START_END, ts in start .. end-1}lol[area, ts]*VALUE_OF_LOSS_LOAD
		+sum{(neighbor1, neighbor2) in INTERCO, (week, start, end) in WEEK_START_END, ts in start .. end-1}1e-4 * p_interco[neighbor1, neighbor2, ts]
	) / N_WEEK
	+
	sum{area in AREA, prod in NEW_PROD}POWER_COST[prod]*p_max_new[area, prod]
	;
	
minimize simulation_cost:
	(
		0
		+sum{area in AREA, prod in PROD, (week, start, end) in WEEK_START_END, ts in start .. end-1}p[area, prod, ts]*ENERGY_COST[prod]
		+sum{area in AREA, prod in NEW_PROD, (week, start, end) in WEEK_START_END, ts in start .. end-1}p_new[area, prod, ts]*ENERGY_COST[prod]
		+sum{area in AREA, (week, start, end) in WEEK_START_END, ts in start .. end-1}lol[area, ts]*VALUE_OF_LOSS_LOAD
		+sum{(neighbor1, neighbor2) in INTERCO, (week, start, end) in WEEK_START_END, ts in start .. end-1}1e-4 * p_interco[neighbor1, neighbor2, ts]
	) / N_WEEK
	;
	
subject to ctr_eod{area in AREA, (week, start, end) in WEEK_START_END, ts in start .. end-1}:
	+lol[area, ts]
	+sum{prod in PROD}p[area, prod, ts]
	+sum{prod in NEW_PROD}p_new[area, prod, ts] 
	-sum{(area, neighbor) in INTERCO}p_interco[area, neighbor, ts]
	+sum{(neighbor, area) in INTERCO}p_interco[neighbor, area, ts]
	= 
	CONSO[area, ts]-EOLIEN[area, ts];
	

subject to ctr_pmax{area in AREA, prod in PROD, ts in USED_TS}:p[area, prod, ts] <= PMAX_AREA[area, prod];

subject to ctr_pmax_new{area in AREA, prod in NEW_PROD, ts in USED_TS}:p_new[area, prod, ts] <= p_max_new[area, prod];