struct holdmaxintervalsum{
	int lowcenter;   /* first three are true centers (actual sequence location */
	int highcenter;  /* times two) */
	int testcenter;
	double matches_in_interval; /* number of matches counted from all centers in interval */
	int tupleindex;          /* g in 1..NTS */
	int intervalindex;       /* h in 1..numberofintervals[g] */
} TestMaxIntervalSum;

