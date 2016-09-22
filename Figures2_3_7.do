
  use "~/Dropbox/Attribution Analysis Feb 2012/attribution.dta", clear
  
  
/* allocating DM was the largest in session 1 and randomly chosen in sessions 2-4 */

gen proposer=0
replace proposer = 1 if vw==53&treat2==1&session==1
replace proposer = 1 if vw==48&treat2==2&session==1
replace proposer = 1 if vw==38&treat2==3&session==1
replace proposer = 1 if vw==23&treat2==4&session==1
replace proposer = adm if session==2
replace proposer = adm if session==3
replace proposer = adm if session==4

gen maj=0
replace maj=1 if vw==53

/* This is Figure 7 in the Online Appendix */
  
  /* first identify those who gave none */
  
	egen overallsubid=group(session period sid)
	sort  session period sid
	egen DMid=seq(), from(1) to(5)
  
	egen totalpunishment=total(ndp), by(overallsubid)
	egen marksubject=tag(overallsubid)
	gen insample=0
	replace insample=1  if  treat1==1& treat3==1&session~=1
	keep if insample==1

	histogram totalpunishment if marksubject==1, freq title(Total Punishment Points Allocated by Subjects) xtitle("Total Punishement Points Allocated")

	tab totalpunishment if marksubject==1
	
	gen didnothing=0
	replace didnothing=1 if totalpunishment==0

	
	
	/* some basic information about patterns of subject allocations. */


keep if treat1==1&(session==2|session==3|session==4)&treat3==1 /* only costless with proposers randomized (sessions=2-4) and subjects
																  know proposer (treat3=1) and see allocation (treat1=1, but really only relbant for a 
																  session not included here (i.e., session 1).*/

egen dmid=group(treat2 vw)  /* create an id, 1-20 for the dms (5 dms in four distributions). 
								This assigns a unique id for each distribution as a whole (treat2)
								crossed with each seat weight. Since there are four distributions
								with five dm's, we have twenty unique seat weight distribution
								combinations, so 20 dm id's.  Notice that this does not respect 
								session, so if you were a dm with weight 29 in distribution = (53,29,10,6,2)
								in session 2 you get the same id as the the same dm in session 3. 

Here are the codes for each DM:

2,6,10,29,53            1          2          3          4          5 
8,11,14,19,48  		    6          7          8          9         10 
11,13,17,21,38 		   11         12         13         14         15
17,19,20,21,23  	   16         17         18         19         20

*/


/* now create an id that is unique to a subject playing a single period (so a facing a single distribution)
This can be used to calculate variables that need to be summed over just the five dm's that the subject
faces in one period. */

	egen superid=group(session period sid)
	egen tagonesub=tag(superid)



	/* identify subject-periods in which the subject allocated 0 to all DMs */


	gen zeromarker=ndp
	recode zeromarker 0=1 .=. *=0

	egen totalzeros=total(zeromarker), by(superid)

	gen allzero=0
	replace allzero=1 if totalzeros==5
	label variable allzero "1=subject allocated all 0's"


/* identify subject-periods in which the subject allocated 6 points (equallity) to all DMs */

	gen equalitymarker=ndp
	recode equalitymarker 6=1 .=. *=0

	egen totalsixes=total(equality), by(superid)

	gen allsixes=0
	replace allsixes=1 if totalsixes==5
	label variable allsixes "1=subject allocated all 6's (equality)"



/* identify subject-periods in which the subject allocated 30 to proposer */

	gen alltoprop=0
	replace alltoprop=1 if ndp==30&proposer==1
	egen totalallprop=total(alltoprop), by(superid)
	replace alltoprop=1 if totalallprop==1
	label variable alltoprop "1=subject allocated all 30 to proposer"
	
	
/*  Check to see how many cases fit these categoreis:

egen tagonesub=tag(superid)
sum  allzero allsixes alltoprop if tagonesub==1

. sum  allzero allsixes alltoprop if tagonesub==1

    Variable |       Obs        Mean    Std. Dev.       Min        Max
-------------+--------------------------------------------------------
     allzero |       672    .3556548    .4790678          0          1
    allsixes |       672    .1264881    .3326462          0          1
   alltoprop |       672     .141369    .3486614          0          1

so, 61% fit into one of these three
*/


/* identify subject-periods in which the subject allocated >0 to proposer and 0 to all other DMs */

	gen onlytopropbutnotall=0
	replace onlytopropbutnotall=1 if ndp>0&ndp~=.&proposer==1&alltoprop~=1&totalzeros==4
	egen totalonlyprop=total(onlytopropbutnotall), by(superid)
	replace onlytopropbutnotall=1 if totalonlyprop==1
	label variable onlytopropbutnotall "1=subject allocated only to proposer, but not all 30"
	
/* explore what remains */

*list superid vw proposer ndp if  allzero==0&allsixes==0&alltoprop==0&onlytopropbutnotall==0

/* this exploration revealed that some people split between proposer and biggest or proposer and majority 
    so create two variabels for these cases */

	gen propmajsplit=0
	replace propmajsplit=1 if ndp>0&ndp~=.&(proposer==1|maj==1)&alltoprop~=1&totalzeros==3
	replace propmajsplit=0 if proposer==1&maj==1
	egen totalpropmajsplit=total(propmajsplit), by(superid)
	replace propmajsplit=0
	replace propmajsplit=1 if totalpropmajsplit==2
	label variable propmajsplit "1=subject allocated to proposer and maj party only"
	
	gen biggest=0
	replace biggest=1 if vw==53|vw==48|vw==38|vw==23

	gen propbigsplit=0
	replace propbigsplit=1 if ndp>0&ndp~=.&(proposer==1|biggest==1)&alltoprop~=1&totalzeros==3
	replace propbigsplit=0 if proposer==1&biggest==1
	egen totalpropbigsplit=total(propbigsplit), by(superid)
	replace propbigsplit=0
	replace propbigsplit=1 if totalpropbigsplit==2
	label variable propbigsplit "1=subject allocated to proposer and biggest party only"
	
	egen TotalPointsUsed=total(ndp), by(superid)


/*** NOTE THAT THESE TWO VARIABLES WILL  BOTH EQUAL 1 FOR CASES OF MAJORITY and that propmajsplit==1 is a
     propoer subset of propbigsplit==1  (so only use propbigsplit in looking at overall number of cases that
	 fit into one of these categories (as below) **/
	 
	 sum  allzero allsixes alltoprop  onlytopropbutnotall propbigsplit
	 sum  allzero allsixes alltoprop  onlytopropbutnotall propbigsplit if tagonesub==1



/* explore what remains */

*list superid vw proposer ndp if  allzero==0&allsixes==0&alltoprop==0&onlytopropbutnotall==0&propbigsplit==0, sepby(superid)

/* identify subject-periods in which the subject allocated equality but less than 30 total  */


	gen PercenttoEach=ndp/TotalPointsUsed  

	gen testmarker=0
	replace testmarker=PercenttoEach if tagonesub==1
	egen marker=total(testmarker), by(superid)
	
	gen samesame=0
	replace samesame=1 if marker==PercenttoEach
	
	egen totalsame=total(samesame), by(superid)
	gen equalitylessthan30=0
	replace equalitylessthan30=1 if totalsame==5&TotalPointsUsed<30
	
	replace equalitylessthan30=0 if TotalPointsUsed==0|TotalPointsUsed==.
	
	label variable equalitylessthan30 "1=subject allocated all equal but less than 30 total"




/* identify equal but most to the proposer */

	gen testmarker2=0
	egen tagonesub2=tag(superid) if proposer==0

	replace testmarker2=PercenttoEach if tagonesub2==1

	egen marker2=total(testmarker2), by(superid)
	
	gen samesame2=0
	replace samesame2=1 if marker2==PercenttoEach&proposer==0
	replace samesame2=1 if marker2<PercenttoEach&proposer==1
	
	egen totalsame2=total(samesame2), by(superid)
	gen equalityexceptprop=0
	replace equalityexceptprop=1 if totalsame2==5&allsixes==0
	
	replace equalityexceptprop=0 if alltoprop==1
		
	replace equalityexceptprop=0 if onlytopropbutnotall==1
	
	label variable equalityexceptprop "1=subject allocated all equal but more to proposer"



/*for cases that meet none of the other criteria, do a regression of vw on ndp with a control for proposer. If the coeffcient on vw is
  positive then mark it */
  
  
  gen suballtoprop=alltoprop+onlytopropbutnotall  /* make a category that combines when you only gave to propsoer */
	gen equalany=allsixes+equalitylessthan30 /* equal to everyone even if that was less that 30 (>0) though */

	gen bbstore=.
	gen sestore=.  
	gen markthecase=1 
	replace markthecase=0 if allzero==1
	replace markthecase=0 if equalany==1
	replace markthecase=0 if suballtoprop==1 
	replace markthecase=0 if propbigsplit==1 
	replace markthecase=0 if equalityexceptprop==1


	levelsof superid if markthecase==1, local(supidloc)
    foreach hh of local supidloc {
		list vw ndp proposer if superid==`hh'
		regress ndp vw proposer if superid==`hh'
		replace bbstore=_b[vw] if superid==`hh'
		replace sestore=_se[vw] if superid==`hh'
	}


gen postivevw=0
replace postivevw=1 if bbstore>0&bbstore~=.

gen tratio=bbstore/sestore
gen positivesigvw=0
replace positivesigvw=1 if bbstore>0&tratio>=2&tratio~=.

/* marks if your in one of my categoreis */


egen oneofthem=rowtotal(allzero equalany suballtoprop propbigsplit equalityexceptpro postivevw)
gen random=oneofthem
recode random 1=0 0=1

egen oneofthem2=rowtotal(allzero equalany suballtoprop propbigsplit equalityexceptpro positivesigvw)
gen random2=oneofthem
recode random2 1=0 0=1


/***** This is FIGURE 2a from Main Text *****/

gen smallalldm=allocdm
recode smallalldm 1/9=1 10/14=2 15/25=3

# delimit ;
graph bar allzero equalany suballtoprop propbigsplit equalityexceptpro postivevw random, over(smallalldm, relabel(1 "Fair" 2 "Mixed" 3 "Unfair"))  
	bar(1,color(gs0)) bar(2, color(gs3)) bar(3, color(gs5)) bar(4, color(gs7)) bar(5, color(gs9)) bar(6, color(gs11)) bar(7, color(gs13))
	stack name(bar1, replace) ytitle("% subject-periods") 
	legend(label(1 "None") label(2 "Equal") label(3 "Prop Only") label(4 "Prop and Largest") label(5 "Equal, more to Prop") label(6 "proportional, except prop") label(7 "Random")) scheme(s1mono);

# delimit cr



/***** FIGURE 2b in the Main Text  *****/

preserve 
drop if random==1
gen sizematters=propbigsplit+postivevw 
gen propbutnosize=suballtoprop+equalityexceptpro
gen nonorequal=allzero+equalany
# delimit ;
graph bar sizematters propbutnosize nonorequal, over(smallalldm, relabel(1 "Fair" 2 "Mixed" 3 "Unfair")) stack name(bar1, replace) ytitle("% subject-periods") 
bar(1,color(gs0)) bar(6, color(gs3)) bar(12, color(gs5)) 
	legend(label(1 "Uses sizes in some way") label(2 "Uses info about Proposer" "but no size info") label(3 "Neither size nor proposer info")) scheme(s1mono);

# delimit cr
restore 


/**********************************/
