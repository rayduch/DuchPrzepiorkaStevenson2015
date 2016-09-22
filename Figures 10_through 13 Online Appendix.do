
/* build an approriate compositional data set that gives the share of the 30 deduction points that
   go to each decisionmaker 
   
   One issue in treating these as compositional data is the treatment of zeros. In many cases, they give
   nothing to any decisionmaker and in others nothing to some decisionmakers.
   
   A problem with the zero inflated negative binomial is that it identifies not giving to a signle dm 
   rather than not giving at all.
   
  */
  
  
 /*** some utilties I need - program is continued below  run before I run the program ****/

capture program drop _simp
program define _simp, rclass
   version 6.0

   di _n "Simulating main parameters.  Please wait...."
   syntax [, B(string) V(string) Sims(int 1000) Genname(string) ANTIsim DRAWT]

   * GENERATE RANDOM NORMAL OR RANDOM T VARIABLES
   if `sims' > _N {                                 /* expand ds to fit sims*/
      di in y _n "Note: Clarify is expanding your dataset from " _N /*
         */ " observations to `sims'" _n "observations in order to " /*
         */ "accommodate the simulations.  This will append" _n "missing " /*
         */ "values to the bottom of your original dataset." _n
      qui set obs `sims'
   }            
   if "`antisim'"~="" {                             /* antithetical sims    */
      local top = int(`sims'/2 + .5)                /*   calculate boundary */
      local bot = `top' + 1                         /*   for top&bottom half*/
   }
   if "`drawt'" ~= "" {                             /* for drawing from T   */
      tempvar u tfactor                             /* rather than Normal   */
      qui g `u' = uniform() in 1/`sims'
      if "`antisim'"~="" { qui replace `u'=1-`u'[_n-`top'] in `bot'/`sims' }
      qui gen `tfactor' = sqrt(e(df_r)/invchi(e(df_r),`u')) in 1/`sims'
   }
   local numpVC = colsof(`v')
   local i 1
   while `i' <= `numpVC' {
      tempvar u c`i'
      qui g `u' = uniform() in 1/`sims'
      if "`antisim'"~="" { qui replace `u'=1-`u'[_n-`top'] in `bot'/`sims' }
      if "`drawt'" == "" { qui gen `c`i''= invnorm(`u') in 1/`sims' }
      else { qui gen `c`i'' = invnorm(`u')*`tfactor' in 1/`sims' }
      local cnames `cnames' `c`i''                  /* collect names of vars*/
      local newvars `newvars' `genname'`i'          /* collect names newvars*/
      local i = `i' + 1
   }

   * SIMULATE BETAS FROM NORMAL OR T DISTRIBUTION
   tempname A row
   _chol `v' `numpVC'                              /* Cholesky decomp of V */
   matrix `A' = r(chol)
   matrix colnames `A' = `cnames'                /* cols to `c1'..`c`numpVC'' */
   matrix colnames `A' = sameeq:                 /* Thx to Randy Stevenson */
   di "% of simulations completed: " _c        
   local i 1
   while `i' <= `numpVC' {
      di int(`i'*100/`numpVC') "% " _c             /* display progress     */
      matrix `row' = `A'[`i',1...]                  /* get i^th row of A    */
      tempvar b`i'                                  /* temporary variable   */
      matrix score `b`i'' = `row'                   /* c(NxK) * row(1xK)'   */
      qui replace `b`i'' = `b`i'' + `b'[1,`i']      /* add mean             */
      local i = `i' + 1
   }

   * SAVE AND LABEL THE PARAMETERS
   local namepVC : colnames(`v')                    /* all parameters in VC */
   local eqnames : coleq(`v')                       /* all equs             */
   tokenize "`eqnames'"                             /* check for distinct   */
   if "`1'" ~= "_" { local haseqnm 1 }              /*   equation names in  */
   else { local haseqnm 0 }                         /*   the var-cov matrix */
   local i 1
   while `i' <= `numpVC' {                          /* for each parameter:  */
      qui gen `genname'`i' = `b`i''                 /*   save sims to dset  */
      local pname : word `i' of `namepVC'           /*   fetch name of param*/
      * if has equation name, add eqname to label
      if `haseqnm' {                                
         local eqname : word `i' of `eqnames'       
         label var `genname'`i' "Simulated `eqname':`pname' parameter"
      }                              
      * otherwise use simple label w/o an eqname
      else { label var `genname'`i' "Simulated `pname' parameter" }
      local i = `i' + 1
   }
   order `newvars'
   di _n
   return local newvars `newvars'  /* names of newvars that were created */

end


************************** COMPUTATION UTILITIES ****************************

*! version 1.3  April 24, 1999  Michael Tomz
* Cholesky decomposition of an arbitrary matrix
* Input: V, the original matrix
*        k, the dimension of the matrix
* Output: chol, the lower-triangle cholesky of V

capture program drop _chol
program define _chol, rclass
   version 6.0
   args V k
   tempname A
   capt matrix `A' = cholesky(`V')             /* square root of VC matrix  */
   if _rc == 506 {                             /* If VC ~ pos definite, ... */
      tempname eye transf varianc vcd tmp
      mat `eye' = I(`k')                       /* identity matrix           */
      mat `transf' = J(1,`k',0)                /* initialize transf matrix  */
      local i 1
      while `i' <= `k' {
         scalar `varianc' = `V'[`i',`i']       /* variance of parameter `i' */
         if `varianc' ~= 0 {                   /* if has variance, add row  */
            mat `tmp' = `eye'[`i',1...]        /* of `eye' to transf matrix */
            mat `transf' = `transf' \ `tmp'
         }
         local i = `i' + 1
      }
      mat `transf' = `transf'[2...,1...]       /* drop 1st row of transf mat*/
      mat `vcd' = `transf'*`V'*`transf''       /* decomposed VC (no 0 rows) */
      mat `A' = cholesky(`vcd')                /* square root of decomp VC  */
      mat `A' = `transf''*`A'*`transf'         /* rebuild full sq root mat  */               
   }
   else if _rc { matrix `A' = cholesky(`V') }  /* redisplay error message   */
   return matrix chol `A'
end


/******** done with utilties ********/ 
 
 

 set more off 
 foreach jj of numlist 1 2 3 4 { 
  
  cd "~/Dropbox/Attribution Analysis Feb 2012/"

    use attribution.dta, clear
  
  egen overallsubid=group(session period sid)
	sort  session period sid vw
	egen DMid=seq(), from(1) to(5)
  
  /**** Note: these are ordered so that when the data set is rearranges varaibles marked "1" are for the smallest dm
        vars marked 2 are the next biggest, etc... up until those marked 5 are for the largest dm in the distribution ***/
		
	egen totalpunishment=total(ndp), by(overallsubid)
	egen marksubject=tag(overallsubid)
	gen insample=0
	replace insample=1  if  treat1==1& treat3==1&session~=1

	keep if insample==1
	
/* now make a composition for those people that did something */


gen percentpunish=ndp/totalpunishment  /* note this generated missing for those who allocated nothing */
replace percentpunish=.00001 if percentpunish==0

egen percentpunishtotal2=total(percentpunish), by(overallsubid)
gen percentpunish2=percentpunish/percentpunishtotal2

egen test=total(percentpunish2), by(overallsubid)

rename vw weight
drop vw*
drop admv
drop ndp
drop dp
drop marksubject percentpunish
reshape wide adm weight  percentpunish2  , i(overallsubid) j(DMid)

gen lrpun1_5=ln(percentpunish21/percentpunish25)
gen lrpun2_5=ln(percentpunish22/percentpunish25)
gen lrpun3_5=ln(percentpunish23/percentpunish25)
gen lrpun4_5=ln(percentpunish24/percentpunish25)

/* the variable names of these dummies refer to the order in the distribution and the distribution so the first number
   is the order in the distribution where 1 is smallest and 5 is biggest adn the second numnber is treatment 1 2 3 or 4, which are different
   distributions */
   

gen weight1_2=(weight1==2)      
gen weight1_8=(weight1==8)
gen weight1_11=(weight1==11)   
gen weight1_17=(weight1==17)   

gen weight2_6=(weight2==6)
gen weight2_11=(weight2==11) 
gen weight2_13=(weight2==13) 
gen weight2_19=(weight2==19)    

gen weight3_10=(weight3==10)
gen weight3_14=(weight3==14) 
gen weight3_17=(weight3==17) 
gen weight3_20=(weight3==20) 

gen weight4_29=(weight4==29) 
gen weight4_19=(weight4==19) 
gen weight4_21a=(weight4==21&treat2==3) 
gen weight4_21b=(weight4==21&treat2==4) 

gen weight5_53=(weight5==53)  
gen weight5_48=(weight5==48)  
gen weight5_38=(weight5==38)  
gen weight5_23=(weight5==23)  
  
   
egen subjectperm=group(subject session) 
	

local eq1vars "lrpun1_5 allocdm weight1_2 weight1_8 weight1_11  adm1 adm2 adm3 adm4  totalpunishment"
local eq2vars "lrpun2_5 allocdm weight2_6 weight2_11 weight2_13  adm1 adm2 adm3 adm4   totalpunishment"
local eq3vars "lrpun3_5 allocdm weight3_10 weight3_14 weight3_17 adm1 adm2 adm3 adm4  totalpunishment"
local eq4vars "lrpun4_5 allocdm weight4_29 weight4_19 weight4_21a  adm1 adm2 adm3 adm4  totalpunishment"



capture drop bbt*

# delimit ;
sureg (lrpun1_5 allocdm weight1_2 weight1_8 weight1_11  adm1 adm2 adm3 adm4  totalpunishment)
	(lrpun2_5 allocdm weight2_6 weight2_11 weight2_13  adm1 adm2 adm3 adm4   totalpunishment)
	(lrpun3_5 allocdm weight3_10 weight3_14 weight3_17 adm1 adm2 adm3 adm4  totalpunishment)
	(lrpun4_5 allocdm weight4_29 weight4_19 weight4_21a  adm1 adm2 adm3 adm4  totalpunishment) if totalpunishment>0;

# delimit cr

	/* now capture some results from when we ran clogit, that we will feed to the functions
        that will generate our simualtions */
	
	
	      tempname b V sig dfsig
	      matrix `b' = e(b)                            /* 1 x k vector         */
	      matrix `V' = e(V)                            /* k x k variance matrix*/
	      local N = `e(N)'                             /* save # observations  */
	
		capture drop bbt
	      _simp, b(`b') v(`V') s(500) g(bbt)		/* run the program to generate simulations */ 



		  
		  		/* generate five cases (different proposers) for each distribution */
	
	set more off	
	capture drop pid
	gen pid=.
	replace pid=1 in 1
	replace pid=2 in 2
	replace pid=3 in 3
	replace pid=4 in 4
	replace pid=5 in 5
		

	foreach gg of numlist 1 2 3 4 5 {
	
		set trace off 
		local allocdm 20
		
		if `jj'==1 {

		
			local weight1 2		/*these are just for labels below */
			local weight2 6
			local weight3 10
			local weight4 29
			local weight5 53
		
			local weight1_2 1
			local weight1_8 0
			local weight1_11 0
			
			local weight2_6 1
			local weight2_11 0
			local weight2_13 0
			
			local weight3_10 1 
			local weight3_14 0
			local weight3_17 0
			
			local weight4_29 1
			local weight4_19 0
			local weight4_21a 0
			
		}
		
		
		else if `jj'==2 {
			local weight1 8
			local weight2 11
			local weight3 14
			local weight4 19
			local weight5 48
			
			local weight1_2 0
			local weight1_8 1
			local weight1_11 0
			
			local weight2_6 0
			local weight2_11 1
			local weight2_13 0
			
			local weight3_10 0 
			local weight3_14 1
			local weight3_17 0
			
			local weight4_29 0
			local weight4_19 1
			local weight4_21a 0
			
			
		}
		
		else if `jj'==3 {
			local weight1 11
			local weight2 13
			local weight3 17
			local weight4 21
			local weight5 38
			
			local weight1_2 0
			local weight1_8 0
			local weight1_11 1
			
			local weight2_6 0
			local weight2_11 0
			local weight2_13 1
			
			local weight3_10 0 
			local weight3_14 0
			local weight3_17 1
			
			local weight4_29 0
			local weight4_19 0
			local weight4_21a 1
		}
		
		else if `jj'==4 {
			local weight1 17
			local weight2 19
			local weight3 20
			local weight4 21
			local weight5 23
			
			local weight1_2 0
			local weight1_8 0
			local weight1_11 0
			
			local weight2_6 0
			local weight2_11 0
			local weight2_13 0
			
			local weight3_10 0 
			local weight3_14 0
			local weight3_17 0
			
			local weight4_29 0
			local weight4_19 0
			local weight4_21a 0
		}	
		
		
		
		local adm1 0
		local adm2 0 
		local adm3 0
		local adm4 0
		local adm5 0
		
		local adm`gg' 1
		
		local totalpunishment 30
		
		
		local index1=0
	# delimit ;	
		foreach inc of numlist 0 10 20 30 {;
			local index1=`index1'+1;
			foreach labnum of numlist 1 2 3 4 5 6 7 8 9 10 {;
				local lab`labnum'=`labnum'+`inc';
			};
			
			if `index1'==1 {;
			capture drop pred_lrpun1_5;
			gen pred_lrpun1_5=exp(bbt`lab1'*`allocdm'
			+bbt`lab2'*`weight1_2'
			+bbt`lab3'*`weight1_8'
			+bbt`lab4'*`weight1_11'
			+bbt`lab5'*`adm1'
			+bbt`lab6'*`adm2'
			+bbt`lab7'*`adm3'
			+bbt`lab8'*`adm4'
			+bbt`lab9'*`totalpunishment'+bbt`lab10');
			};
			
			if `index1'==2 {;
			capture drop pred_lrpun2_5;
			gen pred_lrpun2_5=exp(bbt`lab1'*`allocdm'
			+bbt`lab2'*`weight2_6'
			+bbt`lab3'*`weight2_11'
			+bbt`lab4'*`weight2_13'
			+bbt`lab5'*`adm1'
			+bbt`lab6'*`adm2'
			+bbt`lab7'*`adm3'
			+bbt`lab8'*`adm4'
			+bbt`lab9'*`totalpunishment'+bbt`lab10');
			};
			
			if `index1'==3 {;
			capture drop pred_lrpun3_5;
			gen pred_lrpun3_5=exp(bbt`lab1'*`allocdm'
			+bbt`lab2'*`weight3_10'
			+bbt`lab3'*`weight3_14'
			+bbt`lab4'*`weight3_17'
			+bbt`lab5'*`adm1'
			+bbt`lab6'*`adm2'
			+bbt`lab7'*`adm3'
			+bbt`lab8'*`adm4'
			+bbt`lab9'*`totalpunishment'+bbt`lab10');
			};
			
			if `index1'==4 {;
			capture drop pred_lrpun4_5;
			gen pred_lrpun4_5=exp(bbt`lab1'*`allocdm'
			+bbt`lab2'*`weight4_29'
			+bbt`lab3'*`weight4_19'
			+bbt`lab4'*`weight4_21a'
			+bbt`lab5'*`adm1'
			+bbt`lab6'*`adm2'
			+bbt`lab7'*`adm3'
			+bbt`lab8'*`adm4'
			+bbt`lab9'*`totalpunishment'+bbt`lab10');
			};
			
		};
	# delimit cr
	
		capture drop realpred*
		capture drop sumLR
		gen sumLR=1+pred_lrpun1_5+pred_lrpun2_5+pred_lrpun3_5+pred_lrpun4_5
		gen realpred1=pred_lrpun1_5/sumLR
		gen realpred2=pred_lrpun2_5/sumLR
		gen realpred3=pred_lrpun3_5/sumLR
		gen realpred4=pred_lrpun4_5/sumLR
		gen realpred5=1-(realpred1+realpred2+realpred3+realpred4)
		

		
		capture drop high`gg' low`gg' mean`gg'
	gen high`gg'=.
	gen low`gg'=.
	gen mean`gg'=.

	foreach num of numlist 1 2 3 4 5 {
	
			local diffname "realpred`num'"

			_pctile `diffname', p(2.5, 97.5)
			replace high`gg'=r(r2) in `num'
			replace low`gg'=r(r1) in `num'					
			sum `diffname'
			replace mean`gg'=r(mean) in `num'
	}
	
	/* technically these data are not needed to be saved in this new version of the program, but I save them
	   touse later to make some special graphs */
	   
	preserve
		keep high* low* mean*
		gen dist=`jj'
		gen prop=`gg'
		drop if mean`gg'==.
		gen dm=_n
		sort dist prop dm 
		save CIdata_dist_`jj'_prop`gg', replace	
	restore 
	
	/* a little adjustment just for the graphs to make very small confidence regions able to be seen. */
	
	gen confdiff=high`gg'-low`gg'
	replace high`gg'=high`gg'+.01 if high`gg'<.5&confdiff<.02
	replace low`gg'=low`gg'-.01 if low`gg'>.5&confdiff<.02
 
	drop confdiff
	
 
	label define gname 1 "`weight1'" 2 "`weight2'" 3 "`weight3'" 4 "`weight4'" 5 "`weight5'", replace
	label values pid gname



	capture drop pid
	gen pid=.
	replace pid=`weight1' in 1
	replace pid=`weight2' in 2
	replace pid=`weight3' in 3
	replace pid=`weight4' in 4
	replace pid=`weight5' in 5
	
	
	/* i need a frequency measure that gives the numner of times in the selected sample that percentpunish2X- a certain number 
		lest take a look at what these frequencies are */
		
		tabstat percentpunish21, by(percentpunish21) stats(n)
		
		/* taking a look at this, we clearly need to do some grouping */
		
		set more off
		foreach hh of numlist 1 2 3 4 5 {
			capture drop perpun2`hh'recode
			gen perpun2`hh'recode=.
			replace  perpun2`hh'recode=0 if percentpunish2`hh'>=0&percentpunish2`hh'<=.0001
			local increment=0
			foreach nn of numlist 1/100 {
				local increment=`increment'+.01
				replace  perpun2`hh'recode=`increment' if percentpunish2`hh'>`increment'-.01&percentpunish2`hh'<=`increment'

			}
		}
	


foreach kkk of numlist 1 2 3 4 5 {
	capture drop freqneed_dm`kkk'_prop1
	egen freqneed_dm`kkk'_prop1=count(perpun21recode) if treat2==`jj'&adm1==1&totalpunishment>0  	/* the number of subjects who punished DM `kkk' */
																								/* when distribution was = 1 and DM1 was proposer */
	capture drop freqneed_dm`kkk'_prop2
	egen freqneed_dm`kkk'_prop2=count(perpun22recode) if treat2==`jj'&adm2==1&totalpunishment>0  	/* the number of subjects who punished DM `kkk' */
																								/* when distribution was = 1 and DM2 was proposer */
	capture drop freqneed_dm`kkk'_prop3
	egen freqneed_dm`kkk'_prop3=count(perpun23recode) if treat2==`jj'&adm3==1&totalpunishment>0  	/* the number of subjects who punished DM `kkk' */
																								/* when distribution was = 1 and DM3 was proposer */																					  
	capture drop freqneed_dm`kkk'_prop4
	egen freqneed_dm`kkk'_prop4=count(perpun24recode) if treat2==`jj'&adm4==1&totalpunishment>0  	/* the number of subjects who punished DM `kkk' */
																								/* when distribution was = 1 and DM4 was proposer */
	capture drop freqneed_dm`kkk'_prop5
	egen freqneed_dm`kkk'_prop5=count(perpun25recode) if treat2==`jj'&adm5==1&totalpunishment>0  	/* the number of subjects who punished DM `kkk' */
}																					 		 	/* when distribution was = 1 and DM5 was proposer */





/* now I need a variable that gives the number of subjects who punished a given amount for each case above, where we devide the possible 
   amounts (between 0 and 1) into 20 bins */


foreach kkk of numlist 1 2 3 4 5 {
	capture drop freqbin_dm`kkk'_prop`gg'
	gen freqbin_dm`kkk'_prop`gg'=.
	
	capture drop onecasefreq_dm`kkk'_prop`gg'
	gen onecasefreq_dm`kkk'_prop`gg'=.
	
	local inc1=0
		foreach nn of numlist 1/20 {
			local inc2=`inc1'+.05
			capture drop temp222		     
			egen temp222=count(percentpunish2`kkk') if treat2==`jj'&adm`gg'==1&totalpunishment>0&(percentpunish2`kkk'>`inc1'&percentpunish2`kkk'<`inc2')
			replace freqbin_dm`kkk'_prop`gg'=temp222 if treat2==`jj'&adm`gg'==1&totalpunishment>0&(percentpunish2`kkk'>`inc1'&percentpunish2`kkk'<`inc2') 
			capture drop temp333
			egen temp333=tag(freqbin_dm`kkk'_prop`gg') if treat2==`jj'&adm`gg'==1&totalpunishment>0&(percentpunish2`kkk'>`inc1'&percentpunish2`kkk'<`inc2')
			replace onecasefreq_dm`kkk'_prop`gg'=temp333 if treat2==`jj'&adm`gg'==1&totalpunishment>0&(percentpunish2`kkk'>`inc1'&percentpunish2`kkk'<`inc2')
			local inc1=`inc1'+.05				
		}
}


foreach kkk of numlist 1 2 3 4 5 {
	capture drop c_`kkk'
	capture drop g_`kkk'
	generate  c_`kkk'=round((freqbin_dm`kkk'_prop`gg'-0)/(freqneed_dm`kkk'_prop`gg'-0)*255)			/* this is creating the color */
	quiet levelsof c_`kkk', loc(cs)	
	local g`kkk'					 /* you are going to be creating a scatter of one point (expanded really big and made into a square and colored) 
								to get the heat map we want for one point, so we will build the graphics comands needed for building the 
								scatter point for each point and store it in a local g that we will use */
	foreach c of loc cs {	
	   local c1=120+round(`c'/2)
	   local c2=255-round(`c'/2)	
	   local m`kkk' mcolor("`c2' `c2' `c1'")	/* sets the RGB color using three numbers, 255 255 255 is black */
	   /* local g`kkk' `g`kkk''||scatter perpun2`kkk'recode weight`kkk' if c_`kkk'==`c'&treat2==`jj'&adm`gg'==1&totalpunishment>0, msymbol(S) msize(vhuge) `m' */
	 }	
	 
	 local g`kkk' `g`kkk''||sc perpun2`kkk'recode weight`kkk' if onecasefreq_dm`kkk'_prop`gg'==1&treat2==`jj'&adm`gg'==1&totalpunishment>0, ms(i) mlab(freqbin_dm`kkk'_prop`gg') mlabposition(3) mlabsize(1.5) mlabcolor(black)	/* this adds the labels in the squares */
	 local g`kkk' `g`kkk'' legend(off) ylabel(0 .1 .2 .3 .4 .5 .6 .7 .8 .9 1, notick)             
	 local g`kkk' `g`kkk'' 
 }
 
 
twoway `g1' || `g2' || `g3' || `g4' || `g5' xsize(4) ysize(4)  ///
		|| rbar high`gg' low`gg' pid, mwidth msize(1.5) color(gray) ///
		|| scatter  mean`gg' pid, msymbol(Oh) mcolor(black) ///
		legend(off) ///  
		ylabel(0 .2 .4 .6 .8 1, labsize(vsmall)) ///
		ytitle(" ")   ///
		xlabel(`weight1' `weight2' `weight3' `weight4' `weight5',  labsize(vsmall)) msize(small)  ///
		xtitle(" " "Voting Weight", size(vsmall)) ///
		plotregion(margin(r+5)) graphregion(fcolor(none)) ///
		/* ytitle(" " "Predicted Share of Punishment and"  "Actual share of Punishment" " ", size(small))*/ ///
		title("Proposer = DM with Weight `weight`gg''", size(vsmall)) ///
		name(gg`gg', replace) nodraw
}


	graph combine gg1 gg2 gg3 gg4 gg5, title("Distribution = {`weight1' , `weight2', `weight3', `weight4', `weight5'}", size(small)) saving(punishmentSUR`jj', replace) commonscheme scheme(s2mono)
	
	graph export punishmentSUR`jj'.png, replace
	
}		
		

		
		*graph combine punishment1.gph punishment2.gph punishment3.gph punishment4.gph, note("A case in which the DMs keeps 'allocdm' for themselves" "and a pubisher uses 'totalpunishment' of his 30 points", size (vsmall)) 
	
		
		
		
		
		/*

 display "`g'"
	
	
	
	twoway rbar high1 low1 pid, barwidth(2) fintensity(inten20) ///
		|| scatter  mean1 pid, legend(off)  ylabel(0 .2 .4 .6 .8 1, labsize(vsmall)) ytitle(" ")   ///
		xlabel(, valuelabel labsize(small)) msize(medium) ///
		xtitle(" " "Seat Weight", size(small)) ///
		ytitle(" " "Predicted Share of Punishment" and "Actual share of Punishment" " ", size(vsmall)) ///
		title("Propser = DM with Weight 6", size(vsmall))  ///
		|| scatter perpun21recode weight1 [w=countfreq1] if treat2==`jj'&adm2==1&onecase==1, msize(tiny) msymbol(Oh) mcolor(black)  /// 
		|| scatter perpun22recode weight2 [w=countfreq2] if treat2==`jj'&adm2==1&onecase==1, msize(tiny) msymbol(Oh) mcolor(black) ///
		|| scatter perpun23recode weight3 [w=countfreq3] if treat2==`jj'&adm2==1&onecase==1, msize(tiny) msymbol(Oh) mcolor(black) ///
		|| scatter perpun24recode weight4 [w=countfreq4] if treat2==`jj'&adm2==1&onecase==1, msize(tiny) msymbol(Oh) mcolor(black) ///
		|| scatter perpun25recode weight5 [w=countfreq5] if treat2==`jj'&adm2==1&onecase==1,  msize(tiny) msymbol(Oh) mcolor(black) name(yy, replace) legend(off) xlabel(2 6 10 29 53) 
		
		
		graph combine uu yy, col(1)
		
		
		
*/
