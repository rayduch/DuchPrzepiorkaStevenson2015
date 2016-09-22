
 set more off 
  
  cd "~/Dropbox/Attribution Analysis Feb 2012/"

    use attribution.dta, clear


gen praldm = allocdm/25 if session==1 | session==4
replace praldm = allocdm/23 if session==2 | session==3
gen trt = 1 if treat1 & treat3
replace trt = 2 if treat1 & !treat3
replace trt = 3 if !treat1 & !treat3
lab def trt 1 "Distribution and proposer known" ///
2 "Distribution known, proposer unknown" ///
3 "Distribution and proposer unknown"
lab val trt trt

hist praldm if subject==1, d frac width(.03) by(trt, total ix iy ///
note("") graphregion(color(white) margin(zero))) color(gs8) ///
ytitle(Relative frequency) ylabel(0(.1).6, angle(horizontal)) ///
xtitle("Fraction kept by DMs") xlabel(0(.1)1)

