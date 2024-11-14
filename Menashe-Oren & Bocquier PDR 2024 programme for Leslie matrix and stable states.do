// Programme for reading in demographic data from csv files to populate matrix used to estimate population at stable state //

/*
Used in article:
	The potential of internal migration to shape rural and urban populations across Africa, Asia and Latin America 
	Menashe-Oren & Bocquier
	Population and Development Review 2024
	Special 50th Anniversary Collection, Looking Backward, Looking Forward: Celebrating 50 Years of Population and Development Review
*/


// csv files are named according to continent, period, sector and sex 
/*
*example of data in csv file
Age	 Population	Survival ratio	Fertility rate	Net migration rate
0	 125370	    0.959504382	    0	            -0.001030169
5	 125563	    0.988043488	    0	            -0.0003489
10	 117202	    0.988858852	    0	            -0.0008895
15	 98686	    0.987073997	    83.9782671	    -0.0070858
20	 79214	    0.985338504	    222.5235371	    -0.0092611
25	 72749	    0.983869297	    246.6928904	    -0.0040135
30	 56367	    0.981526369	    167.0663101	    -0.0019156
35	 46896	    0.977183464	    102.7605412	    -0.0012439
40	 43361	    0.97132354	    47.86816909	    -0.0010974
45	 39518	    0.960772213	    12.97516208	    -0.0007638
50	 34131	    0.944219711	    0	            -0.0007122
55	 28685	    0.915922814	    0	            -0.0005665
60	 23956	    0.870458468	    0	            -0.0006559
65	 18794	    0.81010166	    0	            -0.0006941
70	 13686	    0.721768961	    0	            -0.0006496
75	 8514	    0.599776515	    0	            -0.0006658
80	 6223	    0	            0	            -0.0007426

*/


*Set working directory
cd "C:\Users\XXX\rates csv"

clear
local file=1
foreach continent in Africa Asia LAC {
	foreach period in 1970_1984 1985_1999 2000_2014 {
	if "`continent'"=="Africa" & "`period'"=="1970_1984" local SRB=102.8/100
	if "`continent'"=="Asia" & "`period'"=="1970_1984" local SRB=105.9/100
	if "`continent'"=="LAC" & "`period'"=="1970_1984" local SRB=103.3/100
	if "`continent'"=="Africa" & "`period'"=="1985_1999" local SRB=102.9/100
	if "`continent'"=="Asia" & "`period'"=="1985_1999" local SRB=108.1/100
	if "`continent'"=="LAC" & "`period'"=="1985_1999" local SRB=103.8/100
	if "`continent'"=="Africa" & "`period'"=="2000_2014" local SRB=102.8/100
	if "`continent'"=="Asia" & "`period'"=="2000_2014" local SRB=109.7/100
	if "`continent'"=="LAC" & "`period'"=="2000_2014" local SRB=104.1/100
		foreach area in Rural Urban {
			foreach sex in Male Female {
				capture import delimited "`continent'_`period'_`area'_`sex'.csv", ///
						varnames(1) case(preserve)
				if _rc==4 | _rc==0 { 
				* Make the Leslie matrix
				quietly replace f=f/1000
				mkmat S, matrix(S)
				mkmat f, matrix(f)
				matrix Sx = J(_N, _N, 0)
				matrix Bx = J(_N, _N, 0)
				local n = _N - 1
				forval i = 1/`n' {
					matrix Sx[`i'+1, `i'] = S[`i',1]
					* CHECK Sx[1,1]:
					if "`sex'"=="Female" {
						matrix Bx[1,`i'] = (1/(1+`SRB')) * (S[1,1] * 5 / 2) * ///
										((f[`i'+1,1]*Sx[`i'+1,`i']) + f[`i',1])
					}
					if "`sex'"=="Male" {
						matrix Bx[1,`i'] = (`SRB'/(1+`SRB')) * (S[1,1] * 5 / 2) * ///
										((f[`i'+1,1]*Sx[`i'+1,`i']) + f[`i',1])
					}
				}
				matrix Sx[_N,_N]=0 // check assumption for last age group
				matrix rate_surv_fert = Sx + Bx
				* Add the net migration vector to the Leslie matrix
				mkmat net, matrix(net)
				matrix net = net*5
				matrix rate_surv_fert_nmig = rate_surv_fert + diag(net)
				matrix rate_surv_fert_nmig[_N,_N]= 0 
				* use "mata" functions to extract eigenvector
				mata: A = st_matrix("rate_surv_fert_nmig")
				mata: X = .
				mata: L = .
				mata: eigensystem(A, X, L)
				mata: index = (1, 1)
				mata: eigensystemselecti(A, index, X=., L=.)
				mata: L = Re(L)
				mata: st_matrix("EigenValue", L)
				mata: X = Re(X)
				mata: st_matrix("EigenVector", X)
				mata: st_matrix("SumEigenVector", colsum(X[.,1]))
				* Intrinsic growth annual rate assuming exponential growth
				matrix IntrinsicGrowth = (sqrt(EigenValue[1,1]^2) - 1) / 5  
				svmat IntrinsicGrowth, names(IntrGrowthNmig)
				* Intrinsic population structure
				matrix IntrinsicDistrib = EigenVector/SumEigenVector[1,1] 
				svmat IntrinsicDistrib, names(IntrDistNmig) 
				
				* Compute age % distribution of data
				mkmat N, matrix(N)
				mata: X = st_matrix("N")
				mata: st_matrix("total_pop", colsum(X[.,1]))
				matrix vector_pc_N=N/total_pop[1,1]
				svmat vector_pc_N, names(InitDistNmig)
				* Compute absolute difference b/w data & intrinsic distributions
				matrix diff_intr_N = vector_pc_N - IntrinsicDistrib
				mata: st_matrix("absval_intr_N",abs(st_matrix("diff_intr_N")))
				* Compute Sum of Absolute Deviation b/w data & intrinsic distributions
				mata: X = st_matrix("absval_intr_N")
				mata: st_matrix("SAD_intr_N", colsum(X[.,1]))
				svmat SAD_intr_N, names(IntrSADNmig)
				* Compute time to stable state
				scalar t = 0
				scalar SAD_intr_temp_N = .1
				matrix temp_matrix = rate_surv_fert_nmig*rate_surv_fert_nmig
				while (SAD_intr_temp_N>0.001) {
					matrix temp_N  = temp_matrix*N
					mata: X = st_matrix("temp_N")
					mata: st_matrix("total_pop", colsum(X[.,1]))
					matrix vector_pc_temp_N=temp_N/total_pop[1,1]
					mata: st_matrix("diff_intr_temp_N", st_matrix("vector_pc_temp_N") - st_matrix("IntrinsicDistrib"))
					mata: st_matrix("absval_intr_temp_N",abs(st_matrix("diff_intr_temp_N")))
					mata: X = st_matrix("absval_intr_temp_N")
					mata: st_matrix("SAD_intr_temp_N", colsum(X[.,1]))
					scalar SAD_intr_temp_N = SAD_intr_temp_N[1,1]
					if t==10 svmat vector_pc_temp_N, names(Dist50yearsNmig)
					scalar t = t+1
					matrix temp_matrix = temp_matrix*rate_surv_fert_nmig
					
				} // end while
				gen byte IntrStepsNmig=cond(_n==1,t,.) // #steps stable state SAD<1/1000 
				svmat temp_N , names(abs_N_stablemig) // to get the population size at stable state
				
				
				* Without migration
				mata: A = st_matrix("rate_surv_fert")
				mata: X = .
				mata: L = .
				mata: eigensystem(A, X, L)
				mata: index = (1, 1)
				mata: eigensystemselecti(A, index, X=., L=.)
				mata: L = Re(L)
				mata: st_matrix("EigenValue", L)
				mata: X = Re(X)
				mata: st_matrix("EigenVector", X)
				mata: st_matrix("SumEigenVector", colsum(X[.,1]))
				* Intrinsic growth annual rate assuming exponential growth
				matrix IntrinsicGrowth = (sqrt(EigenValue[1,1]^2) - 1) / 5  
				svmat IntrinsicGrowth, names(IntrGrowthNat)
				* Intrinsic population structure
				matrix IntrinsicDistrib = EigenVector/SumEigenVector[1,1] 
				svmat IntrinsicDistrib, names(IntrDistNat) 
				* Compute absolute difference b/w data & intrinsic distributions
				matrix diff_intr_N = vector_pc_N - IntrinsicDistrib
				mata: st_matrix("absval_intr_N",abs(st_matrix("diff_intr_N")))
				* Compute Sum of Absolute Deviation b/w data & intrinsic distributions
				mata: X = st_matrix("absval_intr_N")
				mata: st_matrix("SAD_intr_N", colsum(X[.,1]))
				svmat SAD_intr_N, names(IntrSADNat)
				* Compute time to stable state
				scalar t = 0
				scalar SAD_intr_temp_N = .1
				matrix temp_matrix = rate_surv_fert*rate_surv_fert
				while (SAD_intr_temp_N>0.001) {
					matrix temp_N  = temp_matrix*N
					mata: X = st_matrix("temp_N")
					mata: st_matrix("total_pop", colsum(X[.,1]))
					matrix vector_pc_temp_N=temp_N/total_pop[1,1]
					mata: st_matrix("diff_intr_temp_N", st_matrix("vector_pc_temp_N") - st_matrix("IntrinsicDistrib"))
					mata: st_matrix("absval_intr_temp_N",abs(st_matrix("diff_intr_temp_N")))
					mata: X = st_matrix("absval_intr_temp_N")
					mata: st_matrix("SAD_intr_temp_N", colsum(X[.,1]))
					scalar SAD_intr_temp_N = SAD_intr_temp_N[1,1]
					if t==10 svmat vector_pc_temp_N, names(Dist50yearsNat)
					scalar t = t+1
					matrix temp_matrix = temp_matrix*rate_surv_fert
				} // end while
				gen byte IntrStepsNat=cond(_n==1,t,.) // #steps stable state SAD<1/1000 
				svmat temp_N , names(abs_N_stableNat) // to get the population size at stable state

				* Compute SAD with and w/o migration
				gen diff_Nmig_Nat = sqrt((IntrDistNmig1 - IntrDistNat1)^2)
				egen IntrSADNmig_Nat=sum(diff_Nmig_Nat)
				replace IntrSADNmig_Nat=. if _n!=1
				drop diff_Nmig_Nat 
				gen c="`continent'"
				gen period=real(substr("`period'",1,4))
				gen a="`area'"
				gen s="`sex'"
				
				save "`continent'_`period'_`area'_`sex'.dta", replace
				if `file'==1 { 
					save LeslieStable.dta, replace 
					}
				else { 
					append using LeslieStable.dta 
					save LeslieStable.dta, replace
					}
				local file=`file'+1
				display "File:" `file'
				clear
				} // end _rc (no file)
			} // end sex
		} // end area
	} // end period
} // end continent

use LeslieStable
encode c, gen(continent)
drop c
encode a, gen(rur_urb)
drop a
encode s, gen(sex)
drop s

lab var N "Population size"
lab var S "Survival rate"
lab var f "Fertility rate"
lab var net "Net migration rate"
lab var IntrGrowthNmig1 "Intrinsic growth annual rate with migration"                
lab var IntrDistNmig1   "Intrinsic population structure with migration"                
lab var InitDistNmig1   "Initial population structure"                
lab var IntrSADNmig1  "Sum of Absolute Deviation b/w initial & intrinsic distributions with migration"                  
lab var Dist50yearsNmig1 "Population structure after 50years, with migration"                 
lab var IntrStepsNmig  "Number of steps (in 5years) to reach stable population, with migration"                 
lab var IntrGrowthNat1 "Intrinsic growth annual rate without migration (Natural)"         
lab var IntrDistNat1  "Intrinsic population structure without migration (Natural)"                  
lab var IntrSADNat1   "Sum of Absolute Deviation b/w initial & intrinsic distributions without migration (Natural)"                  
lab var Dist50yearsNat1  "Population structure after 50years, without migration (Natural)"                
lab var IntrStepsNat   "Number of steps (in 5years) to reach stable population, without migration"                  
lab var IntrSADNmig_Nat  "Sum of Absolute Deviation b/w intrinsic distributions with & without migration"
lab var abs_N_stableNat1 "Absolute population size at Stable state (Natural)"
lab var abs_N_stablemig1 "Absolute population size at Stable state, with migration"


save, replace


*******************************************************************************************





