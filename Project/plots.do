/*******************************************************************************
  
  ECON281-Computational Final Project
  
  Plot real PCE and CPI of durable and services over time
  
*******************************************************************************/
 
 
********************************************************************************
* Setup
********************************************************************************
 
clear all
set more off
set matsize 8000
set maxvar 32000

cd "C:\Users\JungHyun\Documents\GitHub\281-Computational\Project"

cap log close
log using plots.log, replace

********************************************************************************
* Import and plot data
********************************************************************************

*** real PCE *** 

import excel "RPCE.xls", case(upper) cellrange(A13) clear

rename (A B C) (date rPCE_dur rPCE_serv)

gen date_monthly = mofd(date)
format date_monthly %tm
tsset date_monthly

tw tsline rPCE_dur rPCE_serv, xtitle("Time") ytitle("Index") ///
       legend(ring(0) pos(5) order(1 "Durable Goods" 2 "Services"))
	   
graph export "rPCE.eps", replace


*** CPI *** 

import excel "CPI.xls", case(upper) cellrange(A13) clear

rename (A B C) (date CPI_dur CPI_serv)
gen date_monthly = mofd(date)
format date_monthly %tm
tsset date_monthly

tw tsline CPI_dur CPI_serv, xtitle("Time") ytitle("% change from previous year") ///
       legend(ring(0) pos(1) order(1 "Durable Goods" 2 "Services"))
	   
graph export "CPI.eps", replace

********************************************************************************
* Plot data
********************************************************************************

log close