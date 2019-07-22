/* This project is about estimating the incremental cost of breast of cancers cases at their end of life*/
/*Declare libraries*/

libname disk21 'D:\Adam\eof-cervix';
libname disk22 'D:\Adam\eof-breast';

/* Select  breast cancer cases that died */
data disk22.eof_br_;
set disk4.br_enrol_025;
if dod_mcr~=. and dod_cl~=.;
run;
proc sort data=disk22.eof_br_ nodupkeys;
by dcn;
run;

/* Create Survival time and Code Dx stages */
data disk22.eof_br_010;
set disk22.eof_br_;
if dod_mcr <'31dec2011'd then
    dxdoddiff = floor(yrdif(dxdate,dod_mcr,'act/365')*12); 
else if dod_cl~=. then 
    dxdoddiff = floor(yrdif(dxdate,dod_cl,'act/365')*12);
    if n_3020 = 0 then dxstage = 1;
    else if n_3020 = 1 then dxstage = 2;
    else if n_3020 in (2,3,4,5) then dxstage = 3;
    else if n_3020 = 7 then dxstage = 4;
    else if n_3020 = 8 then dxstage = 5;
    else dxstage = 6;
run;

/* Check the Distribution of the cases by age group */
proc univariate data=disk22.eof_br_010;
class dxstage;
var dxdoddiff;
HISTOGRAM / NORMAL (COLOR=RED W=5) NROWS=2; 
run;

/* Run continuous enrollment */
%include 'D:\Adam\breast_cancer\eof_ce.sas';
%selectparticipants(data1=disk22.eof_br_010 /*Cases Characteristics from Cancer Registery*/,
                    data2=disk1.br_all_cl_030 /*Claims data*/, 
                    breakcutoff=15 /*Number of days of allowed gap*/, 
                    enrollrequired=6*30  /* Months*Days 'Number of months of continuous enrollment' */, 
                    k=disk22 /*  */);

/* Check the characteritics of those who died but not continuously enrolled */
proc sql;
create table eof_br_non_ce as
select distinct a.* from disk22.eof_br_010 a
left join disk22._mcr_final_6_eof b
on substr(a.dcn,1,8)=substr(b.dcn,1,8)
where b.dcn is null;
quit;

proc univariate data=eof_br_non_ce;
class dxstage;
var dxdoddiff;
HISTOGRAM / NORMAL (COLOR=RED W=5) NROWS=2; 
run;


/* Compute the comorbidity index */ 
%include 'D:\Adam\breast_cancer\analysis_macro_eof.sas';
%analysiscaseeof(data1=disk22._cl_mcr_final_6_eof_icn, group=6, cancer=eof_br,med=all,disk=disk22);

/* Run descriptive statistics*/
%include 'd:\adam\breast_cancer\desc_stat.sas';
%descrCase(data1=disk22._mcr_final_6_eof,data2=disk22._eof_br_6_all);
/* Prepare Breast Cancers' cases */
proc sql;
create table disk22.br_eof_case as
    select a.* from disk11.mcbrcse_030 a
        inner join disk22._mcr_final_6_eof b
            on substr(a.dcn_cse,1,8)=substr(b.dcn,1,8);
%runquit;

/* Prepare the case-control matching data */

proc sql;
create table disk22.brcasectrleof( bufsize=65536 )
  (
   dcn char(8) format=$8. informat=$8. label='DCN',
   PAYCOUNTY char(3) format=$3. informat=$3. label='PAYCOUNTY',
   MECODE char(2) format=$2. informat=$2. label='MECODE',
   FINANCIALGRANTIND char(1) format=$1. informat=$1. label='FINANCIALGRANTIND',
   DAYSPECIFICIND char(1) format=$1. informat=$1. label='DAYSPECIFICIND',
   ELIGIBILITYSTARTDATE num format=DATETIME20. informat=DATETIME20.
    label='ELIGIBILITYSTARTDATE',
   ELIGIBILITYSTOPDATE num format=DATETIME20. informat=DATETIME20. label='ELIGIBILITYSTOPDATE',
   dxdate num format=DATE9.,
   dcn_cse char(12),
   dod num format=DATE9.,
   start num,
   stop num,
   dod_cl num,
   _prev_stop num format=DATE9.,
   enroll_cumsum num,
   gap_count num,
   startingap num,
   startgaplength num,
   enrollbreak num,
   keeper num,
   r num
  );
%runquit;
*****run the control matching algorithm;
%include 'D:\Adam\BREAsT_CANCER\con_enrol_ctrl.sas';
%macro dolist;
%do i= 1 %to 4;
%matchCtrlCase(data1=disk22.br_eof_case/*case data*/, 
                data2=disk22.eof_ctrl_br/*control data*/,
                data3=disk22.brcasectrleof/*case&control data*/,
                breakcutoff=15, 
                enrollrequired=6*30,yrdiff=5);
%end;
%mend;
%dolist();
proc sort data=disk22.brcasectrleof;
by dcn_cse;
%runquit;

/* Prepare the data for cervical cancer patient */
/* Select Cases that die */
data disk21.eof_cr_;
set disk6.cr_enrol_030;
if dod_mcr~=. and dod_cl~=.;
run;

/* Compute the survival time */
data disk21.eof_cr_010;
set disk21.eof_cr_;
if dod_mcr <'31dec2011'd then
    dxdoddiff = floor(yrdif(dxdate,dod_mcr,'act/365')*12); 
else if dod_cl~=. then 
    dxdoddiff = floor(yrdif(dxdate,dod_cl,'act/365')*12);
run;

proc univariate data=disk21.eof_cr_010;
var dxdoddiff;
HISTOGRAM / NORMAL (COLOR=RED W=5) NROWS=2; 
run;


/* Compute the Continuous enrollment */
%include 'D:\Adam\breast_cancer\eof_ce.sas';
%selectparticipants(data1=disk21.eof_cr_010,data2=disk6.cr_all_cl_040,breakcutoff=15, enrollrequired=6*30,disk=disk21);
/* Compute the comorbidity index */ 
%include 'D:\Adam\breast_cancer\analysis_macro_eof.sas';
%analysiscaseeof(data1=disk21._cl_mcr_final_6_eof_icn, group=6, cancer=eof_cr,med=all,disk=disk21);
/* Descriptive statistics */
%include 'd:\adam\breast_cancer\desc_stat.sas';
%descrCase(data1=disk21._mcr_final_6_eof,data2=disk21._eof_cr_6_all);

/* prpepare eof case group */
proc sql;
create table disk21.cr_eof_case as
    select distinct a.* from disk15.mcrcse_020 a
        inner join disk21._mcr_final_6_eof b
            on substr(a.dcn_cse,1,8)=substr(b.dcn,1,8);
%runquit;

*****prepare the case-control matching data;
*****check whether non-disabled breast cancer patients are matched;
proc sql;
create table disk21.crcasectrleof( bufsize=65536 )
  (
   DCN char(8) format=$8. informat=$8. label='DCN',
   PAYCOUNTY char(3) format=$3. informat=$3. label='PAYCOUNTY',
   MECODE char(2) format=$2. informat=$2. label='MECODE',
   FINANCIALGRANTIND char(1) format=$1. informat=$1. label='FINANCIALGRANTIND',
   DAYSPECIFICIND char(1) format=$1. informat=$1. label='DAYSPECIFICIND',
   ELIGIBILITYSTARTDATE num format=DATETIME20. informat=DATETIME20.
    label='ELIGIBILITYSTARTDATE',
   ELIGIBILITYSTOPDATE num format=DATETIME20. informat=DATETIME20. label='ELIGIBILITYSTOPDATE',
   dxdate num format=DATE9.,
   dcn_cse char(12),
   dod num format=DATE9.,
   start num,
   stop num,
   dod_cl num,
   _prev_stop num format=DATE9.,
   enroll_cumsum num,
   gap_count num,
   startingap num,
   startgaplength num,
   enrollbreak num,
   keeper num,
   r num
  );
%runquit;
*****run the control matching algorithm;
%macro dolist();
%do i= 1 %to 4; %input cl&i._lline;
%include 'D:\Adam\BREAsT_CANCER\con_enrol_ctrl.sas';
%matchCtrlCase(data1=disk21.cr_eof_case/*case data*/, 
                data2=disk21.eof_ctrl_cr/*control data*/,
                data3=disk21.crcasectrleof/*case&control data*/,
                breakcutoff=15, 
                enrollrequired=6*30,yrdiff=5);
%end;
%mend  dolist;
%dolist;
*****frequency analysis for case control match for cervical cancer;
proc sql;
create table crcountcasectrl as
select count(distinct dcn) as n_ctrl,
count(distinct dcn_cse) as n_case from disk21.crcasectrleof;
create table crctrlpercase as
select distinct dcn_cse, count(distinct dcn) as n_ctrl from disk21.crcasectrleof
group by dcn_cse; 
quit;

proc freq data=crctrlpercase;
table n_ctrl;
run;
*****frequency analysis for case control match for breast cancer;
proc sql;
create table brcountcasectrl as
select count(distinct dcn) as n_ctrl,
count(distinct dcn_cse) as n_case from disk22.brcasectrleof;
create table brctrlpercase as
select distinct dcn_cse, count(distinct dcn) as n_ctrl from disk22.brcasectrleof
group by dcn_cse; 
quit;
proc freq data=brctrlpercase;
table n_ctrl;
run;

/* conduct descriptive analysis on the case and control the case of breast cancer */
%include 'd:\adam\breast_cancer\analysis_macro.sas';
%analysisctrleof(data1=disk22.brcasectrleof,data2=disk14.ctrl_cls_030, group=6, cancer=br_eof,med=ctrl,disk=disk22);
%include 'D:\adam\breast_cancer\desc_stat_eof.sas';
%DataRegAnacase(data1=disk22._mcr_final_6_eof /*all cases*/,data2=disk22._eof_br_6_all /*claims cases*/, 
data3=disk22.brcasectrleof /*control-case match*/, 
data4=disk11.mcbrcse_020 /*all cases regardless of the window*/);


/* prepare the control data for the regression */
%include 'D:\adam\breast_cancer\desc_stat_eof.sas';
%DataRegAnaCtrl(data1=disk22.eof_ctrl_br,
data2= disk22._br_eof_6_ctrl,
data3=disk11.mcbrcse_020,
data4=disk22.brcasectrleof);

%include 'D:\adam\breast_cancer\desc_stat_eof.sas';
%appendcasectrl(data1=disk22._br_eof_6_ctrl_cost /*control groupdata*/, 
data2=_mcr_final_6_eof_cost /*case group data*/, outds=eof_br_twopm_6);
%include 'D:\adam\breast_cancer\desc_stat_eof.sas';
%casectrldescr(data1=eof_br_twopm_6,grp=0);
%casectrldescr(data1=eof_br_twopm_6,grp=1);
%include 'D:\adam\breast_cancer\desc_stat_eof.sas';
%casectrlcost(data1=eof_br_twopm_6,grp=0);
%casectrlcost(data1=eof_br_twopm_6,grp=1);

proc export data= eof_br_twopm_ttl outfile='D:\Adam\Stata\eof\eof_br_twopm_ttl.dta'; %runquit;
proc export data=eof_br_twopm_6 outfile='D:\Adam\Stata\eof\eof_br_twopm_6.dta'; %runquit;

****export the data into stata:;
proc sql;
create table eof_br_twopm_ttl as
select distinct dcn, dcn_cse,grp,cci,age,dxstage,race, sum(sum_amt) as sum_amt
from eof_br_twopm_6
group by grp, dcn order grp, dcn_cse;
quit;
proc univariate data=eof_br_twopm_6;
class grp claim;
var sum_amt;
HISTOGRAM / NORMAL (COLOR=RED W=5) NROWS=2; 
run;
******conduct descriptive analysis on the case and control the case of cervical cancer;
%include 'd:\adam\breast_cancer\analysis_macro.sas';
%analysisctrleof(data1=disk21.crcasectrleof,data2=disk14.ctrl_cls_030, group=6, cancer=cr_eof,med=ctrl,disk=disk21);
****Compute breast cancer time to death from dx date;
proc format;
 value dxstage
    1 = 'in situ'
    2 = 'local'
    3 = 'regional'
    4 = 'distant'
    5 = 'not applicable'
    6 = 'unstaged';
run;

proc sql;
create table MonthToDeathALL as
select a.*,b.dxstage from disk22.brcasectrleof a
inner join disk22.eof_case b
on substr(left(a.dcn_cse),1,8)=substr(left(b.dcn_cse),1,8);
run;quit;

proc format;
value MonthToDeath
    1 ='One Year or Less'
    2 ='Two years'
    3= 'More than 2 years';
run;

proc sql;
select count(distinct dcn) from disk22.brcasectrleof; 
select count(distinct dcn_cse) from disk22.brcasectrleof;
run;quit;

data MonthToDeathALL;
set MonthToDeathALL;
NumMonthToDeathCtrl=intck('month',dxdate,dod) ;
NumMonthToDeathCase=intck('month',dxdate,dod_cl) ;
if NumMonthToDeathCase<=12 then MonthToDeathCase=1;
else if  NumMonthToDeathCase>12 and NumMonthToDeathCase<=24 then MonthToDeathCase=2;
else MonthToDeathCase=3;

if NumMonthToDeathCtrl<=12 then MonthToDeathCtrl=1;
else if  NumMonthToDeathCtrl>12 and NumMonthToDeathCtrl<=24 then MonthToDeathCtrl=2;
else MonthToDeathCtrl=3;

run;quit;


proc univariate data=MonthToDeathALL;
*class grp claim;
var NumMonthToDeathCtrl;
HISTOGRAM / NORMAL (COLOR=RED W=5) NROWS=2; 
run;


proc freq data=MonthToDeathALL;
table (MonthToDeathCase MonthToDeathCtrl)*dxstage;
format MonthToDeathCtrl MonthToDeath. dxstage dxstage.;
run;

proc sql;
select dcn,dod,dxdate from MonthToDeathALL where NumMonthToDeathCtrl<0; run;quit;
urvival time and Code Dx stages */
data disk22.eof_br_010;
set disk22.eof_br_;
if dod_mcr <'31dec2011'd then
    dxdoddiff = floor(yrdif(dxdate,dod_mcr,'act/365')*12); 
else if dod_cl~=. then 
    dxdoddiff = floor(yrdif(dxdate,dod_cl,'act/365')*12);
    if n_3020 = 0 then dxstage = 1;
    else if n_3020 = 1 then dxstage = 2;
    else if n_3020 in (2,3,4,5) then dxstage = 3;
    else if n_3020 = 7 then dxstage = 4;
    else if n_3020 = 8 then dxstage = 5;
    else dxstage = 6;
run;

/* Check the Distribution of the cases by age group */
proc univariate data=disk22.eof_br_010;
class dxstage;
var dxdoddiff;
HISTOGRAM / NORMAL (COLOR=RED W=5) NROWS=2; 
run;

/* Run continuous enrollment */
%include 'D:\Adam\breast_cancer\eof_ce.sas';
%selectparticipants(data1=disk22.eof_br_010 /*Cases Characteristics from Cancer Registery*/,
                    data2=disk1.br_all_cl_030 /*Claims data*/, 
                    breakcutoff=15 /*Number of days of allowed gap*/, 
                    enrollrequired=6*30  /* Months*Days 'Number of months of continuous enrollment' */, 
                    k=disk22 /*  */);

/* Check the characteritics of those who died but not continuously enrolled */
proc sql;
create table eof_br_non_ce as
select distinct a.* from disk22.eof_br_010 a
left join disk22._mcr_final_6_eof b
on substr(a.dcn,1,8)=substr(b.dcn,1,8)
where b.dcn is null;
quit;

proc univariate data=eof_br_non_ce;
class dxstage;
var dxdoddiff;
HISTOGRAM / NORMAL (COLOR=RED W=5) NROWS=2; 
run;


/* Compute the comorbidity index */ 
%include 'D:\Adam\breast_cancer\analysis_macro_eof.sas';
%analysiscaseeof(data1=disk22._cl_mcr_final_6_eof_icn, group=6, cancer=eof_br,med=all,disk=disk22);

/* Run descriptive statistics*/
%include 'd:\adam\breast_cancer\desc_stat.sas';
%descrCase(data1=disk22._mcr_final_6_eof,data2=disk22._eof_br_6_all);
/* Prepare Breast Cancers' cases */
proc sql;
create table disk22.br_eof_case as
    select a.* from disk11.mcbrcse_030 a
        inner join disk22._mcr_final_6_eof b
            on substr(a.dcn_cse,1,8)=substr(b.dcn,1,8);
%runquit;

/* Prepare the case-control matching data */

proc sql;
create table disk22.brcasectrleof( bufsize=65536 )
  (
   dcn char(8) format=$8. informat=$8. label='DCN',
   PAYCOUNTY char(3) format=$3. informat=$3. label='PAYCOUNTY',
   MECODE char(2) format=$2. informat=$2. label='MECODE',
   FINANCIALGRANTIND char(1) format=$1. informat=$1. label='FINANCIALGRANTIND',
   DAYSPECIFICIND char(1) format=$1. informat=$1. label='DAYSPECIFICIND',
   ELIGIBILITYSTARTDATE num format=DATETIME20. informat=DATETIME20.
    label='ELIGIBILITYSTARTDATE',
   ELIGIBILITYSTOPDATE num format=DATETIME20. informat=DATETIME20. label='ELIGIBILITYSTOPDATE',
   dxdate num format=DATE9.,
   dcn_cse char(12),
   dod num format=DATE9.,
   start num,
   stop num,
   dod_cl num,
   _prev_stop num format=DATE9.,
   enroll_cumsum num,
   gap_count num,
   startingap num,
   startgaplength num,
   enrollbreak num,
   keeper num,
   r num
  );
%runquit;
*****run the control matching algorithm;
%include 'D:\Adam\BREAsT_CANCER\con_enrol_ctrl.sas';
%macro dolist;
%do i= 1 %to 4;
%matchCtrlCase(data1=disk22.br_eof_case/*case data*/, 
                data2=disk22.eof_ctrl_br/*control data*/,
                data3=disk22.brcasectrleof/*case&control data*/,
                breakcutoff=15, 
                enrollrequired=6*30,yrdiff=5);
%end;
%mend;
%dolist();
proc sort data=disk22.brcasectrleof;
by dcn_cse;
%runquit;

/* Prepare the data for cervical cancer patient */
/* Select Cases that die */
data disk21.eof_cr_;
set disk6.cr_enrol_030;
if dod_mcr~=. and dod_cl~=.;
run;

/* Compute the survival time */
data disk21.eof_cr_010;
set disk21.eof_cr_;
if dod_mcr <'31dec2011'd then
    dxdoddiff = floor(yrdif(dxdate,dod_mcr,'act/365')*12); 
else if dod_cl~=. then 
    dxdoddiff = floor(yrdif(dxdate,dod_cl,'act/365')*12);
run;

proc univariate data=disk21.eof_cr_010;
var dxdoddiff;
HISTOGRAM / NORMAL (COLOR=RED W=5) NROWS=2; 
run;


/* Compute the Continuous enrollment */
%include 'D:\Adam\breast_cancer\eof_ce.sas';
%selectparticipants(data1=disk21.eof_cr_010,data2=disk6.cr_all_cl_040,breakcutoff=15, enrollrequired=6*30,disk=disk21);
/* Compute the comorbidity index */ 
%include 'D:\Adam\breast_cancer\analysis_macro_eof.sas';
%analysiscaseeof(data1=disk21._cl_mcr_final_6_eof_icn, group=6, cancer=eof_cr,med=all,disk=disk21);
/* Descriptive statistics */
%include 'd:\adam\breast_cancer\desc_stat.sas';
%descrCase(data1=disk21._mcr_final_6_eof,data2=disk21._eof_cr_6_all);

/* prpepare eof case group */
proc sql;
create table disk21.cr_eof_case as
    select distinct a.* from disk15.mcrcse_020 a
        inner join disk21._mcr_final_6_eof b
            on substr(a.dcn_cse,1,8)=substr(b.dcn,1,8);
%runquit;

*****prepare the case-control matching data;
*****check whether non-disabled breast cancer patients are matched;
proc sql;
create table disk21.crcasectrleof( bufsize=65536 )
  (
   DCN char(8) format=$8. informat=$8. label='DCN',
   PAYCOUNTY char(3) format=$3. informat=$3. label='PAYCOUNTY',
   MECODE char(2) format=$2. informat=$2. label='MECODE',
   FINANCIALGRANTIND char(1) format=$1. informat=$1. label='FINANCIALGRANTIND',
   DAYSPECIFICIND char(1) format=$1. informat=$1. label='DAYSPECIFICIND',
   ELIGIBILITYSTARTDATE num format=DATETIME20. informat=DATETIME20.
    label='ELIGIBILITYSTARTDATE',
   ELIGIBILITYSTOPDATE num format=DATETIME20. informat=DATETIME20. label='ELIGIBILITYSTOPDATE',
   dxdate num format=DATE9.,
   dcn_cse char(12),
   dod num format=DATE9.,
   start num,
   stop num,
   dod_cl num,
   _prev_stop num format=DATE9.,
   enroll_cumsum num,
   gap_count num,
   startingap num,
   startgaplength num,
   enrollbreak num,
   keeper num,
   r num
  );
%runquit;
*****run the control matching algorithm;
%macro dolist();
%do i= 1 %to 4; %input cl&i._lline;
%include 'D:\Adam\BREAsT_CANCER\con_enrol_ctrl.sas';
%matchCtrlCase(data1=disk21.cr_eof_case/*case data*/, 
                data2=disk21.eof_ctrl_cr/*control data*/,
                data3=disk21.crcasectrleof/*case&control data*/,
                breakcutoff=15, 
                enrollrequired=6*30,yrdiff=5);
%end;
%mend  dolist;
%dolist;
*****frequency analysis for case control match for cervical cancer;
proc sql;
create table crcountcasectrl as
select count(distinct dcn) as n_ctrl,
count(distinct dcn_cse) as n_case from disk21.crcasectrleof;
create table crctrlpercase as
select distinct dcn_cse, count(distinct dcn) as n_ctrl from disk21.crcasectrleof
group by dcn_cse; 
quit;

proc freq data=crctrlpercase;
table n_ctrl;
run;
*****frequency analysis for case control match for breast cancer;
proc sql;
create table brcountcasectrl as
select count(distinct dcn) as n_ctrl,
count(distinct dcn_cse) as n_case from disk22.brcasectrleof;
create table brctrlpercase as
select distinct dcn_cse, count(distinct dcn) as n_ctrl from disk22.brcasectrleof
group by dcn_cse; 
quit;
proc freq data=brctrlpercase;
table n_ctrl;
run;

/* conduct descriptive analysis on the case and control the case of breast cancer */
%include 'd:\adam\breast_cancer\analysis_macro.sas';
%analysisctrleof(data1=disk22.brcasectrleof,data2=disk14.ctrl_cls_030, group=6, cancer=br_eof,med=ctrl,disk=disk22);
%include 'D:\adam\breast_cancer\desc_stat_eof.sas';
%DataRegAnacase(data1=disk22._mcr_final_6_eof /*all cases*/,data2=disk22._eof_br_6_all /*claims cases*/, 
data3=disk22.brcasectrleof /*control-case match*/, 
data4=disk11.mcbrcse_020 /*all cases regardless of the window*/);


/* prepare the control data for the regression */
%include 'D:\adam\breast_cancer\desc_stat_eof.sas';
%DataRegAnaCtrl(data1=disk22.eof_ctrl_br,
data2= disk22._br_eof_6_ctrl,
data3=disk11.mcbrcse_020,
data4=disk22.brcasectrleof);

%include 'D:\adam\breast_cancer\desc_stat_eof.sas';
%appendcasectrl(data1=disk22._br_eof_6_ctrl_cost /*control groupdata*/, 
data2=_mcr_final_6_eof_cost /*case group data*/, outds=eof_br_twopm_6);
%include 'D:\adam\breast_cancer\desc_stat_eof.sas';
%casectrldescr(data1=eof_br_twopm_6,grp=0);
%casectrldescr(data1=eof_br_twopm_6,grp=1);
%include 'D:\adam\breast_cancer\desc_stat_eof.sas';
%casectrlcost(data1=eof_br_twopm_6,grp=0);
%casectrlcost(data1=eof_br_twopm_6,grp=1);

proc export data= eof_br_twopm_ttl outfile='D:\Adam\Stata\eof\eof_br_twopm_ttl.dta'; %runquit;
proc export data=eof_br_twopm_6 outfile='D:\Adam\Stata\eof\eof_br_twopm_6.dta'; %runquit;

****export the data into stata:;
proc sql;
create table eof_br_twopm_ttl as
select distinct dcn, dcn_cse,grp,cci,age,dxstage,race, sum(sum_amt) as sum_amt
from eof_br_twopm_6
group by grp, dcn order grp, dcn_cse;
quit;
proc univariate data=eof_br_twopm_6;
class grp claim;
var sum_amt;
HISTOGRAM / NORMAL (COLOR=RED W=5) NROWS=2; 
run;
******conduct descriptive analysis on the case and control the case of cervical cancer;
%include 'd:\adam\breast_cancer\analysis_macro.sas';
%analysisctrleof(data1=disk21.crcasectrleof,data2=disk14.ctrl_cls_030, group=6, cancer=cr_eof,med=ctrl,disk=disk21);
****Compute breast cancer time to death from dx date;
proc format;
 value dxstage
    1 = 'in situ'
    2 = 'local'
    3 = 'regional'
    4 = 'distant'
    5 = 'not applicable'
    6 = 'unstaged';
run;

proc sql;
create table MonthToDeathALL as
select a.*,b.dxstage from disk22.brcasectrleof a
inner join disk22.eof_case b
on substr(left(a.dcn_cse),1,8)=substr(left(b.dcn_cse),1,8);
run;quit;

proc format;
value MonthToDeath
    1 ='One Year or Less'
    2 ='Two years'
    3= 'More than 2 years';
run;

proc sql;
select count(distinct dcn) from disk22.brcasectrleof; 
select count(distinct dcn_cse) from disk22.brcasectrleof;
run;quit;

data MonthToDeathALL;
set MonthToDeathALL;
NumMonthToDeathCtrl=intck('month',dxdate,dod) ;
NumMonthToDeathCase=intck('month',dxdate,dod_cl) ;
if NumMonthToDeathCase<=12 then MonthToDeathCase=1;
else if  NumMonthToDeathCase>12 and NumMonthToDeathCase<=24 then MonthToDeathCase=2;
else MonthToDeathCase=3;

if NumMonthToDeathCtrl<=12 then MonthToDeathCtrl=1;
else if  NumMonthToDeathCtrl>12 and NumMonthToDeathCtrl<=24 then MonthToDeathCtrl=2;
else MonthToDeathCtrl=3;

run;quit;


proc univariate data=MonthToDeathALL;
*class grp claim;
var NumMonthToDeathCtrl;
HISTOGRAM / NORMAL (COLOR=RED W=5) NROWS=2; 
run;


proc freq data=MonthToDeathALL;
table (MonthToDeathCase MonthToDeathCtrl)*dxstage;
format MonthToDeathCtrl MonthToDeath. dxstage dxstage.;
run;

proc sql;
select dcn,dod,dxdate from MonthToDeathALL where NumMonthToDeathCtrl<0; run;quit;
