libname esrd 'D:\建豐\新案\PreESRD';

/*********整理預測模型的語法*********/
/***Effect size***/
/*ESRD*/
%macro KFRE (name1, name2, name3, name4, name5, name6);
proc sort data=&name1.;
by _imp;
run;

proc phreg data=&name1.;
by _imp;
class &name6. gender_dw(ref="女");
model &name2.*&name3.(0)=&name5. age_entry gender_dw eGFR_CKDEPI PCRtotal_log/risklimits; 
ods output ParameterEstimates=&name4._KFRE_survival;	
run;
proc sort data=&name4._KFRE_survival; by parameter ClassVal0 _imp;run;
proc mianalyze data=&name4._KFRE_survival;
	by parameter ClassVal0;
	modeleffects estimate;
	stderr stderr;
	ods output ParameterEstimates =&name4._KFRE_survival_1;
run;
data &name4._KFRE_survival_1_HR;
	set &name4._KFRE_survival_1;
	Log_HR_comb=Estimate;
	HR_comb=exp(Estimate);
	HR_LCL_comb=exp(LCLMean);
	HR_UCL_comb=exp(UCLMean);
	keep Parameter ClassVal0 Log_HR_comb HR_comb HR_LCL_comb HR_UCL_comb
			Probt;
	rename Probt=HR_pval_comb;
run;

data &name4._KFRE_survival_1_HR; set &name4._KFRE_survival_1_HR;
Result=trim(left(round(HR_comb, 0.001)))||" ("||trim(left(round(HR_LCL_comb, 0.001)))||", "||trim(left(round(HR_UCL_comb, 0.001)))||")";
run;

proc print data=&name4._KFRE_survival_1_HR;
run;
%mend;
%KFRE(esrd.Ctrgroupfollowup_mi, dialysistime, event, dialysis, , ) /*KFRE*/
%KFRE(esrd.Ctrgroupfollowup_mi, dialysistime, event, dialysis, CCTR, ) /*KFRE+Continuous CTR*/
%KFRE(esrd.Ctrgroupfollowup_mi, dialysistime, event, dialysis, CTRgroup_New, CTRgroup_New (ref="1")) /*KFRE+Categorical CTR*/

/*CV/All-cause mortality*/
%macro RMPM (name1, name2, name3, name4, name5, name6);
proc phreg data=&name1.;
by _imp;
class &name6. gender_dw(ref="女") havedm(ref='0') havehtn(ref='0') Amemia(ref='0');
model &name2.*&name3.(0)=&name5. age_entry gender_dw havedm havehtn Amemia eGFR_CKDEPI/ risklimits; 
ods output ParameterEstimates=&name4._RMPM_survival;	
run;
proc sort data=&name4._RMPM_survival; by parameter ClassVal0 _imp;run;
proc mianalyze data=&name4._RMPM_survival;
	by parameter ClassVal0;
	modeleffects estimate;
	stderr stderr;
	ods output ParameterEstimates =&name4._RMPM_survival_1;
run;
data &name4._RMPM_survival_1_HR;
	set &name4._RMPM_survival_1;
	Log_HR_comb=Estimate;
	HR_comb=exp(Estimate);
	HR_LCL_comb=exp(LCLMean);
	HR_UCL_comb=exp(UCLMean);
	keep Parameter ClassVal0 Log_HR_comb HR_comb HR_LCL_comb HR_UCL_comb
			Probt;
	rename Probt=HR_pval_comb;
run;

data &name4._RMPM_survival_1_HR; set &name4._RMPM_survival_1_HR;
Result=trim(left(round(HR_comb, 0.001)))||" ("||trim(left(round(HR_LCL_comb, 0.001)))||", "||trim(left(round(HR_UCL_comb, 0.001)))||")";
run;

proc print data=&name4._RMPM_survival_1_HR;
run;
%mend;
%RMPM(esrd.Ctrgroupfollowup_mi, CVTime, CV_event, CV, , ) /*RMPM*/
%RMPM(esrd.Ctrgroupfollowup_mi, CVTime, CV_event, CV, CCTR, ) /*RMPM+Continuous CTR*/
%RMPM(esrd.Ctrgroupfollowup_mi, CVTime, CV_event, CV, CTRgroup_New , CTRgroup_New (ref="1")) /*RMPM+Categorical CTR*/
%RMPM(esrd.Ctrgroupfollowup_mi, deathtime, death, death, , ) /*RMPM*/
%RMPM(esrd.Ctrgroupfollowup_mi, deathtime, death, death, CCTR, ) /*RMPM+Continuous CTR*/
%RMPM(esrd.Ctrgroupfollowup_mi, deathtime, death, death, CTRgroup_New , CTRgroup_New (ref="1")) /*RMPM+Categorical CTR*/


/***Model Fitness (AIC, BIC)***/
/*ESRD*/
%macro KFRE(name1, name2, name3, name4, name5);
proc sort data=&name1.;
by _imp;
run;

proc phreg data=&name1.;
by _imp;
class &name5. gender_dw(ref="女");
model &name2.*&name3.(0)=&name4. age_entry gender_dw eGFR_CKDEPI PCRtotal_log/ risklimits;
ods output  FitStatistics=Dialysis_KFRE_Fit;	
run;

proc sql;
select mean(WithCovariates) as AIC
from Dialysis_KFRE_Fit
where Criterion='AIC';
quit;

proc sql;
select mean(WithCovariates) as BIC
from Dialysis_KFRE_Fit
where Criterion=';SBC ;';
quit;
%mend;
%KFRE(esrd.Ctrgroupfollowup_mi, dialysistime, event, , ) /*KFRE*/
%KFRE(esrd.Ctrgroupfollowup_mi, dialysistime, event, CCTR, ) /*KFRE+Continuous CTR*/
%KFRE(esrd.Ctrgroupfollowup_mi, dialysistime, event, CTRgroup_New, CTRgroup_New(ref="1")) /*KFRE+Categorical CTR*/

/*CV/All-cause Mortality*/
%macro Death (name1, name2, name3, name4, name5);
proc sort data=&name1.;
by _imp;
run;

proc phreg data=&name1.;
by _imp;
class &name5. gender_dw(ref="女") havedm(ref='0') havehtn(ref='0') Amemia(ref='0');
model &name2.*&name3.(0)=&name4. age_entry gender_dw havedm havehtn Amemia eGFR_CKDEPI/ risklimits; 
ods output  FitStatistics=Death_RMPM_Fit;
run;

proc sql;
select mean(WithCovariates) as AIC
from Death_RMPM_Fit
where Criterion='AIC';
quit;

proc sql;
select mean(WithCovariates) as BIC
from Death_RMPM_Fit
where Criterion=';SBC ;';
quit;
%mend;
%Death(esrd.Ctrgroupfollowup_mi, CVTime, CV_event, , ) /*RMPM*/
%Death(esrd.Ctrgroupfollowup_mi, CVTime, CV_event,  CCTR, ) /*RMPM +Continuous CTR*/
%Death(esrd.Ctrgroupfollowup_mi, CVTime, CV_event,  CTRgroup_New, CTRgroup_New(ref="1")) /*RMPM +Categorical CTR*/
%Death(esrd.Ctrgroupfollowup_mi, deathtime, death, , ) /*RMPM*/
%Death(esrd.Ctrgroupfollowup_mi, deathtime, death, CCTR, ) /*RMPM +Continuous CTR*/
%Death(esrd.Ctrgroupfollowup_mi, deathtime, death, CTRgroup_New, CTRgroup_New(ref="1")) /*RMPM +Categorical CTR*/

/***Single model's c-statistics***/
/*ESRD*/
%macro DialysisC (name1, name2);
proc sort data=esrd.Ctrgroupfollowup_mi;
by _imp;
run;

proc phreg data=esrd.Ctrgroupfollowup_mi concordance=uno(se seed=1234 iter=100);
by _imp;
class  &name2. gender_dw(ref="女");
model dialysistime*event(0)=&name1. age_entry gender_dw eGFR_CKDEPI PCRtotal_log;
ods output  ParameterEstimates=param_origi  Concordance=cstatistic_origi;
run;

proc mianalyze data=cstatistic_origi;
	modeleffects estimate;
	stderr stderr;
run;
%mend;
%DialysisC( , )
%DialysisC(CCTR, )
%DialysisC(CTRgroup_New, CTRgroup_New(ref="1"))

/*CV/All-cause mortality*/
%macro DeathC(name1, name2, name3, name4);
proc sort data=esrd.Ctrgroupfollowup_mi;
by _imp;
run;

proc phreg data=esrd.Ctrgroupfollowup_mi concordance=uno(se seed=1234 iter=100);
by _imp;
class  &name4. gender_dw(ref="女") havedm(ref='0') havehtn(ref='0') Amemia(ref='0');
model &name1.*&name2.(0)=&name3. age_entry gender_dw havedm havehtn Amemia eGFR_CKDEPI;
ods output  ParameterEstimates=param_origi  Concordance=cstatistic_origi;
run;

proc mianalyze data=cstatistic_origi;
	modeleffects estimate;
	stderr stderr;
run;
%mend;
%DeathC(CVTime, CV_event, , )
%DeathC(CVTime, CV_event, CCTR, )
%DeathC(CVTime, CV_event, CTRgroup_New, CTRgroup_New(ref="1"))
%DeathC(deathtime, death, , )
%DeathC(deathtime, death, CCTR, )
%DeathC(deathtime, death, CTRgroup_New, CTRgroup_New(ref="1"))

/***Compare different model's c-statistics***/
%macro DiffC(name1, name2, name3, name4, name5, name6, name7, name8);
proc sort data=esrd.Ctrgroupfollowup_mi;
by _imp;
run;

proc phreg data=esrd.Ctrgroupfollowup_mi concordance=uno (diff se seed=1234 iter=100);
    by _imp;
    class CTRgroup_New (ref="1") gender_dw(ref="女") &name4.;
    model &name1.*&name2.(0)=CCTR CTRgroup_New age_entry gender_dw &name3./nofit;
    roc "&name5." &name6. age_entry gender_dw &name3.;
    roc "&name5.+baseline CTR" &name7. age_entry gender_dw &name3.;
    roc "&name5.+baseline CTR group" &name8. age_entry gender_dw &name3.;
    ods output  ConcordanceDiff= Cdiff_statistic;
run;
proc sort data=Cdiff_statistic;by Source Source2;run;
proc mianalyze data=Cdiff_statistic;
	by Source Source2;
	modeleffects Estimate;
	stderr StdErr;
run;
%mend;
%DiffC(dialysistime, event, eGFR_CKDEPI PCRtotal_log, , KFRE,  , CCTR, CTRgroup_New)
%DiffC(CVTime, CV_event, havedm havehtn Amemia eGFR_CKDEPI, havedm(ref='0') havehtn(ref='0') Amemia(ref='0'), RMPM,  , CCTR, CTRgroup_New)
%DiffC(deathtime, death, havedm havehtn Amemia eGFR_CKDEPI, havedm(ref='0') havehtn(ref='0') Amemia(ref='0'), RMPM,  , CCTR, CTRgroup_New)

