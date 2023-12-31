## Background
[Tempus](https://github.com/tempuslabs) has developed a gene expression profiler (GEP) to predict the risk of distant recurrence in early stage endometrioid endometrial carcinoma (EEC) patients with a specific focus on high-intermediate risk (HIR) patients.  The GEP is an RNA-Sequencing based assay that was trained with TCGA data and independently validated on Tempus data and an external cohort of early EEC patients.  These validation studies demonstrated that the Tempus GEP molecular risk (MR) score correlates with known pathologic features and shows significant stratification of distant recurrence free survival (DRFS) in the MR-High versus MR-Low risk groups in all early stage EEC patients, the high-intermediate risk subgroup, and the TP53wt/POLEwt subgroup. 

## Primary Endpoint
To clinically validate this novel GEP MR risk predictor, Tempus will use samples and outcome data from HIR patients in the GOG-210 study to determine the difference in 4-year DRFS between GEP status (MR-High/MR-Low). Below is the methodology used to determine the number of samples required from the GOG-210 study and the R code to make these calculations these are provided herein.

## Sample Size Determination
Sample size and power for the primary endpoint were calculated using 10,000 simulations under the following assumptions:
* [PORTEC-2](https://doi.org/10.1016/S0140-6736(09)62163-2) and [GOG-99](https://doi.org/10.1016/j.ygyno.2003.11.048) reported distant recurrence rates between 5.7% to 10%. A 4-year distant recurrence rate of 5.7% was assumed in order to conservatively estimate power.
* Based on the Stanford-Howitt training dataset, a sensitivity of 54.5% and specificity of 78.5% was considered and an additional scenario with sensitivity of 49.5% and specificity of 73.5% was considered
* Yearly attrition rate of 15%

Under these assumptions, power and the sample size required to test the primary endpoint at the one-sided 5% significance can be determined.

## Time to recurrence data generation method
1. Based on the assumed 4-year recurrence rate (4YR), the hazard function λ is obtained by -log(1-4YR)/48 and the time to recurrence (time<sub>R</sub>) is generated by the exponential distribution Exp( λ ).
2. The attrition time (time<sub>D</sub>) for each patient is generated in a similar way based on the assumed yearly attrition rate.
3. 4-year recurrence event (event) indicator is obtained by:
    * event = 1 if time<sub>R</sub> ≤ time<sub>D</sub> and time<sub>R</sub> < 48,
    * event = 0  otherwise.
4. Follow-up time (FU) for each patient is obtained by:
    * FU = time<sub>R</sub>  if event=1,
    * FU = min(time<sub>D</sub> , 48)  if event=0.
5. MR-High patients for the 4-year recurrent patients (TP) will be generated by the binomial distribution B(N<sub>R</sub>, Se), where N<sub>R</sub> is the number of patients with event=1, and Se is the assumed sensitivity for the simulation.
6. MR-Low patients for the 4-year recurrence free patients (TN) will be generated by the binomial distribution B(N<sub>RF</sub>, Sp), where N<sub>RF</sub> is the number of patients with event=0, and Sp is the assumed specificity for the simulation.
7. GEP status can be generated by: 
    * GEP = MR-High if (event=1 and TP=1) or (event=0 and TN=0),
    * GEP = MR-Low otherwise.
8. This will create a simulated data frame with FU, event, and GEP.
