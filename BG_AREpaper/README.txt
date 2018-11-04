
--------------------------------------------------------

Title: Selection of composite binary endpoints in clinical trials
Authors: Marta Bofill Roig, Guadalupe Gómez Melis
E-mail address: marta.bofill.roig@upc.edu, lupe.gomez@upc.edu

Software version: R version 3.3.3 (2017-03-06)
package ‘plyr’ version 1.8.4


--------------------------------------------------------
File list
--------------------------------------------------------

BofillGomez_BiometricalJournalScript.R		R code
DATABASE_Scenarios.csv				Data for the efficiency guidelines


--------------------------------------------------------
R code: BofillGomez_BiometricalJournalScript.R
--------------------------------------------------------

This is the R code used to carry out the figures, tables and results in this paper. 
This code is structured as follows:

(I) FUNCTIONS
R functions defined:
-  Bahadur.composite
Description: this function returns the probability of a composite binary endpoint given the probabilities of its components and the correlation between them.

- OR.composite
Description: this function returns the odds ratio for a composite binary endpoint given  the probabilities under control group and the odds ratio for the components of the composite 	endpoint and the correlation between them. 

- OR
Description: this function returns the odds ratio given the probability of observing the event in the control and treatment group. 

- Correlation bounds (correlation.min.function and correlation.max.function)
Description: these functions return the lower and upper bound of the correlation given the probabilities of the single events. 

- ARE.betaOR.binary.endpoints (beta parametrization) 
Description: this calculates the ARE method given the probabilities under control group and the odds ratio for the components of the composite endpoint, and the correlation between them. 

(II) ASYMPTOTIC RELATIVE EFFICIENCY
Code to reproduce Figure 1 in Section 4.

(III) CASE STUDY
Code used to create Figure 2 in Section 5.

(IV) STATISTICAL EFFICIENCY GUIDELINES
Code used to perform the efficiency guidelines described in Section 6.


Since the time required to create the considered scenarios for efficiency guidelines was 16.58h, the resulting database of these scenarios is attached (DATABASE_Scenarios.csv).


--------------------------------------------------------


