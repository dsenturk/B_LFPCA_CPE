CONTENTS OF THIS FOLDER ——————————————

B_LFPCA_CPE_tutorial.Rmd: A step-by-step implementation of obtaining functional depth based central posterior envelopes (CPEs) for posterior functional components in the Bayesian longitudinal functional principal component analysis (B-LFPCA) model detailed in "Central posterior envelopes for Bayesian longitudinal functional principal component analysis" by Boland et al. (2024). This procedure assumes that the data is comprised of densely observed longitudinal functional observations. 

B_LFPCA_CPE_simulateData.R: Function for simulating longitudinal functional data using the product FPCA model, as introduced in Section 5 of the main paper and Supplemental Materials Appdendix D of "Central posterior envelopes for Bayesian longitudinal functional principal component analysis".

B_LFPCA_CPE_MCMC.R: Function including MCMC steps in the Bayesian longitudinal functional principal component analysis (B-LFPCA) model to obtain posterior estimates, as described in Section 2 of "Central posterior envelopes for Bayesian longitudinal functional principal component analysis".

B_LFPCA_CPE_postSumms.R: Function to calculate traditional and proposed depth-based posterior summaries described in Section 3 and 4 of "Central posterior envelopes for Bayesian longitudinal functional principal component analysis".

B_LFPCA_CPE_muCPEcontour.R: Function to obtain the visualization of CPE contours for the mean function described in Section 3 of "Central posterior envelopes for Bayesian longitudinal functional principal component analysis".

B_LFPCA_CPE_CPEcontour.R: Function to obtain the visualization of CPE contours for the marginal longitudinal/functional eigenfunctions described in Section 3 of "Central posterior envelopes for Bayesian longitudinal functional principal component analysis".

INTRODUCTION ——————————————

The contents of this folder allow for implementation of obtaining functional depth based central posterior envelopes (CPEs) for posterior functional components in the Bayesian longitudinal functional principal component analysis (B-LFPCA) model detailed in "Central posterior envelopes for Bayesian longitudinal functional principal component analysis" by Boland et al. (2024). 
Users can simulate a sample dataset comprised of densely observed longitudinal functional observations (B_LFPCA_CPE_simulateData.R) and apply the Bayesian longitudinal functional principal component analysis (B-LFPCA) model to obtain posterior estimates of model parameters (B_LFPCA_CPE_MCMC.R). We include tools to calculate traditional and proposed depth-based posterior summaries (B_LFPCA_CPE_postSumms.R), and obtain the visualization of the proposed CPE contours for the mean function and marginal longitudinal/functional eigenfunctions (B_LFPCA_CPE_muCPEcontour.R, B_LFPCA_CPE_CPEcontour.R). Detailed instructions on how to perform the aforementioned procedures, including making CPE contour plots are included in B_LFPCA_CPE_tutorial.Rmd.

REQUIREMENTS ——————————————

The included R programs require R 4.4.3 (R Core Team, 2024) and the packages listed in B_LFPCA_CPE_tutorial.Rmd.

INSTALLATION ——————————————

Load the R program files into the global environment and install required packages using commands in B_LFPCA_CPE_tutorial.Rmd.
