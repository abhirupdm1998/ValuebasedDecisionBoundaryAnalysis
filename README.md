# Introduction
Welcome to our GitHub repository for DBA! 
This GitHub repository has been created for easy communication of codes with likeminded researchers who wish to simulate standard errors of any target TPP hazard ratios using a anchor study checklist early on the health technology lifecycle. Further, this repository was also created as an accompanying deliverable to the companion paper by Dutta Majumdar et al. 
This repository comprises the reproducible R code and parameters of the fitted models for OS and PFS for all the 37 HR scenarios between (0.1 to 1.00 by-step = 0.025) mentioned in the accompanying paper. Below we have mentioned the each of the files and pre-requisites and any special notes or comments that we may have. Again we note, the target product profile (TPP) for our illustrative example is for wonderumab versus sunitinib for determining value-based thresholds. 

**Copyright**
The copyright and ownership of the code completely vests with the authors of this code and the actual paper. Please cite this page whenever reusing this code properly along with the authors: Abhirup Dutta Majumdar, Sekhar Dutta, Anns Thomas, Ronan Mahon and Subrata Bhattacharrya and any reproduction of this code is strictly prohibited. 

# Choosing the Anchor Study
For any time-to-event be it OS or PFS, the anchor study from the baseline risk is adjusted for the TPP study (novel HT) is followed as per the steps mentioned in the accompanying paper by Dutta Majumdar et al.
After the scoring, in our case, for both OS and PFS, [Checkmate-214 Motzer et al](https://pubmed.ncbi.nlm.nih.gov/35383908/) was selected as best suitable for modeling sunitinib's baseline risk. This step is conducted using summary data/information reported in respective trial/studies (should be the latest one across all candidate studies, see the companion paper) 
Most importantly, as the next step, the individual data for the sunitinib arm is extracted using [WebplotDigitizer app](https://automeris.io/) and reconstructed using the algorithm by [Guyot et al 2012](https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/1471-2288-12-9)


# How to use the code scripts
The user should first download the scripts and data in their own computer and save. Then the users can choose to copy the code in their local script and run the code in the following sequence mentioned below after setting their working directory correctly. The current analysis has been done using the [R software (version 4.4.2)](https://cran.r-project.org/bin/windows/base/old/4.4.2/) which was accessed through the [R studio IDE by Posit](https://posit.co/download/rstudio-desktop/) on a Windows 64-bit computing software having 16 GB memory and core i7 processor. 

Before, we proceed the following R packages were used in the code: [tidyverse](https://www.tidyverse.org/packages/), [survival](https://cran.r-project.org/web/packages/survival/index.html), [survminer](https://cran.r-project.org/web/packages/survminer/index.html), [flexsurv](https://cran.r-project.org/web/packages/flexsurv/index.html), [simsurv](https://cran.r-project.org/web/packages/simsurv/index.html), [ggplo2](https://cran.r-project.org/web/packages/ggplot2/index.html) and [MASS](https://cran.r-project.org/web/packages/MASS/index.html). It should be noted that any copyrights of these packages (if applicable) completely vests with the original creators and developers and not us. 

# For simulation of TPP SE of log-HR
After identifying the anchor study, follow the TPP simulaion [code](https://github.com/abhirupdm1998/Decision-boundary-analysis-/blob/main/TPP%20log-HR%20SE%20simulation%20code.R). This exercise has to be repeated for each time-to-event outcome to be modelled. The KM data used for this analysis can be found [here](https://github.com/abhirupdm1998/Decision-boundary-analysis-/blob/main/IPD_checkmate_214_PFS_Sunitinib.csv)

