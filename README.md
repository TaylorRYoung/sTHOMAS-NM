# sTHOMAS-NM
Code and trained normative models of thalamic nuclei and subcortical volumes generated with sTHOMAS (Saranathan et al. 2025). Models were trained on a cohort of 2,619 healthy controls ages 18-95 from 7 neuroimaging cohorts. Per Young et al. 2024, a multivariable fractional polynomial (MFP) normative model was fitted for each structure bilaterally with the R package mfp (Young TR et al. 2024).

Figure 1 below demonstrates a shift method for applying the models to a new site (Balakrishnan et al. 2025). New site data provided courtesy of Theodorus Van Erp and Annalisa Lella at University of California Irvine.

Full reference in preparation.

<img width="3456" height="2112" alt="image" src="https://github.com/user-attachments/assets/b4fbdda3-9b66-457a-b933-28728d17ec43" />
__Figure 1__ - Z-scores were generated for the control samples of a new site and demonstrate potential site effects notably in L13 and R13. Site effects were accounted for by the shift method or retraining the models incorporating controls from the new site. Both methods appear to address site effects equally and have the same impact on a population of study participants with schizophrenia from the same site. 


References:
1. Saranathan M, Cogliandro G, Hicks T, et al. Comprehensive Segmentation of Deep Grey Nuclei From Structural MRI Data. Human Brain Mapping. 2025;46(14):e70350. doi:10.1002/hbm.70350
2. Young TR, Kumar VJ, Saranathan M. Normative Modeling of Thalamic Nuclear Volumes and Characterization of Lateralized Volume Alterations in Alzheimer’s Disease Versus Schizophrenia. Biological Psychiatry: Cognitive Neuroscience and Neuroimaging. 2025;10(7):726-739. doi:10.1016/j.bpsc.2024.08.006
3. S. Balakrishnan, A. Banerjee, T. Young, M. Saranathan. Applications of Normative Modeling: A Shift Method for Site-Independent Normative
Modeling and Brain Age Prediction. International Society of Magnetic Resonance in Medicine Annual Meeting. Honolulu, HI. 2025. DOI: https://doi.org/10.58530/2025/4820


