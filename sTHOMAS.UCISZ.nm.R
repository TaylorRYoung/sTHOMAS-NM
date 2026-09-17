#### UCI SZ
library(data.table)
library(mfp)
library(ggplot2)
library(ggridges)

## Source the nm functions
source("~/Documents/imaging_scripts/bpcnni.nm.functions.R")

## Vector of region names corresponding to a normative model .rds file
## Create a vector of regions as they are named in the dataset, there's no L3 or R3
roi_num <- c(1,2,4:14,26:34)
roi_lst <- c(paste0("L",roi_num), paste0("R",roi_num))

## Model covariates
cv <- c("age","male","female","eTIV","site")

## Directory to the normative models
mod_dir <- "~/Documents/imaging_data/hcp_ep/sthomas_nm/"

## Read in existing nm data
straw <- fread("~/Documents/imaging_data/hcp_ep/sTHOMAS_raw.csv", stringsAsFactors = FALSE)
stscr <- fread("~/Documents/imaging_data/hcp_ep/sTHOMAS_mfp_z.csv", stringsAsFactors = FALSE)

## Prepare UCI SZ data
ucisz <- fread("~/Documents/imaging_data/uci_sz/ucisz_hipsthomas_wide_edit.csv", stringsAsFactors = FALSE)
ucisz <- ucisz[, `:=` (phenotype=ifelse(Dx==0,"Control","Scz"),
                       group=ifelse(Dx==0,"Control","Case"),
                       male=ifelse(sex==1,1,0),
                       female=ifelse(sex==2,1,0),
                       site="HCP Aging")]
ucisz <- ucisz[, c("id","phenotype","group",cv,roi_lst), with = FALSE]

## Predict
ucizz <- mfp_get_z_bpcnni(roi_lst,ucisz,mod_dir)
ucizm <- melt(ucizz, id.vars = names(ucizz)[!names(ucizz) %in% roi_lst],
              value.name = "z", variable.name = "roi")

## Predict and shift
ucism <- mfp_shift_z_bpcnni(roi_lst,ucisz,mod_dir)
ucism <- melt(ucism, id.vars = names(ucizz)[!names(ucizz) %in% roi_lst],
              value.name = "z", variable.name = "roi")
ucism <- ucism[, phenotype := paste(phenotype,"shift")]

#### Fit new models
uci_dir <- "~/Documents/imaging_data/uci_sz/sthomas_nm/"
sm.fit.lst <- lapply(roi_lst, mfp_train_save_bpcnni, 
                     dat=rbind(ucisz[, site := "UCI SZ"],straw, fill=TRUE), out_dir=uci_dir)
ucinm <- mfp_get_z_bpcnni(roi_lst,ucisz,uci_dir)
ucinm <- melt(ucinm, id.vars = names(ucisz)[!names(ucisz) %in% roi_lst],
              value.name = "z", variable.name = "roi")
ucinm <- ucinm[, phenotype := paste(phenotype,"retrain")]

ucicmb <- rbind(ucizm,ucism,ucinm, fill=TRUE)
ucicmb <- ucicmb[, phenotype:=factor(phenotype, levels = c("Control","Control shift","Control retrain",
                                                           "Scz","Scz shift","Scz retrain"))]

uciwde <- dcast(ucicmb,id+age+male+female+eTIV+site+roi ~ phenotype, value.var = "z")

#### Plot distributions
ggplot(ucicmb, aes(x = z, y = roi, fill = phenotype)) + 
  geom_density_ridges(alpha = 0.8, show.legend = FALSE, color="grey25",
                      rel_min_height = 0.05, scale=1.5, quantile_lines = TRUE, quantiles = 2) +
  scale_x_continuous(name = "Z-score", limits = c(-4,4)) +
  ylab("") +
  scale_y_discrete(limits = rev(roi_lst)) +
  scale_fill_viridis_d(option = "C", direction = 1) + # colourblind-safe colours
  facet_grid( ~ phenotype) +
  labs(title = "Z-score Distributions") +
  geom_vline(xintercept = 0, alpha=0.2) +
  theme_minimal()

