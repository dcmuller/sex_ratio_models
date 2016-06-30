sink("./analysis/output/l00_get_data.rlog", split=TRUE)
#! p01_get_data.r

sessionInfo()
################################################
## case counts and rates are stored in separate datasets by year
casefiles <- lapply(paste0("./data/",
                           list.files(path="./data/",
                                      pattern="group.summary.*number.cases")),
                    FUN=read.csv, stringsAsFactors=FALSE, strip.white=TRUE)
popfiles <- lapply(paste0("./data/",
                          list.files(path="./data/",
                                     pattern="group.summary.*number.population")),
                   FUN=read.csv, stringsAsFactors=FALSE, strip.white=TRUE)
ci5_cases <- do.call(rbind.data.frame, casefiles)
ci5_pop <- do.call(rbind.data.frame, popfiles)
rm(casefiles, popfiles)
ci5_cases$REGISTRY <- NULL
ci5_pop$REGISTRY <- NULL
ci5_cases$N_unk<- NULL
ci5_pop$P_unk <- NULL
names(ci5_cases) <- tolower(names(ci5_cases))
names(ci5_pop) <- tolower(names(ci5_pop))
ci5_pop$icd <- NULL
ci5_pop$label <- NULL
ci5_pop <- ci5_pop[!duplicated(ci5_pop),]

## get rid of whitespace
ci5_cases$label <- trimws(ci5_cases$label)
ci5_pop$registry <- trimws(ci5_pop$registry)
ci5_cases$registry <- trimws(ci5_cases$registry)
ci5_cases$icd <- trimws(ci5_cases$icd)


## reshape to long format
widevars <- grep("^n[0-9]", names(ci5_cases), value=TRUE)
ci5_cases_long <- reshape(ci5_cases, varying=list(widevars), direction="long", timevar="age_group")
ci5_cases_long$age_group_start <- 5*as.numeric(ci5_cases_long$age_group) - 5
ci5_cases_long$n_ca <- ci5_cases_long$n0_4
ci5_cases_long$id <- NULL
ci5_cases_long$n0_4 <- NULL
ci5_cases_long$total <- NULL
widevars <- gsub("n", "", widevars)
widevars <- gsub("_", "-", widevars)
ci5_cases_long$age_group <- factor(ci5_cases_long$age_group, labels=widevars)

widevars <- grep("^p[0-9]", names(ci5_pop), value=TRUE)
ci5_pop_long <- reshape(ci5_pop, varying=list(widevars), direction="long", timevar="age_group")
ci5_pop_long$age_group_start <- 5*as.numeric(ci5_pop_long$age_group) - 5
ci5_pop_long$n_pop <- ci5_pop_long$p0_4
ci5_pop_long$id <- NULL
ci5_pop_long$p0_4 <- NULL
ci5_pop_long$total <- NULL
ci5_pop_long$n_agr <- NULL

widevars <- gsub("p", "", widevars)
widevars <- gsub("_", "-", widevars)
ci5_pop_long$age_group <- factor(ci5_pop_long$age_group, labels=widevars)

##
ci5_long <- merge(ci5_cases_long, ci5_pop_long)

## Age restriction: counts above 80 get sketchy, and below 40 will be
## contaminated by cancers related to genetic syndromes
ci5_long <- ci5_long[ci5_long$age_group_start>=40 & ci5_long$age_group_start<80, ]
ci5_long$age_group <- droplevels(ci5_long$age_group)

## drop registries with potential overlap
t(t(table(ci5_long$registry)))
todrop <- c("Italy, Florence", "Singapore: Chinese", "Singapore: Malay",
            "The Netherlands, Eindhoven",
            "USA, California, Los Angeles: Black",
            "USA, California, Los Angeles: Chinese",
            "USA, California, Los Angeles: Filipino",
            "USA, California, Los Angeles: Hispanic White",
            "USA, California, Los Angeles: Japanese",
            "USA, California, Los Angeles: Korean",
            "USA, California, Los Angeles: Non-Hispanic White",
            "USA, California, San Francisco: Black",
            "USA, California, San Francisco: White",
            "USA, Connecticut: Black",
            "USA, Connecticut: White",
            "USA, Georgia, Atlanta: Black",
            "USA, Georgia, Atlanta: White",
            "USA, Hawaii",
            "USA, Iowa",
            "USA, Michigan, Detroit: Black",
            "USA, Michigan, Detroit: White",
            "USA, New Jersey",
            "USA, New Jersey: Black",
            "USA, New Jersey: White",
            "USA, New Mexico",
            "USA, New York State",
            "USA, Utah",
            "USA, Washington, Seattle",
            "USA, SEER (9 registries): Black",
            "USA, SEER (9 registries): White")

ci5_long <- ci5_long[!(ci5_long$registry %in% todrop), ]

## look at registries by year
registries <- ci5_long[, c("year", "registry")]
registries <- registries[!duplicated(registries), ]
ci5_long$registry_cde <- factor(ci5_long$registry)

## start in 1978 -- 26 registries at least and gives us 30 years
rm(registries)
ci5_long <- ci5_long[ci5_long$year>1977, ]
ci5_long$year_cde <- factor(ci5_long$year)

## calculate raw rates
ci5_long$rate_raw <- with(ci5_long, 100000*n_ca/n_pop)

## recode sex
ci5_long$male <- as.numeric(ci5_long$sex==1)
with(ci5_long, table(sex, male))
ci5_long$sex <- factor(ci5_long$male, labels=c("Women", "Men"))
with(ci5_long, table(sex, male))

## get rid of irrelevant cancers (breast, prostate, etc), and total
todrop <- c("Breast", "Cervix uteri", "Corpus uteri",
            "Ovary and other uterine adnexa", "Prostate", "Testis",
            "All sites but non-melanoma skin")
ci5_long <- ci5_long[!(ci5_long$label %in% todrop), ]

## get rid of Bladder cancer, which has had weird and changing classification
## over periods and regions
todrop <- c("Bladder")
ci5_long <- ci5_long[!(ci5_long$label %in% todrop), ]

## Factor for cancer type
ci5_long$cancer_type <- factor(ci5_long$label)

## set up key-value data frames
registry_key <- data.frame(registry_id=as.numeric(ci5_long$registry_cde),
                           registry_cde=ci5_long$registry_cde)
registry_key <- registry_key[!duplicated(registry_key), ]
registry_key <- registry_key[order(registry_key$registry_id), ]
rownames(registry_key) <- NULL

age_key <- data.frame(age_id=as.numeric(ci5_long$age_group),
                      age_group=ci5_long$age_group)
age_key <- age_key[!duplicated(age_key), ]
age_key <- age_key[order(age_key$age_id), ]
rownames(age_key) <- NULL

year_key <- data.frame(year_id=as.numeric(ci5_long$year_cde),
                       year=ci5_long$year_cde)
year_key <- year_key[!duplicated(year_key), ]
year_key <- year_key[order(year_key$year_id), ]
rownames(year_key) <- NULL

cancer_key <- data.frame(cancer_id=as.numeric(ci5_long$cancer_type),
                         cancer_type=ci5_long$cancer_type)
cancer_key <- cancer_key[!duplicated(cancer_key), ]
cancer_key <- cancer_key[order(cancer_key$cancer_id), ]
rownames(cancer_key) <- NULL

## split data by cancer type
ci5_bycancer <- split(ci5_long, ci5_long$cancer_type)

################################
save(ci5_long, ci5_bycancer, age_key, year_key, registry_key, cancer_key, file="./data/d00_ci5.Rdta")
sink()