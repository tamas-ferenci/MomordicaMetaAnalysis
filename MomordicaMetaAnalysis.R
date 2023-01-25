library(data.table)

RawData <- as.data.table(readxl::read_excel("STAT momo 2022-10-15.xlsx", sheet = 1,
                                            .name_repair = "universal"))

table(is.na(RawData$Before), is.na(RawData$After), is.na(RawData$Change))

RawData[!is.na(Change)&!is.na(Before)&!is.na(After),
        .((BeforeSD^2+AfterSD^2-ChangeSD^2)/(2*BeforeSD*AfterSD)), .(Study, Arm, Outcome)]

EstCorrs <- RawData[!is.na(Change)&!is.na(Before)&!is.na(After),
                    .((BeforeSD^2+AfterSD^2-ChangeSD^2)/(2*BeforeSD*AfterSD)),
                    .(Study, Arm, Outcome)]

dcast(EstCorrs, Outcome ~ Study + Arm, value.var = "V1")

ImputedCorr <- EstCorrs[, .(ImputedCorr = if(diff(range(V1, na.rm = TRUE))<0.4 & all(V1<1 & V1>-1)) mean(V1) else NA_real_),
                        .(Outcome, Arm)]
ImputedCorr

write.csv2(merge(EstCorrs, ImputedCorr), "EstCorrs.csv", row.names = FALSE)

RawData <- merge(RawData, ImputedCorr, by = c("Outcome", "Arm"), sort = FALSE, all = TRUE)

RawData$ImputedChange <- RawData$After-RawData$Before
RawData$ImputedChangeSD <- sqrt(RawData$BeforeSD^2 + RawData$AfterSD^2 -
                                  (2*RawData$ImputedCorr*RawData$BeforeSD*RawData$AfterSD))

RawData$ChangeFinal <- ifelse(!is.na(RawData$Change), RawData$Change, RawData$ImputedChange)
RawData$ChangeSDFinal <- ifelse(!is.na(RawData$Change), RawData$ChangeSD, RawData$ImputedChangeSD)

RawData$ImputedAfter <- RawData$Before + RawData$Change
RawData$ImputedAfterSD1 <- RawData$ImputedCorr*RawData$BeforeSD +
  sqrt(RawData$ImputedCorr^2*RawData$BeforeSD^2-RawData$BeforeSD^2+RawData$ChangeSD^2)
RawData$ImputedAfterSD2 <- RawData$ImputedCorr*RawData$BeforeSD -
  sqrt(RawData$ImputedCorr^2*RawData$BeforeSD^2-RawData$BeforeSD^2+RawData$ChangeSD^2)

RawData[is.na(After)]
RawData[Outcome=="ALT"]
RawData[Outcome=="AST"]

RawData$ImputedAfterSD <- NA_real_
RawData[Study=="Trakoon- osot"&(Outcome=="ALT"|(Outcome=="AST"&Arm=="active"))]$ImputedAfterSD <-
  RawData[Study=="Trakoon- osot"&(Outcome=="ALT"|(Outcome=="AST"&Arm=="active"))]$ImputedAfterSD1

RawData$AfterFinal <- ifelse(!is.na(RawData$After), RawData$After, RawData$ImputedAfter)
RawData$AfterSDFinal <- ifelse(!is.na(RawData$After), RawData$AfterSD, RawData$ImputedAfterSD)

write.csv2(RawData, "STAT_processed.csv", row.names = FALSE)

###

for(outc in unique(RawData[!is.na(AfterFinal)&!is.na(AfterSDFinal)]$Outcome)) {
  dat <- dcast(RawData[Outcome==outc&!is.na(AfterFinal)&!is.na(AfterSDFinal),
                       .(Study, Arm, N.pre, N.post, AfterFinal, AfterSDFinal)], Study ~ Arm,
               value.var = c("N.pre", "N.post", "AfterFinal", "AfterSDFinal"))
  dat <- dat[!is.na(N.pre_active)&!is.na(N.pre_placebo)&!is.na(N.post_active)&!is.na(N.post_placebo)]
  res <- meta::metacont(n.e = dat$N.post_active, mean.e = dat$AfterFinal_active,
                        sd.e = dat$AfterSDFinal_active, n.c = dat$N.post_placebo,
                        mean.c = dat$AfterFinal_placebo, sd.c = dat$AfterSDFinal_placebo,
                        studlab = dat$Study, random = TRUE, sm = "MD")
  
  cairo_pdf(paste0("./results/Forest_MD_After_", outc, ".pdf"), width = 12, height = 9/2)
  if(outc=="HDL")
    meta::forest(res, label.left = "Favours control", label.right = "Favours Momordica",
                 label.e = "Momordica", digits.sd = 1, digits.mean = 1) else
                   meta::forest(res, label.right = "Favours control", label.left = "Favours Momordica",
                                label.e = "Momordica", digits.sd = 1, digits.mean = 1)
  dev.off()
  
  res <- meta::update.meta(res, sm = "SMD")
  cairo_pdf(paste0("./results/Forest_SMD_After_", outc, ".pdf"), width = 12, height = 9/2)
  if(outc=="HDL")
    meta::forest(res, label.left = "Favours control", label.right = "Favours Momordica",
                 label.e = "Momordica", digits.sd = 1, digits.mean = 1) else
                   meta::forest(res, label.right = "Favours control", label.left = "Favours Momordica",
                                label.e = "Momordica", digits.sd = 1, digits.mean = 1)
  dev.off()
}

for(outc in unique(RawData[!is.na(ChangeFinal)&!is.na(ChangeSDFinal)]$Outcome)) {
  dat <- dcast(RawData[Outcome==outc&!is.na(ChangeFinal)&!is.na(ChangeSDFinal),
                       .(Study, Arm, N.pre, N.post, ChangeFinal, ChangeSDFinal)],
               Study ~ Arm, value.var = c("N.pre", "N.post", "ChangeFinal", "ChangeSDFinal"))
  dat <- dat[!is.na(N.pre_active)&!is.na(N.pre_placebo)&!is.na(N.post_active)&!is.na(N.post_placebo)]
  res <- meta::metacont(n.e = dat$N.post_active, mean.e = dat$ChangeFinal_active,
                        sd.e = dat$ChangeSDFinal_active, n.c = dat$N.post_placebo,
                        mean.c = dat$ChangeFinal_placebo, sd.c = dat$ChangeSDFinal_placebo,
                        studlab = dat$Study, random = TRUE, sm = "MD")
  
  cairo_pdf(paste0("./results/Forest_MD_Change_", outc, ".pdf"), width = 12, height = 9/2)
  if(outc=="HDL")
    meta::forest(res, label.left = "Favours control", label.right = "Favours Momordica",
                 label.e = "Momordica", digits.sd = 1, digits.mean = 1) else
                   meta::forest(res, label.right = "Favours control", label.left = "Favours Momordica",
                                label.e = "Momordica", digits.sd = 1, digits.mean = 1)
  dev.off()
  
  res <- meta::update.meta(res, sm = "SMD")
  cairo_pdf(paste0("./results/Forest_SMD_Change_", outc, ".pdf"), width = 12, height = 9/2)
  if(outc=="HDL")
    meta::forest(res, label.left = "Favours control", label.right = "Favours Momordica",
                 label.e = "Momordica", digits.sd = 1, digits.mean = 1) else
                   meta::forest(res, label.right = "Favours control", label.left = "Favours Momordica",
                                label.e = "Momordica", digits.sd = 1, digits.mean = 1)
  dev.off()
}