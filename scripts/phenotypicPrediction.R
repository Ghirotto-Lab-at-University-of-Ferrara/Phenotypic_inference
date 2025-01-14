args<-commandArgs(trailingOnly=TRUE)

sample_id<-as.character(args[1])
sample_csv<-read.csv(as.character(args[2]),header=T)

## Creazione data.frame ##
df_prediction <- sample_csv[, c("sampleid", "PBlueEye", "PIntermediateEye", "PBrownEye", "PBlondHair", "PBrownHair", "PRedHair", "PBlackHair", "PLightHair", "PDarkHair", "PVeryPaleSkin", "PPaleSkin", "PIntermediateSkin", "PDarkSkin", "PDarktoBlackSkin")]
rownames(df_prediction) <- df_prediction$sampleid
df_prediction$sampleid <- NULL
df_prediction <- round(df_prediction, 3)

##FUNZIONE PER PREDIRE COLORE DEGLI OCCHI##
predict_eye_colour <- function(PBlueEye, PIntermediateEye, PBrownEye) {
  if (PBlueEye > PIntermediateEye & 
      PBlueEye > PBrownEye & 
      PBlueEye >= 0.7) {
    return("Blue")
  } else if (PIntermediateEye > PBlueEye & 
             PIntermediateEye > PBrownEye & 
             PIntermediateEye >= 0.7) {
    return("Intermediate")
  } else if (PBrownEye > PIntermediateEye & 
             PBrownEye > PBlueEye & 
             PBrownEye >= 0.7) {
    return("Brown")
  } else {
    return("Not predicted")
  }
}

##PREDIZIONE OCCHI## 
df_prediction$predicted_eye_colour <- apply(df_prediction[, c("PBlueEye", "PIntermediateEye", "PBrownEye")], 1, function(x) {
  predict_eye_colour(x["PBlueEye"], x["PIntermediateEye"], x["PBrownEye"])
})

##FUNZIONE PER PREDIRE COLORE DEI CAPELLI##
predict_hair_colour <- function(PBlondHair, PBrownHair, PRedHair, PBlackHair, PLightHair, PDarkHair) {
  if (PBlackHair > PBlondHair & 
      PBlackHair > PBrownHair & 
      PBlackHair > PRedHair & 
      PBlackHair >= 0.7 & 
      PDarkHair >= 0.5) {
    return("Black")
  } else if (PBlackHair > PBlondHair & 
             PBlackHair > PBrownHair & 
             PBlackHair > PRedHair & 
             PBlackHair >= 0.7 & 
             PDarkHair < 0.5) {
    return("Black")
  } else if (PBlackHair > PBlondHair & 
             PBlackHair > PBrownHair & 
             PBlackHair > PRedHair & 
             PBlackHair < 0.7 & 
             PDarkHair >= 0.5) {
    return("Black")
  } else if (PBlackHair > PBlondHair & 
             PBlackHair > PBrownHair & 
             PBlackHair > PRedHair & 
             PBlackHair < 0.7 & 
             PDarkHair < 0.5) {
    return("Dark brown/black")
  } else if (PBlondHair > PBlackHair & 
             PBlondHair > PBrownHair & 
             PBlondHair > PRedHair & 
             PBlondHair >= 0.7 & 
             PLightHair >= 0.95) {
    return("Blond")
  } else if (PBlondHair > PBlackHair & 
             PBlondHair > PBrownHair & 
             PBlondHair > PRedHair & 
             PBlondHair >= 0.7 & 
             PLightHair < 0.95) {
    return("Blond")
  } else if (PBlondHair > PBlackHair & 
             PBlondHair > PBrownHair & 
             PBlondHair > PRedHair & 
             PBlondHair < 0.7 & 
             PLightHair >= 0.9) {
    return("Blond")
  } else if (PBlondHair > PBlackHair & 
             PBlondHair > PBrownHair & 
             PBlondHair > PRedHair & 
             PBlondHair < 0.7 & 
             PLightHair < 0.9) {
    return("Dark blond")
  } else if (PBrownHair > PBlackHair & 
             PBrownHair > PBlondHair & 
             PBrownHair > PRedHair & 
             PBrownHair >= 0.7 & 
             PLightHair >= 0.8) {
    return("Brown")
  } else if (PBrownHair > PBlackHair & 
             PBrownHair > PBlondHair & 
             PBrownHair > PRedHair & 
             PBrownHair >= 0.7 & 
             PLightHair < 0.8) {
    return("Brown/dark brown")
  } else if (PBrownHair > PBlackHair & 
             PBrownHair > PBlondHair & 
             PBrownHair > PRedHair & 
             PBrownHair < 0.7 & 
             PLightHair >= 0.8) {
    return("Brown/dark brown")
  } else if (PBrownHair > PBlackHair & 
             PBrownHair > PBlondHair & 
             PBrownHair > PRedHair & 
             PBrownHair < 0.7 & 
             PLightHair < 0.8) {
    return("Dark brown/black")
  } else if (PRedHair > PBlackHair & 
             PRedHair > PBlondHair & 
             PRedHair > PBrownHair & 
             PRedHair >= 0.7 & 
             PLightHair >= 0.9) {
    return("Red")
  } else if (PRedHair > PBlackHair & 
             PRedHair > PBlondHair & 
             PRedHair > PBrownHair & 
             PRedHair >= 0.7 & 
             PLightHair < 0.9) {
    return("Red")
  } else if (PRedHair > PBlackHair & 
             PRedHair > PBlondHair & 
             PRedHair > PBrownHair & 
             PRedHair < 0.7 & 
             PLightHair >= 0.9) {
    return("Red")
  } else if (PRedHair > PBlackHair & 
             PRedHair > PBlondHair & 
             PRedHair > PBrownHair & 
             PRedHair < 0.7 & 
             PLightHair < 0.9) {
    return("Red")
  } else {
    return("Not predicted")
  }
}

##PREDIZIONE CAPELLI## 
df_prediction$predicted_hair_colour <- apply(df_prediction[, c("PBlondHair", "PBrownHair", "PRedHair", "PBlackHair", "PLightHair", "PDarkHair")], 1, function(x) {
  predict_hair_colour(x["PBlondHair"], x["PBrownHair"], x["PRedHair"], x["PBlackHair"], x["PLightHair"], x["PDarkHair"])
})

#df_prediction$predicted_eye_colour <- "Unable to predict"

##FUNZIONE PER PREDIRE COLORE DELLA PELLE##
predict_skin_colour <- function(PVeryPaleSkin, PPaleSkin, PIntermediateSkin, PDarkSkin, PDarktoBlackSkin) {
  if (PVeryPaleSkin > PPaleSkin & 
      PVeryPaleSkin > PIntermediateSkin & 
      PVeryPaleSkin > PDarkSkin & 
      PVeryPaleSkin > PDarktoBlackSkin & 
      PVeryPaleSkin >= 0.9) {
    return("Very pale")
  } else if (PVeryPaleSkin > PPaleSkin & 
             PVeryPaleSkin > PIntermediateSkin & 
             PVeryPaleSkin > PDarkSkin & 
             PVeryPaleSkin > PDarktoBlackSkin & 
             PVeryPaleSkin >= 0.7 & 
             (PPaleSkin >= 0.15 | PIntermediateSkin >= 0.15 | PDarkSkin >= 0.15 | PDarktoBlackSkin >= 0.15)) {
    return("Pale")
  } else if (PVeryPaleSkin > PPaleSkin & 
             PVeryPaleSkin > PIntermediateSkin & 
             PVeryPaleSkin > PDarkSkin & 
             PVeryPaleSkin > PDarktoBlackSkin & 
             PVeryPaleSkin >= 0.7 & 
             (PPaleSkin < 0.15 | PIntermediateSkin < 0.15 | PDarkSkin < 0.15 | PDarktoBlackSkin < 0.15)) {
    return("Very pale")
  } else if (PVeryPaleSkin > PPaleSkin & 
             PVeryPaleSkin > PIntermediateSkin & 
             PVeryPaleSkin > PDarkSkin & 
             PVeryPaleSkin > PDarktoBlackSkin & 
             PVeryPaleSkin >= 0.5) {
    return("Very pale")
  } else if (PPaleSkin > PVeryPaleSkin & 
             PPaleSkin > PIntermediateSkin & 
             PPaleSkin > PDarkSkin & 
             PPaleSkin > PDarktoBlackSkin & 
             PPaleSkin >= 0.9) {
    return("Pale")
  } else if (PPaleSkin > PVeryPaleSkin & 
             PPaleSkin > PIntermediateSkin & 
             PPaleSkin > PDarkSkin & 
             PPaleSkin > PDarktoBlackSkin & 
             PPaleSkin >= 0.7 & 
             PVeryPaleSkin >= 0.15) {
    return("Pale")
  } else if (PPaleSkin > PVeryPaleSkin & 
             PPaleSkin > PIntermediateSkin & 
             PPaleSkin > PDarkSkin & 
             PPaleSkin > PDarktoBlackSkin & 
             PPaleSkin >= 0.7 & 
             (PIntermediateSkin >= 0.15 | PDarkSkin >= 0.15 | PDarktoBlackSkin >= 0.15)) {
    return("Intermediate")
  } else if (PPaleSkin > PVeryPaleSkin & 
             PPaleSkin > PIntermediateSkin & 
             PPaleSkin > PDarkSkin & 
             PPaleSkin > PDarktoBlackSkin & 
             PPaleSkin >= 0.7 & 
             (PVeryPaleSkin < 0.15 | PIntermediateSkin < 0.15 | PDarkSkin < 0.15 | PDarktoBlackSkin < 0.15)) {
    return("Pale")
  } else if (PPaleSkin > PVeryPaleSkin & 
             PPaleSkin > PIntermediateSkin & 
             PPaleSkin > PDarkSkin & 
             PPaleSkin > PDarktoBlackSkin & 
             PPaleSkin >= 0.5 & 
             (PVeryPaleSkin > PIntermediateSkin & PVeryPaleSkin > PDarkSkin & PVeryPaleSkin > PDarktoBlackSkin)) {
    return("Pale")
  } else if (PPaleSkin > PVeryPaleSkin & 
             PPaleSkin > PIntermediateSkin & 
             PPaleSkin > PDarkSkin & 
             PPaleSkin > PDarktoBlackSkin & 
             PPaleSkin >= 0.5 & 
             (PIntermediateSkin > PVeryPaleSkin | PDarkSkin > PVeryPaleSkin | PDarktoBlackSkin > PVeryPaleSkin)) {
    return("Pale")
  } else if (PIntermediateSkin > PVeryPaleSkin & 
             PIntermediateSkin > PPaleSkin & 
             PIntermediateSkin > PDarkSkin & 
             PIntermediateSkin > PDarktoBlackSkin & 
             PIntermediateSkin >= 0.9 & 
             PDarktoBlackSkin > PVeryPaleSkin & 
             PDarktoBlackSkin > PPaleSkin & 
             PDarktoBlackSkin > PDarkSkin) {
    return("Dark")
  } else if (PIntermediateSkin > PVeryPaleSkin & 
             PIntermediateSkin > PPaleSkin & 
             PIntermediateSkin > PDarkSkin & 
             PIntermediateSkin > PDarktoBlackSkin & 
             PIntermediateSkin >= 0.9 & 
             (PDarktoBlackSkin < PVeryPaleSkin | PDarktoBlackSkin < PPaleSkin | PDarktoBlackSkin < PDarkSkin)) {
    return("Intermediate")
  } else if (PIntermediateSkin > PVeryPaleSkin & 
             PIntermediateSkin > PPaleSkin & 
             PIntermediateSkin > PDarkSkin & 
             PIntermediateSkin > PDarktoBlackSkin & 
             PIntermediateSkin >= 0.7 & 
             (PVeryPaleSkin >= 0.15 | PPaleSkin >= 0.15) & 
             (PDarkSkin < 0.15 | PDarktoBlackSkin < 0.15)) {
    return("Pale")
  } else if (PIntermediateSkin > PVeryPaleSkin & 
             PIntermediateSkin > PPaleSkin & 
             PIntermediateSkin > PDarkSkin & 
             PIntermediateSkin > PDarktoBlackSkin & 
             PIntermediateSkin >= 0.7 & 
             (PVeryPaleSkin < 0.15 | PPaleSkin < 0.15) & 
             (PDarkSkin >= 0.15 | PDarktoBlackSkin >= 0.15)) {
    return("Dark")
  } else if (PIntermediateSkin > PVeryPaleSkin & 
             PIntermediateSkin > PPaleSkin & 
             PIntermediateSkin > PDarkSkin & 
             PIntermediateSkin > PDarktoBlackSkin & 
             PIntermediateSkin >= 0.7 & 
             (PVeryPaleSkin < 0.15 | PPaleSkin < 0.15 | PDarkSkin < 0.15 | PDarktoBlackSkin < 0.15)) {
    return("Intermediate")
  } else if (PIntermediateSkin > PVeryPaleSkin & 
             PIntermediateSkin > PPaleSkin & 
             PIntermediateSkin > PDarkSkin & 
             PIntermediateSkin > PDarktoBlackSkin & 
             PIntermediateSkin >= 0.5 & 
             (PDarktoBlackSkin > PVeryPaleSkin & PDarktoBlackSkin > PPaleSkin | PDarkSkin > PVeryPaleSkin & PDarkSkin > PPaleSkin)) {
    return("Dark")
  } else if (PIntermediateSkin > PVeryPaleSkin & 
             PIntermediateSkin > PPaleSkin & 
             PIntermediateSkin > PDarkSkin & 
             PIntermediateSkin > PDarktoBlackSkin & 
             PIntermediateSkin >= 0.5 & 
             (PDarktoBlackSkin < PVeryPaleSkin | PDarktoBlackSkin < PPaleSkin | PDarkSkin < PVeryPaleSkin | PDarkSkin < PPaleSkin)) {
    return("Intermediate")
  } else if (PDarkSkin > PVeryPaleSkin & 
             PDarkSkin > PPaleSkin & 
             PDarkSkin > PIntermediateSkin & 
             PDarkSkin > PDarktoBlackSkin & 
             PDarkSkin >= 0.9) {
    return("Dark")
  } else if (PDarkSkin > PVeryPaleSkin & 
             PDarkSkin > PPaleSkin & 
             PDarkSkin > PIntermediateSkin & 
             PDarkSkin > PDarktoBlackSkin & 
             PDarkSkin >= 0.7 & 
             PDarktoBlackSkin > PVeryPaleSkin & 
             PDarktoBlackSkin > PPaleSkin & 
             PDarktoBlackSkin > PIntermediateSkin) {
    return("Dark to black")
  } else if (PDarkSkin > PVeryPaleSkin & 
             PDarkSkin > PPaleSkin & 
             PDarkSkin > PIntermediateSkin & 
             PDarkSkin > PDarktoBlackSkin & 
             PDarkSkin >= 0.5 & 
             PDarktoBlackSkin > PVeryPaleSkin & 
             PDarktoBlackSkin > PPaleSkin & 
             PDarktoBlackSkin > PIntermediateSkin) {
    return("Dark to black")
  } else if (PDarkSkin > PVeryPaleSkin & 
             PDarkSkin > PPaleSkin & 
             PDarkSkin > PIntermediateSkin & 
             PDarkSkin > PDarktoBlackSkin & 
             PDarkSkin >= 0.5 & 
             PVeryPaleSkin > PPaleSkin & 
             PVeryPaleSkin > PIntermediateSkin & 
             PVeryPaleSkin > PDarktoBlackSkin) {
    return("Intermediate")
  } else if (PDarkSkin > PVeryPaleSkin & 
             PDarkSkin > PPaleSkin & 
             PDarkSkin > PIntermediateSkin & 
             PDarkSkin > PDarktoBlackSkin & 
             PDarkSkin >= 0.5 & 
             PPaleSkin > PVeryPaleSkin & 
             PPaleSkin > PIntermediateSkin & 
             PPaleSkin > PDarktoBlackSkin) {
    return("Dark")
  } else if (PDarkSkin > PVeryPaleSkin & 
             PDarkSkin > PPaleSkin & 
             PDarkSkin > PIntermediateSkin & 
             PDarkSkin > PDarktoBlackSkin & 
             PDarkSkin >= 0.5 & 
             PIntermediateSkin > PVeryPaleSkin & 
             PIntermediateSkin > PPaleSkin & 
             PIntermediateSkin > PDarktoBlackSkin) {
    return("Dark")
  } else if (PDarktoBlackSkin > PVeryPaleSkin & 
             PDarktoBlackSkin > PPaleSkin & 
             PDarktoBlackSkin > PIntermediateSkin & 
             PDarktoBlackSkin > PDarkSkin & 
             PDarktoBlackSkin >= 0.9) {
    return("Dark to black")
  } else if (PDarktoBlackSkin > PVeryPaleSkin & 
             PDarktoBlackSkin > PPaleSkin & 
             PDarktoBlackSkin > PIntermediateSkin & 
             PDarktoBlackSkin > PDarkSkin & 
             PDarktoBlackSkin >= 0.7) {
    return("Dark to black")
  } else if (PDarktoBlackSkin > PVeryPaleSkin & 
             PDarktoBlackSkin > PPaleSkin & 
             PDarktoBlackSkin > PIntermediateSkin & 
             PDarktoBlackSkin > PDarkSkin & 
             PDarktoBlackSkin >= 0.5 & 
             PVeryPaleSkin > PPaleSkin & 
             PVeryPaleSkin > PIntermediateSkin & 
             PVeryPaleSkin > PDarkSkin) {
    return("Intermediate")
  } else if (PDarktoBlackSkin > PVeryPaleSkin & 
             PDarktoBlackSkin > PPaleSkin & 
             PDarktoBlackSkin > PIntermediateSkin & 
             PDarktoBlackSkin > PDarkSkin & 
             PDarktoBlackSkin >= 0.5 & 
             PPaleSkin > PVeryPaleSkin & 
             PPaleSkin > PIntermediateSkin & 
             PPaleSkin > PDarkSkin) {
    return("Intermediate")
  } else if (PDarktoBlackSkin > PVeryPaleSkin & 
             PDarktoBlackSkin > PPaleSkin & 
             PDarktoBlackSkin > PIntermediateSkin & 
             PDarktoBlackSkin > PDarkSkin & 
             PDarktoBlackSkin >= 0.5 & 
             PIntermediateSkin > PVeryPaleSkin & 
             PIntermediateSkin > PPaleSkin & 
             PIntermediateSkin > PDarkSkin) {
    return("Dark")
  } else if (PDarktoBlackSkin > PVeryPaleSkin & 
             PDarktoBlackSkin > PPaleSkin & 
             PDarktoBlackSkin > PIntermediateSkin & 
             PDarktoBlackSkin > PDarkSkin & 
             PDarktoBlackSkin >= 0.5 & 
             PDarkSkin > PVeryPaleSkin & 
             PDarkSkin > PPaleSkin & 
             PDarkSkin > PIntermediateSkin) {
    return("Dark to black")
  } else {
    return("Not predicted")
  }
}

##PREDIZIONE PELLE## 
df_prediction$predicted_skin_colour <- apply(df_prediction[, c("PVeryPaleSkin", "PPaleSkin", "PIntermediateSkin", "PDarkSkin", "PDarktoBlackSkin")], 1, function(x) {
  predict_skin_colour(x["PVeryPaleSkin"], x["PPaleSkin"], x["PIntermediateSkin"], x["PDarkSkin"], x["PDarktoBlackSkin"])
})

write.table(df_prediction, paste(sample_id,"_phenotypicPrediction.csv",sep=""), quote=F, sep=";")