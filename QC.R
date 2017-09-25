# QC

qc <- getQC(mSetSw)
plotQC(qc, id=colnames(mSetSw))
densityPlot(mSetSw, sampGroups = subtype)
densityBeanPlot(mSet, sampGroups = subtype)
controlStripPlot(rgSet, controls="BISULFITE CONVERSION II")
qcReport(rgSet, pdf= "~/Koop_Klinghammer/Results/qcReport.pdf")

GRset <- mapToGenome(mSetSw)
predictedSex <- getSex(GRset, cutoff = -2)$predictedSex
head(predictedSex)
plotSex(getSex(GRset, cutoff = -2), id = s_names)
