context("cellMarkerToGmt")

library("GSEABase")

tmp <- tempdir()
infile1 = paste0(tmp, "/Human_cell_markers.txt")
infile2 = paste0(tmp, "/Mouse_cell_markers.txt")
infile3 = paste0(tmp, "/Single_cell_markers.txt")
infile4 = paste0(tmp, "/all_cell_markers.txt")

outfile1_1 = paste0(tmp, "/Human_cell_markers_1.gmt")
outfile1_2 = paste0(tmp, "/Human_cell_markers_2.gmt")
outfile1_3 = paste0(tmp, "/Human_cell_markers_3.gmt")
outfile1_4 = paste0(tmp, "/Human_cell_markers_4.gmt")
outfile2_1 = paste0(tmp, "/Mouse_cell_markers_1.gmt")
outfile2_2 = paste0(tmp, "/Mouse_cell_markers_2.gmt")
outfile2_3 = paste0(tmp, "/Mouse_cell_markers_3.gmt")
outfile2_4 = paste0(tmp, "/Mouse_cell_markers_4.gmt")
outfile3_1 = paste0(tmp, "/Single_cell_markers_1.gmt")
outfile3_2 = paste0(tmp, "/Single_cell_markers_2.gmt")
outfile3_3 = paste0(tmp, "/Single_cell_markers_3.gmt")
outfile3_4 = paste0(tmp, "/Single_cell_markers_4.gmt")
outfile4_1 = paste0(tmp, "/all_cell_markers_1.gmt")
outfile4_2 = paste0(tmp, "/all_cell_markers_2.gmt")
outfile4_3 = paste0(tmp, "/all_cell_markers_3.gmt")
outfile4_4 = paste0(tmp, "/all_cell_markers_4.gmt")

# Human_cell_markers.txt
sink(infile1)
cat("speciesType\ttissueType\tUberonOntologyID\tcancerType\tcellType\tcellName\tCellOntologyID\tcellMarker\tgeneSymbol\tgeneID\tproteinName\tproteinID\tmarkerResource\tPMID\tCompany\n")
cat("Human\tKidney\tUBERON_0002113\tNormal\tNormal cell\tProximal tubular cell\tNA\tIntestinal Alkaline Phosphatase\tALPI\t248\tPPBI\tP09923\tExperiment\t9263997\tNA\n")
cat("Human\tLiver\tUBERON_0002107\tNormal\tNormal cell\tIto cell (hepatic stellate cell)\tCL_0000632\tSynaptophysin\tSYP\t6855\tSYPH\tP08247\tExperiment\t10595912\tNA\n")
cat("Human\tEndometrium\tUBERON_0001295\tNormal\tNormal cell\tTrophoblast cell\tCL_0000351\tCEACAM1\tCEACAM1\t634\tCEAM1\tP13688\tExperiment\t10751340\tNA\n")
cat("Human\tGerm\tUBERON_0000923\tNormal\tNormal cell\tPrimordial germ cell\tCL_0000670\tVASA\tDDX4\t54514\tDDX4\tQ9NQI0\tExperiment\t10920202\tNA\n")
cat("Human\tCorneal epithelium\tUBERON_0001772\tNormal\tNormal cell\tEpithelial cell\tCL_0000066\tKLF6\tKLF6\t1316\tKLF6\tQ99612\tExperiment\t12407152\tNA\n")
cat("Human\tPlacenta\tUBERON_0001987\tNormal\tNormal cell\tCytotrophoblast\tCL_0000351\tFGF10\tFGF10\t2255\tFGF10\tO15520\tExperiment\t15950061\tNA\n")
cat("Human\tPeriosteum\tUBERON_0002515\tNormal\tNormal cell\tPeriosteum-derived progenitor cell\tNA\tCD166, CD45, CD9, CD90\tALCAM, PTPRC, CD9, THY1\t214, 5788, 928, 7070\tCD166, PTPRC, CD9, THY1\tQ13740, P08575, P21926, P04216\tExperiment\t15977065\tNA\n")
cat("Human\tAmniotic membrane\tUBERON_0009742\tNormal\tNormal cell\tAmnion epithelial cell\tCL_0002536\tNANOG, OCT3/4\tNANOG, POU5F1\t79923, 5460\tNANOG, PO5F1\tQ9H9S0, Q01860\tExperiment\t16081662\tNA\n")
cat("Human\tPrimitive streak\tUBERON_0004341\tNormal\tNormal cell\tPrimitive streak cell\tNA\tLHX1, MIXL1\tLHX1, MIXL1\t3975, 83881\tLHX1, MIXL1\tP48742, Q9H2W2\tExperiment\t16258519\tNA\n")
sink()

# Mouse_cell_markers.txt
sink(infile2)
cat("speciesType\ttissueType\tUberonOntologyID\tcancerType\tcellType\tcellName\tCellOntologyID\tcellMarker\tgeneSymbol\tgeneID\tproteinName\tproteinID\tmarkerResource\tPMID\tCompany\n")
cat("Mouse\tBone marrow\tUBERON_0002371\tNormal\tNormal cell\tBlastema cell\tCL_0000354\tCD31, Sca-1, Vim\tPecam1, Ly6a, Vim\t18613, 110454, 22352\tPECA1, LY6A, VIME\tQ08481, P05533, P20152\tExperiment\t29105393\tNA\n")
cat("Mouse\tTaste bud\tUBERON_0001727\tNormal\tNormal cell\tType II taste bud cell\tCL_0002285\tGustducin, PLC beta2\tPlcb2, Gnat3\t18796, 242851\tGNAT3, PLCB2\tA3KGF7, Q3V3I2\tExperiment\t29525940NA\n")
cat("Mouse\tTaste bud\tUBERON_0001727\tNormal\tNormal cell\tType III taste bud cell\tCL_0002285\tAADC, CA4, GAD67, NCAM, SNAP25\tDdc, Car4, Gad1, Ncam1, Snap25\t13195, 12351, 14451, 17967, 20614\tDDC, CAH4, DCE1, NCAM1, SNP25\tO88533, P22748, P48318, P13595, P60879\tExperiment\t29525940\tNA\n")
cat("Mouse\tAdipose tissue\tUBERON_0001013\tNormal\tNormal cell\tAdipose-derived stem cell\tCL_0000034\tSca-1\tLy6a\t110454\tLY6A\tP05533\tExperiment\t24726953\tNA\n")
cat("Mouse\tUndefined\tNA\tNormal\tNormal cell\tRegulatory T (Treg) cell\tCL_0000815\tCD49b, LAG3\tItga2, Lag3\t16398, 16768\tITA2, LAG3\tQ62469, P18627\tExperiment\t28929191\tNA\n")
cat("Mouse\tAdipose tissue\tUBERON_0001013\tNormal\tNormal cell\tWhite fat cell\tCL_0000448\tAsc-1\tSlc7a10\t53896\tAAA1\tP63115\tExperiment\t25080478\tNA\n")
cat("Mouse\tMeniscus\tUBERON_0000387\tNormal\tNormal cell\tMeniscus-derived stem cell\tNA\tCD105, CD34, CD44, CD73, CD90, Sca-1\tEng, Cd34, Cd44, Nt5e, Thy1, Ly6a\t13805, 12490, 12505, 23959, 21838, 110454\tEGLN, CD34, CD44, 5NTD, THY1, LY6A\tQ63961, Q64314, P15379, Q61503, P01831, P05533\tExperiment\t28005443\tNA\n")
cat("Mouse\tMeniscus\tUBERON_0000387\tNormal\tNormal cell\tMeniscus-derived progenitor cell\tNA\tCD105, CD34, CD44, CD73, CD90, Sca-1\tEng, Cd34, Cd44, Nt5e, Thy1, Ly6a\t13805, 12490, 12505, 23959, 21838, 110454\tEGLN, CD34, CD44, 5NTD, THY1, LY6A\tQ63961, Q64314, P15379, Q61503, P01831, P05533\tExperiment\t28005443\tNA\n")
cat("Mouse\tPeyer patch\tUBERON_0001211\tNormal\tNormal cell\tB cell\tCL_0000236\tCD45R\tPtprc\t19264\tPTPRC\tP06800\tExperiment\t24182520\tNA\n")
sink()

# Single_cell_markers.txt
sink(infile3)
cat("speciesType\ttissueType\tUberonOntologyID\tcancerType\tcellType\tcellName\tCellOntologyID\tcellMarker\tgeneSymbol\tgeneID\tproteinName\tproteinID\tmarkerResource\tPMID\tCompany\n")
cat("Human\tLung\tUBERON_0002048\tNormal\tNormal cell\tType II pneumocyte\tCL_0002063\tABCA3, SLC34A2\tABCA3, SLC34A2\t21, 10568\tABCA3, NPT2B\tQ99758, O95436\tSingle-cell sequencing\t27942595\tNA\n")
cat("Human\tLung\tUBERON_0002048\tNormal\tNormal cell\tAirway secretory cell\tCL_1000272\tAQP3, MUC5AC, MUC5B, PIGR, SCGB1A1\tAQP3, MUC5AC, MUC5B, PIGR, SCGB1A1\t360, 4586, 727897, 5284, 7356\tAQP3, MUC5A, MUC5B, PIGR, UTER\tQ92482, P98088, Q9HC84, P01833, P11684\tSingle-cell sequencing\t27942595\tNA\n")
cat("Human\tPancreas\tUBERON_0001264\tNormal\tNormal cell\tAlpha cell\tCL_0000171\tARX, CLIM1, CRYBA2, FEV, GBA, HMGB3, ILTMP, IRX2, KRP1, LOXL4, MAFB, PGE2-R, PGR, PLCE1, RFX6, RGP4, SMARCA1\tUBA2, LDB2, CRYBA2, FEV, GBA, HMGB3, TM4SF4, IRX2, KLHL41, LOXL4, MAFB, PTGER3, PGR, PLCE1, RFX6, RGPD4, SMARCA1\t10054, 9079, 1412, 54738, 2629, 3149, 7104, 153572, 10324, 84171, 9935, 5733, 5241, 51196, 222546, 285190, 6594\tARX, PDLI1, CRBA2, FEV, GLCM, HMGB3, T4S4, IRX2, KLH41, LOXL4, MAFB, PE2R3, PRGR, PLCE1, RFX6, RGS4, SMCA1\tQ96QS3, O00151, P53672, Q99581, P04062, O15347, P48230, Q9BZI1, O60662, Q96JB6, Q9Y5Q3, P43115, P06401, Q9P212, Q8HWS3, P49798, P28370\tSingle-cell sequencing\t27693023\tNA\n")
cat("Human\tPancreas\tUBERON_0001264\tNormal\tNormal cell\tBeta cellCL_0000169\tBMP-5, CDKN1C, CRTR1, DLK1, NPTX2, PACAP, PDX1, PFKFB2, PIR, SIX2, SIX3, SMAD9, SYT13, TGFBR3\tBMP5, CDKN1C, TFCP2L1, DLK1, NPTX2, MZB1, PDX1, PFKFB2, PIR, SIX2, SIX3, SMAD9, SYT13, TGFBR3\t653, 1028, 29842, 8788, 4885, 51237, 3651, 5208, 8544, 10736, 6496, 4093, 57586, 7049\tBMP5, CDN1C, TF2L1, DLK1, NPTX2, PACA, ODPX, F262, PIR, SIX2, SIX3, SMAD9, SYT13, TGBR3\tP22003, P49918, Q9NZI6, P80370, P47972, P18509, O00330, O60825, O00625, Q9NPC8, O95343, O15198, Q7L8C5, Q03167\tSingle-cell sequencing\t27693023\tNA\n")
cat("Human\tPancreas\tUBERON_0001264\tNormal\tNormal cell\tDelta cell\tCL_0000173\tCHE1, ESE3B, ETV1, GABRG2, HER4, ISL1, LCORL, LEDGF, LEPR, NEC1, PDLIM4, POU3F1, PRBP, PRG4, RGS2, SFRP3, SHARP1\tAATF, EHF, ETV1, GABRG2, ERBB4, ISL1, LCORL, PSIP1, LEPROT, PCSK1, PDLIM4, POU3F1, RBP4, PRG4, RGS2, FRZB, BHLHE41\t26574, 26298, 2115, 2566, 2066, 3670, 254251, 11168, 54741, 5122, 8572, 5453, 5950, 10216, 5997, 2487, 79365\tAATF, EHF, ETV1, GBRG2, ERBB4, ISL1, LCORL, PSIP1, LEPR, NEC1, PDLI4, PO3F1, RET4, PLPR2, RGS2, SFRP3, BHE41\tQ9NY61, Q9NZC4, P50549, P18507, Q15303, P61371, Q8N3X6, O75475, P48357, P29120, P50479, Q03052, P02753, Q96GM1, P41220, Q92765, Q9C0J9\tSingle-cell sequencing\t27693023NA\n")
cat("Human\tPancreas\tUBERON_0001264\tNormal\tNormal cell\tPancreatic polypeptide cell\tCL_0002275\tAQP3, ARHGAP3, ARX, BHLHB26, BHLHB27, CARTPT, EGR3, ENTPD2, ETV1, MEIS1, MEIS2, PAX6, PTGFR, RBTN3, SERTM1, SLITRK6, THSD7A, ZNF506\tAQP3, CHN2, UBA2, ID2, ID4, CARTPT, EGR3, ENTPD2, ETV1, MEIS1, MEIS2, PAX6, PTGFR, LMO3, SERTM1, SLITRK6, THSD7A, ZNF506\t360, 1124, 10054, 3398, 3400, 9607, 1960, 954, 2115, 4211, 4212, 5080, 5737, 55885, 400120, 84189, 221981, 440515\tAQP3, CHIO, ARX, ID2, ID4, CART, EGR3, ENTP2, ETV1, MEIS1, MEIS2, PAX6, PF2R, LMO3, SRTM1, SLIK6, THS7A, ZN506\tQ92482, P52757, Q96QS3, Q02363, P47928, Q16568, Q06889, Q9Y5L3, P50549, O00470, O14770, P26367, P43088, Q8TAP4, A2A2V5, Q9H5Y7, Q9UPZ6, Q5JVG8\tSingle-cell sequencing\t27693023\tNA\n")
cat("Human\tBrain\tUBERON_0000955\tNormal\tNormal cell\tAstrocyte\tCL_0000127\tGFAP\tGFAP\t2670\tGFAP\tP14136\tSingle-cell sequencing\t27587997\tNA\n")
cat("Human\tBrain\tUBERON_0000955\tNormal\tNormal cell\tNeuron\tCL_0000540\tENO2, SYP\tENO2, SYP\t2026, 6855\tENOG, SYPH\tP09104, P08247\tSingle-cell sequencing\t27587997\tNA\n")
cat("Human\tUndefined\tNA\tNormal\tNormal cell\tNeuronal progenitor cell\tNA\tMAP2, PAX6, SOX2\tMAP2, PAX6, SOX2\t4133, 5080, 6657\tMAP2, PAX6, SOX2\tP50579, P26367, P48431\tSingle-cell sequencing\t27534536\tNA\n")
sink()


# all_cell_markers.txt
sink(infile4)
cat("speciesType\ttissueType\tUberonOntologyID\tcancerType\tcellType\tcellName\tCellOntologyID\tcellMarker\tgeneSymbol\tgeneID\tproteinName\tproteinID\tmarkerResource\tPMID\tCompany\n")
cat("Human\tKidney\tUBERON_0002113\tNormal\tNormal cell\tProximal tubular cell\tNA\tIntestinal Alkaline Phosphatase\tALPI\t248\tPPBI\tP09923\tExperiment\t9263997\tNA\n")
cat("Human\tLiver\tUBERON_0002107\tNormal\tNormal cell\tIto cell (hepatic stellate cell)\tCL_0000632\tSynaptophysin\tSYP\t6855\tSYPH\tP08247\tExperiment\t10595912\tNA\n")
cat("Human\tEndometrium\tUBERON_0001295\tNormal\tNormal cell\tTrophoblast cell\tCL_0000351\tCEACAM1\tCEACAM1\t634\tCEAM1\tP13688\tExperiment\t10751340\tNA\n")
cat("Human\tGerm\tUBERON_0000923\tNormal\tNormal cell\tPrimordial germ cell\tCL_0000670\tVASA\tDDX4\t54514\tDDX4\tQ9NQI0\tExperiment\t10920202\tNA\n")
cat("Human\tCorneal epithelium\tUBERON_0001772\tNormal\tNormal cell\tEpithelial cell\tCL_0000066\tKLF6\tKLF6\t1316\tKLF6\tQ99612\tExperiment\t12407152\tNA\n")
cat("Human\tPlacenta\tUBERON_0001987\tNormal\tNormal cell\tCytotrophoblast\tCL_0000351\tFGF10\tFGF10\t2255\tFGF10\tO15520\tExperiment\t15950061\tNA\n")
cat("Human\tPeriosteum\tUBERON_0002515\tNormal\tNormal cell\tPeriosteum-derived progenitor cell\tNA\tCD166, CD45, CD9, CD90\tALCAM, PTPRC, CD9, THY1\t214, 5788, 928, 7070\tCD166, PTPRC, CD9, THY1\tQ13740, P08575, P21926, P04216\tExperiment\t15977065\tNA\n")
cat("Human\tAmniotic membrane\tUBERON_0009742\tNormal\tNormal cell\tAmnion epithelial cell\tCL_0002536\tNANOG, OCT3/4\tNANOG, POU5F1\t79923, 5460\tNANOG, PO5F1\tQ9H9S0, Q01860\tExperiment\t16081662NA\n")
cat("Human\tPrimitive streak\tUBERON_0004341\tNormal\tNormal cell\tPrimitive streak cell\tNA\tLHX1, MIXL1\tLHX1, MIXL1\t3975, 83881\tLHX1, MIXL1\tP48742, Q9H2W2\tExperiment\t16258519\tNA\n")
sink()


# Human_cell_markers.txt
cellMarkerToGmt(infile1, outfile1_1, uniq.column=c("tissueType"),
	geneid.type=c("geneID"))
cellMarkerToGmt(infile1, outfile1_2, uniq.column=c("tissueType"),
	geneid.type=c("geneSymbol"))
cellMarkerToGmt(infile1, outfile1_3, uniq.column=c("cellName"),
	geneid.type=c("geneID"))
cellMarkerToGmt(infile1, outfile1_4, uniq.column=c("cellName"),
	geneid.type=c("geneSymbol"))

# Mouse_cell_markers.txt
cellMarkerToGmt(infile2, outfile2_1, uniq.column=c("tissueType"),
	geneid.type=c("geneID"))
cellMarkerToGmt(infile2, outfile2_2, uniq.column=c("tissueType"),
	geneid.type=c("geneSymbol"))
cellMarkerToGmt(infile2, outfile2_3, uniq.column=c("cellName"),
	geneid.type=c("geneID"))
cellMarkerToGmt(infile2, outfile2_4, uniq.column=c("cellName"),
	geneid.type=c("geneSymbol"))

# Single_cell_markers.txt
cellMarkerToGmt(infile3, outfile3_1, uniq.column=c("tissueType"),
	geneid.type=c("geneID"))
cellMarkerToGmt(infile3, outfile3_2, uniq.column=c("tissueType"),
	geneid.type=c("geneSymbol"))
cellMarkerToGmt(infile3, outfile3_3, uniq.column=c("cellName"),
	geneid.type=c("geneID"))
cellMarkerToGmt(infile3, outfile3_4, uniq.column=c("cellName"),
	geneid.type=c("geneSymbol"))

# all_cell_markers.txt
cellMarkerToGmt(infile4, outfile4_1, uniq.column=c("tissueType"),
	geneid.type=c("geneID"))
cellMarkerToGmt(infile4, outfile4_2, uniq.column=c("tissueType"),
	geneid.type=c("geneSymbol"))
cellMarkerToGmt(infile4, outfile4_3, uniq.column=c("cellName"),
	geneid.type=c("geneID"))
cellMarkerToGmt(infile4, outfile4_4, uniq.column=c("cellName"),
	geneid.type=c("geneSymbol"))


# Human_cell_markers.txt
gmt1_1 <- getGmt(outfile1_1)
gmt1_2 <- getGmt(outfile1_2)
gmt1_3 <- getGmt(outfile1_3)
gmt1_4 <- getGmt(outfile1_4)

# Mouse_cell_markers.txt
gmt2_1 <- getGmt(outfile2_1)
gmt2_2 <- getGmt(outfile2_2)
gmt2_3 <- getGmt(outfile2_3)
gmt2_4 <- getGmt(outfile2_4)

# Single_cell_markers.txt
gmt3_1 <- getGmt(outfile3_1)
gmt3_2 <- getGmt(outfile3_2)
gmt3_3 <- getGmt(outfile3_3)
gmt3_4 <- getGmt(outfile3_4)

# all_cell_markers.txt
gmt4_1 <- getGmt(outfile4_1)
gmt4_2 <- getGmt(outfile4_2)
gmt4_3 <- getGmt(outfile4_3)
gmt4_4 <- getGmt(outfile4_4)

# Test
expect_true(length(gmt1_1) == 9)
expect_true(length(gmt1_2) == 9)
expect_true(length(gmt1_3) == 9)
expect_true(length(gmt1_4) == 9)

expect_true(length(gmt2_1) == 6)
expect_true(length(gmt2_2) == 6)
expect_true(length(gmt2_3) == 9)
expect_true(length(gmt2_4) == 9)

expect_true(length(gmt3_1) == 4)
expect_true(length(gmt3_2) == 4)
expect_true(length(gmt3_3) == 9)
expect_true(length(gmt3_4) == 9)

expect_true(length(gmt4_1) == 9)
expect_true(length(gmt4_2) == 9)
expect_true(length(gmt4_3) == 9)
expect_true(length(gmt4_4) == 9)
