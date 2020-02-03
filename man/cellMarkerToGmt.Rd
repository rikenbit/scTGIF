\name{cellMarkerToGmt}
\alias{cellMarkerToGmt}
\docType{methods}
\title{
A function to convert the CellMarker data to GMT files.
}
\description{
The GMT (Gene Matrix Transposed file format : *.gmt) file is formatted by the Broad Institute (https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29). The data can be downloaded from the website of CellMarker (http://biocc.hrbmu.edu.cn/CellMarker).
}
\usage{
cellMarkerToGmt(infile, outfile,
  uniq.column=c("tissueType", "cellName"),
  geneid.type=c("geneID", "geneSymbol"))
}
\arguments{
\item{infile}{
The input file downloaded from CellMarker website
}
\item{outfile}{
The output GMT file converted from the CellMarker data
}
\item{uniq.column}{
The duplicated terms in the specified column are aggrgated as a row of GMT file (Default: geneID)
}
\item{geneid.type}{
Output gene identifier. (Default: geneID)
}
}
\value{
\item{output}{
A GMT file is generated.
}
}
\author{
\packageAuthor{scTGIF}
}
\examples{
  library("GSEABase")

  tmp <- tempdir()
  infile1 = paste0(tmp, "/Human_cell_markers.txt")
  outfile1_1 = paste0(tmp, "/Human_cell_markers_1.gmt")
  outfile1_2 = paste0(tmp, "/Human_cell_markers_2.gmt")
  outfile1_3 = paste0(tmp, "/Human_cell_markers_3.gmt")
  outfile1_4 = paste0(tmp, "/Human_cell_markers_4.gmt")

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

  cellMarkerToGmt(infile1, outfile1_1, uniq.column=c("tissueType"),
    geneid.type=c("geneID"))
  cellMarkerToGmt(infile1, outfile1_2, uniq.column=c("tissueType"),
    geneid.type=c("geneSymbol"))
  cellMarkerToGmt(infile1, outfile1_3, uniq.column=c("cellName"),
    geneid.type=c("geneID"))
  cellMarkerToGmt(infile1, outfile1_4, uniq.column=c("cellName"),
    geneid.type=c("geneSymbol"))

  gmt1_1 <- getGmt(outfile1_1)
  gmt1_2 <- getGmt(outfile1_2)
  gmt1_3 <- getGmt(outfile1_3)
  gmt1_4 <- getGmt(outfile1_4)
}

\keyword{methods}