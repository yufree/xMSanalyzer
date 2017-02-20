\name{feat.batch.annotation.KEGG}
\alias{feat.batch.annotation.KEGG}
\title{feat.batch.annotation.KEGG}
\description{Batch annotation of features using KEGG REST for one or more adducts.}
\usage{feat.batch.annotation.KEGG(dataA, max.mz.diff=10, queryadductlist=c("M+H"),
xMSanalyzer.outloc,numnodes=1,syssleep=1)}
\arguments{
  \item{dataA}{Feature table from apLCMS or XCMS. The first column should be m/z.}
  \item{max.mz.diff}{Metlin matching m/z threshold in ppm. eg: 5}
  \item{queryadductlist}{List of adducts to be used for searching. eg: c("M+H","M+Na","M+K"),
  c("positive") for all positive ion adducts, or c("negative") for all negative ion adducts, or 
  c("all") for all adducts as defined in Metlin.}
  \item{xMSanalyzer.outloc}{Output location where the HTML and text reports will be written. 
  eg: "C:/experiment1/xMSanalyzeroutput/"}
  \item{numnodes}{Number of computing nodes to use. Default: 1}
  \item{syssleep}{Wait time (seconds) in between queries. Default: 1}
}
\details{
This utility uses the readHTMLTable() function in the XML package in R and the KEGG REST 
interface to get the list of compounds and pathways IDs from KEGG.
respectively. The output is generated as an HTML report and a text file that includes
pathway and compound annotations with links to external databases such as 
KEGG Compound, KEGG Pathway, PubChem Substance, HMDB,
ChEBI, CAS, and LipidMAPS. METLIN is no longer supported due to their updated
terms of use. The function takes as input a data frame with a list of
input m/z, a user-defined m/z threshold (ppm) to define the minimum and maximum
mass range, list of adducts; eg: c(“M+H”, “M+H-H2O”), and the output folder
location. A sample annotation report is available at the software homepage:
SampleAnnotation.html
}
\value{A list is returned.  \item{html.res }{Annotation report in HTML format}
\item{text.res }{Text delimited annotation report}
}
\author{Karan Uppal <kuppal2@emory.edu>}
\keyword{ ~annotation }
\keyword{ ~KEGG }
\keyword{ ~pathway }