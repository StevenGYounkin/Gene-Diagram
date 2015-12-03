snpplot <- function(genesnps=genesnps,protsnps=protsnps,wd=70,cex=0.7,start=0,end=0,genome="hg19",
                    domains=domains,exin=exin,domprop=domprop,col="royalblue1",intfact=10,
                    utrfact=1,interval=20,strand="+") {
  
  # load packages
  library(seqinr)
  library(Gviz)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(grid)
  
  # read in exons/introns and domains
  exin <- read.fasta(exin) #exin must be DNA-FASTA with alternating exons/introns separated by '>', last entry is counted as 3'UTR
  dom <- read.fasta(domains) #domains must be DNA-FASTA with domains separated by '>'
  
  # extract start/end from FASTA if not set
  ann <- getAnnot(exin)
  strand <- sub(".*?strand=(.*?) .*", "\\1", ann[[1]]) #extracts strand direction (+ or -) from FASTA
  if (strand=="+") {
    if (start==0) {
      start <- as.numeric(sub(".*?:(.*?)-.*", "\\1", ann[[1]])) } #extracts start of gene from FASTA 
    if (end==0) {
      end <- as.numeric(sub(".*?-(.*?) .*", "\\1", ann[[length(ann)]])) } #extracts end of gene from FASTA
  }
  if (strand=="-") {
    if (start==0) {
      start <- as.numeric(sub(".*?-(.*?) .*", "\\1", ann[[1]])) } #extracts start of gene from FASTA 
    if (end==0) {
      end <- as.numeric(sub(".*?:(.*?)-.*", "\\1", ann[[length(ann)]])) } #extracts end of gene from FASTA
  }
  
  # SNP position relative to start of gene
  snps <- abs(genesnps$bp - start) 
 
  
  ## chromosome
  #create variables needed for Gviz
  M <- genesnps # 'M' needed by IdeogramTrack in Gviz
  chrnum <- paste("chr",genesnps$Chr[1],sep="")
  #run ChromosomeTrack in Gviz
  chrTr <- IdeogramTrack(genome = genome, chromosome = chrnum, showFeatureId = "showLab.chr", cex = 2*cex)
  
  
  ## gene
  
  # calculate lengths of exons and introns
  exinlength <- vector()
  for (i in 1:length(exin)) {
    exinlength[i] <- length(exin[[i]]) 
  }
  
  # split into lengths of exons only and lengths of introns only
  exlength <- exinlength[1:length(exinlength)%%2 == 1]
  inlength <- exinlength[1:length(exinlength)%%2 == 0]
  
  # calculate start of each exon
  exst <- vector()
  exst[1] <- 1
  for (i in 2:length(exlength)) {
    exst[i] <- exst[i-1] + exlength[i-1] + inlength[i-1] 
  }
  exst[length(exlength)+1] <- abs(start-end)+1 # add end of gene as last exst entry
  
  # calculate exon number for each SNP
  snpexon <- vector()
  snpexon <- rep(0,length(snps))
  
  for (i in 1:length(snpexon)) {
    j <- 1
    while (snps[i]>=exst[j])  {
      snpexon[i] <- snpexon[i]+1
      j <- j+1  }
  }
  snpexes <- sort(unique(snpexon))
  
  # calculate SNP position within exon
  snpexpos <- vector()
  for (i in 1:length(snpexon)) {
    snpexpos[i] <- snps[i] - exst[snpexon[i]]
  }
  
  # calculate start of each exon and intron with adjusted proportions
  exstart <- vector()
  instart <- vector()
  exstart[1] <- 1
  for (i in 2:length(exlength)) {
    instart[i-1] <- exstart[i-1] + exlength[i-1] 
    exstart[i] <- exstart[i-1] + exlength[i-1] + (1/intfact)*inlength[i-1] 
  }
  instart[length(inlength)] <- exstart[length(inlength)] + exlength[length(inlength)]
  
  # calculate total length of gene with adjusted proportions of introns and 3'UTR
  geneend <- instart[length(inlength)] + 1/utrfact*inlength[length(inlength)] #adjust size of 3'UTR
  
  
  ## exons containing snps of interest
  
  # calculate end of displayed exons plus intervals
  length <- vector()
  length[1] <- 0
  h <- 1
  for (i in snpexes) {
    length[h+1] <- length[h] + exlength[i] + interval
    h <- h +1
  }
  length[length(length)] <- length[length(length)] - interval
  
  
  ## protein domains
  
  # calculate length of each domain
  domlength <- vector()
  for (i in 1:length(dom))
    domlength[i] <- length(dom[[i]])
  names(domlength) <- domprop[1,]
  
  # calculate end of each domain
  domend <- vector()
  domend[1] <- length(dom[[1]])
  for (i in 2:length(dom))
    domend[i] <- length(dom[[i]]) + domend[[i-1]]
  
  # calculate protein length
  protlength <- sum(domlength)
  
  # define colors
  domcol <- domprop[2,]
  
  ###plot
  
  # plot chromosome
  pdf('chrgenesnps.pdf') # save all figures in pdf
  par(no.readonly=TRUE)
  grid.newpage()
  pushViewport(viewport(height=0.07, y=0.5, just="bottom"))
  plotTracks(chrTr, from = start ,to = end, panel.only = TRUE, add=TRUE)
  popViewport()
  
  ## plot other figures
  par(no.readonly=TRUE)
  
  par(mfrow=c(3,1),mfg=c(3,1),mar=c(2,1,2,1))
  
  # initiate gene plot 
  plot(c(0, geneend), c(-250, 250), type= "n", xaxt='n', ylab="",xlab="",yaxt='n',frame.plot=FALSE)
  
  # plot exons and introns
  for (i in 1:(length(exlength)-1)) {
    rect(exstart[i], -wd, instart[i]-1, wd, col=col, border="black") #exons
    rect(instart[i], -0.2, exstart[i+1]-1, 0.2, col="grey", border="black") #introns
  } 
  rect(exstart[length(exlength)], -wd, instart[length(exlength)], wd, col=col, border="black") #last exon
  rect(instart[length(exlength)],-0.7*wd, geneend, 0.7*wd, col="lightskyblue", border="black" ) #3'UTR
  
  # plot SNPs
  for (i in 1:(length(snps))) {
    rect(exstart[snpexon[i]]+snpexpos[i]-0.1, -wd, exstart[snpexon[i]]+snpexpos[i]+0.1, wd+20, col="black", border="black") 
  }
  # plot SNP names
  for (i in 1:(length(snps))) {
    text(exstart[snpexon[i]]+snpexpos[i],wd+30, labels = genesnps$id[i], cex=0.9*cex, srt = 45, pos =4, offset=0.1) 
  }
  text(instart[length(inlength)],-0.5*wd, labels = "3'UTR", cex=0.9*cex, srt = 90, pos =4) #plot 3'UTR
  
  
  # initiate exons plot
  plot(c(length[1], length[length(length)]), c(-250, 250), xaxt='n', ylab="", xlab="", yaxt='n', type= "n",frame.plot=FALSE)
  
  # plot exons 
  for (i in 1:length(snpexes)) {
    rect(length[i], -wd, length[i] + exlength[snpexes[i]], wd, col=col, border="black") #exons
  }
  
  # plot SNPs
  for (i in 1:(length(snps))) {
    rect(length[snpexes==snpexon[i]]+snpexpos[i]-0.01, -wd, length[snpexes==snpexon[i]]+snpexpos[i]+0.01, wd + 20, col="black", border="black") }
  # plot SNP names
  for (i in 1:(length(snps))) {
    text(length[snpexes==snpexon[i]]+snpexpos[i],wd +30, labels = genesnps$id[i], cex=0.9*cex, srt = 45, pos =4,offset=0.1) }
  # plot exon names
  exonnum <- paste("exon ", snpexes,"/",length(exlength),sep="")
  for (i in 1:(length(snpexes))) {
    text(length[i],-0.5*wd, labels = exonnum[i], cex=1.5*cex, pos =4) }
  
  
  # initiate protein plot
  plot(c(0, protlength), c(-250, 250), type= "n", xaxt='n', ylab="",xlab="",yaxt='n',frame.plot=FALSE)
  
  # plot domains in different colors
  rect(0, -wd, domlength[i], wd, col=domcol[1], border="black") 
  for (i in 2:length(dom))
    rect(domend[i-1], -wd, domend[i], wd, col=domcol[i], border="black") 
  
  # plot SNPs 
  for (i in 1:(length(snps))) {
    rect(3*protsnps[i]-0.1, -wd, 3*protsnps[i]+0.1, wd + 20, col="black", border="black") 
  }
  # plot domain names
  text(domend[1]/2,- wd - 10, labels = names(domlength[1]), cex=1*cex, pos=1, col=domcol[1])
  for (i in 2:(length(domend))) {
    text(domend[i-1]+(domend[i]-domend[i-1])/2, -wd - 10, labels = names(domlength[i]), cex=1*cex, pos=1, col=domcol[i]) 
  }
  # plot SNP names
  for (i in 1:(length(snps))) {
    text(3*protsnps[i],wd + 30, labels = genesnps$id[i], cex=0.9*cex, srt = 45, pos =4,offset=0.1) }
  
  dev.off()
}

#TLR5
#genesnps <- data.frame(Chr=c(1,1,1,1),id=c("rs56243703","rs5744174","rs5744175","rs5744177"),bp=c(223284069,223284528,223284444,223283837))
#domprop <- rbind(c("SS","LRRNT","LRR 1-22","LRRCT","TM","TIR"),c("red4","orange","yellow","green","cyan","blue"))
#protsnps <- c(2303/3,1844/3,1928/3,2535/3)
#domains <- "TLR5 domains complete.txt"
#exin <- "TLR5genome.txt"

#snpplot(genesnps=genesnps,exin=exin,protsnps=protsnps,domprop=domprop,domains=domains)

#SorCS2
genesnps <- data.frame(Chr=c(4,4,4,4),id=c("rs34058821","rs16840892","rs35935435","rs6816604"),bp=c(7666160,7716061,7666180,7716932))
domprop <- rbind(c("SS","Propeptide","VPS10P","PKD","LRR","TM","C-tail"),c("red4","orange","yellow","green","cyan","blue","magenta4"))
protsnps <- c(345,351,695,716)
domains <- "SorCS2 domains.txt"
exin <- "SorCS2 EXON intron 3'UTR.txt"

snpplot(genesnps=genesnps,protsnps=protsnps,domains=domains,exin=exin,domprop=domprop,intfact=100,
        utrfact=10)