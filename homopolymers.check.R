library(Biostrings)

# Change working dir
#setwd("~/storage/QIIME_DATABASES/")
setwd(".")

# Load sequence multifasta file as object

#file<-file.path("~/storage//QIIME_DATABASES//maarjAM//maarjAM.5.fasta")
#file<-file.path("~/storage//QIIME_DATABASES//sh_qiime_release_s_09.02.2014/sh_refs_qiime_ver6_dynamic_s_09.02.2014.fasta")

#file<-file.path("~/DATA/RISINNOVA-18S//risinnova_AMF_soil_all_libraries.fasta")
#file<-file.path("~/DATA/RISINNOVA-18S//risinnova_AMF_soil_denoised_inflated_filtered.fasta")

#file<-file.path("~/DATA/RISINNOVA-ITS//risinnova_ITS_all_libraries.fasta")
#file<-file.path("~/DATA/RISINNOVA-ITS//risinnova_ITS_left_inflated_filtered.fasta")
#file<-file.path("~/DATA/RISINNOVA-ITS//risinnova_ITS_right_inflated_filtered.fasta")
file<-file.path("~/storage//ORCHIDEE//ANALYSIS//GZL-1/GZL-1.pear28.fastq.assembled.fna")

seqsObject <- readDNAStringSet(file, "fasta")
basename(file)

# Convert object into vector
seqsVector<-toString(seqsObject)

# Declare some variables
col_name=NULL
myList <- list() 


for (nucl in c("A","T","G","C")) {
  
  # Calculate max homopolymer extension
  lC_nucl=longestConsecutive(seqsVector, nucl)

  summa_temp=NULL
  
  for (n in 7:lC_nucl)
  {
    pattern = paste0(rep(nucl,times= n), collapse="")
    
    temp =sum(vcountPattern(pattern, seqsObject,
                            max.mismatch=0, min.mismatch=0,
                            with.indels=FALSE, fixed=TRUE,
                            algorithm="auto"))
    
    summa_temp = c(summa_temp,temp)
    col_name<-c(col_name,n)
    myList[[paste("vector_",nucl,sep="")]] <- summa_temp
      }
}

myList

col_name.f<-factor(col_name)

# Define colors to be used for A, T, G, C
plot_colors <- rainbow("4")

# Calculate range from 0 to max value of count
g_range <- range(0, myList)

# Graph first vecto using rainbow points overlayed by a line
plot(myList$vector_A, type="b", pch=22, lty=2, col=plot_colors[1], ylim=g_range, xaxt="n", ann=FALSE)

# Graph other vectors with dashed line and square points
lines(myList$vector_T, type="b", pch=22, lty=2, col=plot_colors[2])
lines(myList$vector_G, type="b", pch=22, lty=2, col=plot_colors[3])
lines(myList$vector_C, type="b", pch=22, lty=2, col=plot_colors[4])

# Make x axis
axis(1, at=1:length(levels(col_name.f)), labels=paste0(levels(col_name.f),"X"))

# Create a title with a red, bold/italic font
title(main="Homopolymers Distribution", col.main="red", font.main=4)
#title(main=, col.main="red", font.main=3)
mtext(paste0("\n\n",basename(file)), side=3, col="red", cex=.9, line = 0.5)

# Label the x and y axes with dark green text
title(xlab="Homopolymers Length", col.lab=rgb(0,0.5,0))
title(ylab="Counts", col.lab=rgb(0,0.5,0))

# Create a legend that usesthe same line colors and points 
# used by the actual plots
#x<-length(levels(col_name.f))-1
legend("topright", c("A","T","G","C"), cex=.9, col=c(plot_colors[1],plot_colors[2],plot_colors[3],plot_colors[4]), bty='n', pch=22, lty=2)
