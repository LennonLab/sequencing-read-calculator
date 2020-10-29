
library(ggplot2)

# Scenario: resequencing an evolved population (single species)
# to charecterize diversity. 

# The following code calculates the number of reads needed
# to insure that all the bases are sequenced by some minimal covrage.

# Calculation is based on footnote #2 in Ilumina techical note 
# "Estimating Sequencing Coverage"
# https://support.illumina.com/content/dam/illumina-marketing/documents/products/technotes/technote_optimizing_coverage_for_targeted_resequencing.pdf
# see pdf file in repo
###############################
# Alter the values of varibles y, L and G 
# to match specific system of intrest.

# The minimal coverage per base
  #  Illumina: "y is the number of times a base is read"
  # in ppois() y is "q"

y=100 

# Read length, based on specifics of sequencing protocol   
  # Illumina: "L is the read length".
    # If paired end enter the sum of the reads (eg. if sequencing 2X150 enter 300)
L=75 # unpaired reads of 75 nt

# Genome size
    # Illumina: "G is the haploid genome length (bp)"
G=132562 #SPO1 genome length

################################

# N is the number of reads
# We strart with a very broad range and will refine in a second step
N=10^seq(4,14,0.1)

# C stands for coverage
#This is the average coverage:
C=L*N/G

  # in ppois() C is "lambda"
# Number of reads needed aprox.
N.needed <- N[1+sum(ppois(y,C)>1/(G+1))]


# Refine the range of N based on aprox.  N.needed
N=seq(from = N.needed/2, to =  N.needed*2,length.out = 200)

# calculate coverage for refined range of N
C=L*N/G

# Number of reads needed refined
N.needed <- N[1+sum(ppois(y,C)>1/(G+1))]

d=data.frame(x=log10(N),y=ppois(y,C))
d <- d[d$y>1e-20,]
#plot 
ggplot(d, aes(x,y))+
  geom_point()+
  #draw line indicatng the probabilty that every base will be sequenced at coverage C
  geom_hline(yintercept = 1/(G+1))+
  #draw line indicatng the number of reads needed
  geom_vline(xintercept = log10(N.needed))+
  scale_y_log10()+
  xlab("Log10 of number of reads")+
  ylab(paste0("P(min. coverage ","\u2265", y,")"))+
  labs(title= paste("Reads needed",round(N.needed),"= 10 ^",round(log10(N.needed),2)),
       subtitle = paste("Average coverage = ",round(N.needed*L/G,2)),
       caption = "horizontal line indicates the frequency lower than a single bp\n(1/genome size+1)")

# output to console
paste("Reads needed",round(N.needed))


