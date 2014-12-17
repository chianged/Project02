#command that takes a primer or pair of primers and tells you how good they are with respect to the taxon they target


data <- scan(file="sample.database", what="list",sep=">")     #Import database split into 3's: (1) -nothing-, (2) descriptors/taxa, (3) sequence
#1 + 3x = -nothing-
#2 + 3x = descriptors/taxa
#3x = sequence

primer.quality <- function(primer,taxa){     #primer = primer sequence, taxa = taxon targeted by primer
  rev.primer <-  sapply(lapply(strsplit(primer, NULL), rev), paste, collapse="")     #Reverse sequence
  comp.primer <- chartr("ATCG","TAGC",rev.primer)     #Switches base pairs to complement
  prime.seq <- grep(comp.primer,data)     #Finds when the primer sequence appears in the database; matches primer with sequence
  tax <- grep(taxa,data)     #Finds when taxon appears in database
  l.prime <- length(prime.seq)     #Sum primer sequence appearance
  l.tax <- length(tax)     #Sum taxon appearance
  taxa.seq <- prime.seq - 1     #Every group comes in threes: (1) -nothing-, (2) Taxa info, (3) Sequence
  #If I have the sequence index, I can subtract by 1 to access the taxa/descriptor index
  new.taxa <- data[taxa.seq]     #Index the descriptor/taxa corresponding to the sequence matches
  primer.target <- grep(taxa,new.taxa)     #Which of the descriptor/taxa actually have the taxon of interest
  l.target <- length(primer.target)     #Number of descriptor/taxa that have taxon of interest
  nonspecific <- 1 - (l.target / l.prime)     #What proportion of sequence matches are nonspecific
  l.nonspecific <- l.prime - l.target     #Calculate number of nonspecific matches
  l.prime <- l.prime - l.nonspecific      #Update number of primer matches by removing the number of nonspecific matches
  target.per <- l.prime/l.tax     #Divide sequence by taxon; what % of taxon is targeted
  print(target.per)     #Proportion of target taxon that are bound
  print(nonspecific)     #Proportion of primer binding that is nonspecific
}

##Reference: rev.primer command from http://stackoverflow.com/questions/13612967/how-to-reverse-a-string-in-r
##chartr function from http://stackoverflow.com/questions/20371854/complement-a-dna-sequence-in-r





################################First attempt that took waaaaaaayyyy too long to run#########################################
# #####Import Silva database#####
# data <- scan("GreenGenes.database",what="")
# 
# 
# #####Create matrix of taxa sequences#####
# 
# a <- grepl("A",data)     #TRUE or FALSE for letter A
# t <- grepl("T",data)     #TRUE or FALSE for letter T
# c <- grepl("C",data)     #TRUE or FALSE for letter C
# g <- grepl("G",data)     #TRUE or FALSE for letter G
# 
# a.true <- which(a)     #Which spots contain letter A
# t.true <- which(t)     #Which spots contain letter A
# c.true <- which(c)     #Which spots contain letter A
# g.true <- which(g)     #Which spots contain letter A
# 
# if(which(a)==which())
# 
# which(c(a,t))
# 
# seqs <- matrix(rep("NA"),ncol=1,nrow=870628)     #Create a matrix of NA's
# 
# 
# 
# 
# 
# seq.spots <- rep(NA,3)     #Create list of spots where seq's are located
# 
# 
#   a.val <- sum(a.true==x)     #Is there an A in spot x?
#   t.val <- sum(t.true==x)     #Is there an T in spot x?
#   c.val <- sum(c.true==x)     #Is there an C in spot x?
#   g.val <- sum(g.true==x)     #Is there an G in spot x?
#   if(a.val==1 & t.val==1 & c.val==1 & g.val==1){     #If all bases exist
#     b <- b+1
#     seqs[b]<- data[x]     #Fill in the matrix with that sequence
#     seq.spots[b] <- x     #Fill in list with spot location of seq's
#   }
# }
# 
# 
# 
# 
# 
# 
# 
# x=0     #Refers to spot in string to evaluate
# b=0     #Refers to spot in matrix to fill in
# seq.spots <- rep(NA,3)     #Create list of spots where seq's are located
# for (i in data[1:8017354]){
#   x <- x+1
#   a.val <- sum(a.true==x)     #Is there an A in spot x?
#   t.val <- sum(t.true==x)     #Is there an T in spot x?
#   c.val <- sum(c.true==x)     #Is there an C in spot x?
#   g.val <- sum(g.true==x)     #Is there an G in spot x?
#   if(a.val==1 & t.val==1 & c.val==1 & g.val==1){     #If all bases exist
#     b <- b+1
#     seqs[b]<- data[x]     #Fill in the matrix with that sequence
#     seq.spots[b] <- x     #Fill in list with spot location of seq's
#   }
# }
# 
# 
# #####Create matrix of taxonomy#####
# 
# ##Fill in matrix by row with taxonomy
# ##If spot # > than spot # for true ATCG in previous part, move to next row
# 
# 
# underscore <- grepl("_",data)     #TRUE or FALSE for underscore
# 
# underscore.true <- which(underscore)     #Which spots contain underscore
# 
# tax <- matrix(rep("NA"),ncol=8,nrow=10)     #Create a matrix of NA's
# 
# y=0     #Refers to spot in string to evaluate
# e=0     #Refers to spot in matrix to fill in
# f=1     #
# for (i in data[1:50]){
#   y <- y+1
#   underscore.val <- sum(underscore.true==y)     #Is there an underscore in spot y?
#   underscore.sum <- sum(y > seq.spots)     #Is the taxonomy located before the seq? (All taxa belongs to the first following sequence)
#   if(underscore.val==1 & underscore.sum<f){     #If there's an underscore that's before a seq
#     e <- e+1     #Matrix column to fill
#     tax[f,e] <- data[y]     #Fill in the matrix with that sequence
#   } else if(underscore.val==1 & underscore.sum==f){     #If there's an underscore that's after a seq
#     f <- f+1     #Move to the next row
#     e <- 1     #Reset column
#     tax[f,e] <- data[y]     #Fill in matrix
#   } else {     #For the icky info that I don't need
#     
#   }
# }
# 
# 
# #####Combine sequence and taxonomy tables in one powerhouse table#####
# tax[,8] <- seqs
# gg.data <- as.data.frame(tax)