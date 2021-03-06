# Primer Quality
Edna Chiang  
December 16, 2014  
  
<br>
      
      
###primer.quality  
This command assesses primer quality in regards to how well a 16S rRNA primer/probe matches its target taxon. It specifically uses the Greengenes 16S rRNA database. The database is first imported & split based upon descriptors (such as taxonomy) and 16S rRNA sequence. The primer.quality command receives 2 inputs: primer sequence and taxon of interest. The command searches the database for matches to the two inputs, and calculates primer quality by comparing the number of sequence matches to taxon matches. The outputs are the proportion of targeted taxon matches, and the proportion of non-targeted taxon sequence matches.
  
  <br>
      
###Import Database  
The Greengenes 16S rRNA database can be downloaded at [http://greengenes.lbl.gov/Download/Sequence_Data/Fasta_data_files/current_GREENGENES_gg16S_unaligned.fasta.gz]  (link provided by [Human Oral Microbiome Database](http://www.homd.org/index.php?name=Article&file=index&sid=28)). The database was originally to be provided in this repository, but the file is 5.6 MB too big. Sad face. Instead, I'm using a subset of the database called sample.database.  
The next step is to organize/split the database:

```r
data <- scan(file="sample.database", what="list",sep=">")     #Import database split into 3's: (1) -nothing, (2) descriptors/taxa, (3) sequence
  #1 + 3x = -nothing-
  #2 + 3x = descriptors/taxa
  #3x = sequence
```
This splits the database into a character list organized as the following clusters:   

Seat (Number) | Content
------------- | -------------
1 + 3x        | ""           (-no information-)
2 + 3x        | Descriptors/taxonomy information
3x            | Sequence  

Each cluster of 3 seats correspond to each other; the descriptors/taxonomy information refers to the sequence, and vice versa.
  
  <br>
      
###Load Command  
Once the database has been prepared, you can load the primer.quality command:  

```r
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
```
**Inputs:**  
<ul>
<li>primer = primer sequence
<ul>
<li>Must be 5' -> 3'</li>
<li>Include only A, T, C, and G's (other IUPAC nucleotide codes will not recognized)</li>
</ul></li>
<li>taxa = taxon of interest
<ul>
<li>Can be kingdom, phylum, class, orer, family, genus, or species/OTU</li>
</ul></li>
</ul>
  
  
  <br>
    
**Outputs:**  
<ol>
<li>Proportion of target taxon to which primer binds</li> 
<li>Proportion of primer binding that is to non-targeted taxa</li>
</ol>
  
  
  
<br>
  
    
**What's happening in the command?**  
1. Calculate the complement of the primer sequence  
2. Determine where and how many times the following appear in the database: Primer complement sequence, taxon of interest  
3. Calculate binding to non-target taxa  
4. Subtract non-targeted binding from total primer binding  
5. Calculate the proportion of targeted taxon which is bound by primer  

<br>
    
      
###Example  
Let's test out the command!   
For this example, we'll use the probe EUB338 that targets Bacteria and whose sequence is "GCTGCCTCCCGTAGGAGT."  
Its sequence is can be found at [probeBase](http://www.microbial-ecology.net/probebase/search.asp).

```r
primer.quality("GCTGCCTCCCGTAGGAGT","Bacteria")
```

```
## [1] 0.6894
## [1] 0
```


 <br>  



###Limitations  
A primer's quality does not depend on only the percentage of the target taxon which it targets; there are many other factors that come into play such as random nonspecific binding and CG-content (if the primer is to be used for PCR). The primer.quality command does not account for these additional influences and requires further development in order to become a better determinant of primer quality. Until then, more stringent evaluations of primer quality can be conducted at [Silva](http://www.arb-silva.de/), [Greengenes](http://greengenes.lbl.gov/cgi-bin/nph-index.cgi), [probeCheck](http://131.130.66.200/cgi-bin/probecheck/content.pl?id=home), and [RDP](http://rdp.cme.msu.edu/probematch/search.jsp)

(Reference: [probeBase](http://131.130.66.201/probebase/))
