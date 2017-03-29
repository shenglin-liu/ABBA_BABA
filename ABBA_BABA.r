#####################################################################################
# Functions for conducting ABBA BABA test. By Shenglin Liu, Oct 12, 2016.
# The inputs are BStringSet objects.
#	The characters in the sequences can be arbitrary.
#	But you need to tell the functions what character you use for missing values.
#####################################################################################

#####################################################################################
# References:
#	Durand EY, Patterson N, Reich D, Slatkin M (2011) Testing for ancient
#		admixture between closely related populations. Mol. Biol. Evol. 28
#		(8): 2239-2252.
#	Martin SH, Davey JW, Jiggins CD (2014) Evaluating the use of ABBA-BABA
#		statistics to locate introgressed loci. Mol. Biol. Evol. 32(1):
#		244-257.
#	Blackmon H, Adams RH (2015) evobiR: Comparative and
#		Population Genetic Analyses. R package version 1.1.
#		http://CRAN.R-project.org/package=evobiR
#####################################################################################

#####################################################################################
# Function list:
#	vcf2ab		: Generate input file for ABBA BABA test from small VCF file.
#	vcf2ab.big	: Generate input file for ABBA BABA test from big VCF files.
#	vcf2ab.core	: Core function for vcf2ab.big.
#	DSeq		: Calculates sequence-based D statistics.
#	DGeno		: Calculates genotype-based D statistics.
#	abba.baba	: Run ABBA BABA test.
#	introgressed	: Calculate introgression rate.
#	fSeq		: Calculate sequence-based fHom from P3 to P2.
#	fGeno		: Calculate genotype-based fHom from P3 to P2.
#	abba.baba2	: My version of ABBA BABA test.
#	slidingWindow	: Run ABBA BABA tests in sliding windows.
#	popBatch	: Run a batch of ABBA BABA tests for population data.
#####################################################################################

# Required libraries.
library(Biostrings)

#######=======================================================================#######
# Generate input file for ABBA BABA test from VCF file.
# vcf2ab is a simplified version of the conversion function.
#	It is easy to use, and is suitable when the VCF file is small
#		e.g., #SNP * #indi < 1e7.
# vcf2ab.big is the complete version of the conversion function.
#	It is a wrapper function of vcf2ab.core.
#	It is efficient when the VCF file is too big
#		or when multiple VCF files are involved.
# vcf2ab.core is the core function of vcf2ab.big.
#	It runs the main conversion algorithm.

# f.vcf: name of the VCF file; first 9 columns are information columns.
# f.popmap: name of the popmap file; two columns;
#	the first column contain names of individuals in the VCF file;
#	the second column specifies the population affinity of each individual.
# f.output: name of the output file.
# randomize: randomly assign the two alleles of a genotype to two chromosomes;
#	applicable to unphased data.
vcf2ab<-function(f.vcf,f.popmap,f.output,randomize=TRUE)
{
	vcf<-read.table(f.vcf,sep="\t",stringsAsFactors=FALSE)
	#refalt<-vcf[,4:5]
	vcf<-as.matrix(vcf[,-c(1:9)])
	nrow.vcf<-nrow(vcf)
	ncol.vcf<-ncol(vcf)
	chrom1<-as.vector(substring(vcf,1,1))
	chrom2<-as.vector(substring(vcf,3,3))
	if(randomize)
	{
		chrom<-rbind(chrom1,chrom2)
		index<-sample(c(TRUE,FALSE),nrow.vcf*ncol.vcf,replace=TRUE)
		index1<-2*(1:(nrow.vcf*ncol.vcf))-1+index
		index2<-2*(1:(nrow.vcf*ncol.vcf))-1+(!index)
		chrom1<-chrom[index1]
		chrom2<-chrom[index2]
	}
	chrom1<-matrix(chrom1,nrow.vcf,ncol.vcf)
	chrom2<-matrix(chrom2,nrow.vcf,ncol.vcf)
	
	popmap<-read.table(f.popmap,sep="\t",stringsAsFactors=FALSE)[,1]
	sink(f.output)
	for(i in 1:length(popmap))
	{
		cat(">",popmap[i],"_A","\n",sep="")
		cat(chrom1[,i],"\n",sep="")
		cat(">",popmap[i],"_B","\n",sep="")
		cat(chrom2[,i],"\n",sep="")
	}
	sink()
}

# fs.vcf: names of the VCF files to be handled.
# interval: the number of lines to be read in each time from the input stream.
vcf2ab.big<-function(fs.vcf,f.popmap,f.output,randomize=TRUE,interval=1e5)
{
	popmap<-read.table(f.popmap,sep="\t",stringsAsFactors=FALSE)[,1]
	fs.temp<-rbind(
		paste("abtemp_",popmap,"_A",sep=""),
		paste("abtemp_",popmap,"_B",sep="")
		)
	temp<-sapply(1:ncol(fs.temp),function(x)
		cat(">",popmap[x],"_A","\n",file=fs.temp[1,x],sep="")
		)
	temp<-sapply(1:ncol(fs.temp),function(x)
		cat(">",popmap[x],"_B","\n",file=fs.temp[2,x],sep="")
		)
	for(f.vcf in fs.vcf)
	{
		c.vcf<-file(f.vcf,"r")
		repeat
		{
			vcf<-try(read.table(c.vcf,sep="\t",stringsAsFactors=FALSE,
				nrows=interval),silent=TRUE)
			if(class(vcf)=="try-error")
			{
				break
			}
			vcf2ab.core(vcf,fs.temp,randomize)
		}
		close(c.vcf)
	}
	c.output<-file(f.output,"w")
	for(f.temp in fs.temp)
	{
		x<-readLines(f.temp)
		writeLines(x,c.output)
	}
	close(c.output)
	file.remove(fs.temp)
}

vcf2ab.core<-function(vcf,fs.temp,randomize)
{
	vcf<-as.matrix(vcf[,-c(1:9)])
	nrow.vcf<-nrow(vcf)
	ncol.vcf<-ncol(vcf)
	chrom1<-as.vector(substring(vcf,1,1))
	chrom2<-as.vector(substring(vcf,3,3))
	if(randomize)
	{
		chrom<-rbind(chrom1,chrom2)
		index<-sample(c(TRUE,FALSE),nrow.vcf*ncol.vcf,replace=TRUE)
		index1<-2*(1:(nrow.vcf*ncol.vcf))-1+index
		index2<-2*(1:(nrow.vcf*ncol.vcf))-1+(!index)
		chrom1<-chrom[index1]
		chrom2<-chrom[index2]
	}
	chrom1<-matrix(chrom1,nrow.vcf,ncol.vcf)
	chrom2<-matrix(chrom2,nrow.vcf,ncol.vcf)
	
	temp<-sapply(1:ncol(fs.temp),function(x)
		cat(chrom1[,x],file=fs.temp[1,x],sep="",append=TRUE)
		)
	temp<-sapply(1:ncol(fs.temp),function(x)
		cat(chrom2[,x],file=fs.temp[2,x],sep="",append=TRUE)
		)
}

#######=======================================================================#######
# Conduct ABBA BABA test according to Durand et al. (2011).
# DSeq calculates sequence-based D statistics.
# DGeno calculates genotype-based D statistics.
# abba.baba runs ABBA BABA test.
# introgressed calculates introgression rate.

# a: matrixed sequence alignment; individuals as rows.
DSeq<-function(a)
{
	abba<-(a[1,]==a[4,])&(a[2,]==a[3,])&(a[1,]!=a[2,])
	abba<-sum(abba)
	baba<-(a[1,]==a[3,])&(a[2,]==a[4,])&(a[1,]!=a[2,])
	baba<-sum(baba)
	d<-(abba-baba)/(abba+baba)
	if(is.nan(d))d<-0
	data.frame(D=d,ABBA=abba,BABA=baba)
}

# a: allele frequency matrix; populations as rows.
DGeno<-function(a)
{
	abba<-(1-a[1,])*a[2,]*a[3,]*(1-a[4,])+a[1,]*(1-a[2,])*(1-a[3,])*a[4,]
	abba<-sum(abba)
	baba<-a[1,]*(1-a[2,])*a[3,]*(1-a[4,])+(1-a[1,])*a[2,]*(1-a[3,])*a[4,]
	baba<-sum(baba)
	d<-(abba-baba)/(abba+baba)
	if(is.nan(d))d<-0
	data.frame(D=d,ABBA=abba,BABA=baba)
}

# a: sequence alignment, all sequences with the same length;
#	BStringSet object of at least 4 elements;
#	elements are arranged according to P1, P2, P3, and O;
#	see Durand et al. (2011) for the meanings of P1, P2, P3, and O.
# genotype: logic;
#	if FALSE (default), D statistics calculation is sequence-based;
#	if TRUE, D statistics calculation is genotype-based.
# ns.seq: 4 integers; number of sequences from P1, P2, P3, and O when genotype=TRUE.
# missing.char: character for missing values in a.
# test: possible values are "n" (no test), "b" (bootstrap), and "j" (jackknife).
abba.baba<-function(a,genotype=FALSE,ns.seq=c(1,1,1,1),
	missing.char=".",test="n",repeats=1000,block.size=1000)
{
	a<-as.matrix(a)
	fun<-DSeq
	
	# Filter for missing and non-biallelic (remove).
	Original<-ncol(a)
	index<-colSums(a==missing.char)==0
	a<-a[,index]
	index<-sapply(1:ncol(a),function(x)length(unique(a[,x]))<=2)
	a<-a[,index]
	AfterFilter<-ncol(a)
	SequenceLen<-data.frame(Original,AfterFilter)
	
	if(genotype)
	{
		# Calculate allele frequency.
		reference<-a[1,]
		freq<-numeric(0)
		j<-0
		for(i in ns.seq)
		{
			temp<-t(matrix(a[1:i+j,],i,AfterFilter))
			temp<-rowSums(temp!=reference)/ncol(temp)
			freq<-rbind(freq,temp)
			j<-j+i
		}
		a<-freq
		fun<-DGeno
	}
	
	# Calculate D statistics.
	DStats<-fun(a)
	
	# Do not implement significance test.
	if(test=="n")
	{
		b<-list(SequenceLen,DStats)
		names(b)<-c("SequenceLen","DStats")
	}
	
	# Test significance using bootstrap.
	if(test=="b")
	{
		bootstrap<-sapply(1:repeats,function(x)
			fun(a[,sample(1:AfterFilter,AfterFilter,replace=TRUE)])$D[1]
			)
		SD<-sd(bootstrap)
		Z<-abs(DStats$D[1]/SD)
		P.val<-2*(1-pnorm(Z))	#######!!!
		Bootstrap<-data.frame(n=repeats,SD,Z,P.val)
		b<-list(SequenceLen,DStats,Bootstrap)
		names(b)<-c("SequenceLen","DStats","Bootstrap")
	}
	
	# Test significance using jackknife.
	if(test=="j")
	{
		max.rep<-AfterFilter-block.size+1
		drop.pos<-as.integer(seq(from=1,to=max.rep,length.out=repeats))
		jackknife<-sapply(1:repeats,function(x)
			fun(a[,-c(drop.pos[x]:(drop.pos[x]+block.size-1))])$D[1]
			)
		SD<-sd(jackknife)
		Z<-abs(DStats$D[1]/SD)
		P.val<-2*(1-pnorm(Z))	#######!!!
		Jackknife<-data.frame(n=repeats,BlockSize=block.size,SD,Z,P.val)
		b<-list(SequenceLen,DStats,Jackknife)
		names(b)<-c("SequenceLen","DStats","Jackknife")
	}
	
	# Output a list.
	return(b)
}

# a: sequence alignment;
#	BStringSet object of at least 5 elements;
#	elements are arranged according to P1,P2,P3.1,P3.2,O;
#	introgression is from P3.1 to P2;
#	see Durand et al. (2011) for the meanings of P1,P2,P3.1,P3.2,O.
# genotype: logic;
#	if FALSE (default), introgression rate calculation is sequence-based;
#	if TRUE, introgression rate calculation is genotype-based.
# ns.seq: 5 integers; number of sequences from P1,P2,P3.1,P3.2,O when genotype=TRUE.
introgressed<-function(a,genotype=FALSE,ns.seq=c(1,1,1,1,1),missing.char=".")
{
	a<-as.matrix(a)
	fun<-DSeq
	
	# Filter for missing and non-biallelic (remove).
	Original<-ncol(a)
	index<-colSums(a==missing.char)==0
	a<-a[,index]
	index<-sapply(1:ncol(a),function(x)length(unique(a[,x]))<=2)
	a<-a[,index]
	AfterFilter<-ncol(a)
	SequenceLen<-data.frame(Original,AfterFilter)
	
	if(genotype)
	{
		# Calculate allele frequency.
		reference<-a[1,]
		freq<-numeric(0)
		j<-0
		for(i in ns.seq)
		{
			temp<-t(matrix(a[1:i+j,],i,AfterFilter))
			temp<-rowSums(temp!=reference)/ncol(temp)
			freq<-rbind(freq,temp)
			j<-j+i
		}
		a<-freq
		fun<-DGeno
	}
	
	# Calculate introgression rate.
	temp<-a
	
	a<-temp[-4,]
	DStats<-fun(a)
	abba<-DStats$ABBA
	baba<-DStats$BABA
	s1<-abba-baba
	b<-data.frame(ABBA1=abba,BABA1=baba)
	
	a<-temp[c(1,4,3,5),]
	DStats<-fun(a)
	abba<-DStats$ABBA
	baba<-DStats$BABA
	s2<-abba-baba
	b<-data.frame(b,ABBA2=abba,BABA2=baba)
	
	f<-s1/s2
	if(is.nan(f))f<-0
	b<-data.frame(f,b)
	b<-list(SequenceLen=SequenceLen,Introgressed=b)
	return(b)
}

#######=======================================================================#######
# My version of ABBA BABA rest.
# fSeq calculates sequence-based fHom from P3 to P2;
#	see Martin et al. (2014) for the meanding of fHom.
# fGeno calculates genotype-based fHom from P3 to P2.
# abba.baba2 runs my ABBA BABA test, knowing the gene flow is from P3 to P2.

# a: matrixed sequence alignment; individuals as rows.
fSeq<-function(a)
{
	abba<-(a[1,]==a[4,])&(a[2,]==a[3,])&(a[1,]!=a[2,])
	abba<-sum(abba)
	baba<-(a[1,]==a[3,])&(a[2,]==a[4,])&(a[1,]!=a[2,])
	baba<-sum(baba)
	aaba<-(a[1,]==a[2,])&(a[2,]==a[4,])&(a[1,]!=a[3,])
	aaba<-sum(aaba)
	f<-(abba-baba)/(aaba+abba)
	if(is.nan(f))f<-0
	data.frame(f,ABBA=abba,BABA=baba,AABA=aaba)
}

# a: allele frequency matrix; populations as rows.
fGeno<-function(a)
{
	abba<-(1-a[1,])*a[2,]*a[3,]*(1-a[4,])+a[1,]*(1-a[2,])*(1-a[3,])*a[4,]
	abba<-sum(abba)
	baba<-a[1,]*(1-a[2,])*a[3,]*(1-a[4,])+(1-a[1,])*a[2,]*(1-a[3,])*a[4,]
	baba<-sum(baba)
	aaba<-(1-a[1,])*(1-a[2,])*a[3,]*(1-a[4,])+a[1,]*a[2,]*(1-a[3,])*a[4,]
	aaba<-sum(aaba)
	f<-(abba-baba)/(aaba+abba)
	if(is.nan(f))f<-0
	data.frame(f,ABBA=abba,BABA=baba,AABA=aaba)
}

# a: sequence alignment;
#	BStringSet object of at least 4 elements;
#	elements are arranged according to P1, P2, P3, and O;
#	see Durand et al. (2011) for the meanings of P1, P2, P3, and O.
# genotype: logic;
#	if FALSE (default), f statistics calculation is sequence-based;
#	if TRUE, f statistics calculation is genotype-based.
# ns.seq: 4 integers; number of sequences from P1, P2, P3, and O when genotype=TRUE.
# missing.char: character for missing values in a.
# test: possible values are "n" (no test), "a" (Poisson),
#	"b" (bootstrap), and "j" (jackknife).
abba.baba2<-function(a,genotype=FALSE,ns.seq=c(1,1,1,1),
	missing.char=".",test="n",repeats=1000,block.size=1000)
{
	a<-as.matrix(a)
	fun<-fSeq
	
	# Filter for missing and non-biallelic (remove).
	Original<-ncol(a)
	index<-colSums(a==missing.char)==0
	a<-a[,index]
	index<-sapply(1:ncol(a),function(x)length(unique(a[,x]))<=2)
	a<-a[,index]
	AfterFilter<-ncol(a)
	SequenceLen<-data.frame(Original,AfterFilter)
	
	if(genotype)
	{
		# Calculate allele frequency.
		reference<-a[1,]
		freq<-numeric(0)
		j<-0
		for(i in ns.seq)
		{
			temp<-t(matrix(a[1:i+j,],i,AfterFilter))
			temp<-rowSums(temp!=reference)/ncol(temp)
			freq<-rbind(freq,temp)
			j<-j+i
		}
		a<-freq
		fun<-fGeno
	}
	
	# Calculate f statistics.
	fStats<-fun(a)
	
	# Do not implement significance test.
	if(test=="n")
	{
		b<-list(SequenceLen,fStats)
		names(b)<-c("SequenceLen","fStats")
	}
	
	# Test significance using half analytical method (Poisson).
	if(test=="a")	####### To be improved!!!
	{
		e.baba<-fStats$BABA/(1-fStats$f)	# Inaccurate when BABA is small.
		e.aaba<-fStats$AABA/(1-fStats$f)	# Inaccurate when AABA is small.
		e.abba<-e.baba
		d.baba<-rpois(repeats,e.baba)
		d.aaba<-rpois(repeats,e.aaba)
		d.abba<-rpois(repeats,e.abba)
		d.f<-(d.abba-d.baba)/(d.aaba+d.abba)
		P.val<-sum(d.f>fStats$f)/repeats
		Analytic<-data.frame(n=repeats,P.val)
		b<-list(SequenceLen,fStats,Analytic)
		names(b)<-c("SequenceLen","fStats","Analytic")
	}
	
	# Test significance using bootstrap.
	if(test=="b")
	{
		bootstrap<-sapply(1:repeats,function(x)
			fun(a[,sample(1:AfterFilter,AfterFilter,replace=TRUE)])$f
			)
		SD<-sd(bootstrap)
		Z<-abs(fStats$f/SD)
		P.val<-2*(1-pnorm(Z))	#######!!!
		Bootstrap<-data.frame(n=repeats,SD,Z,P.val)
		b<-list(SequenceLen,fStats,Bootstrap)
		names(b)<-c("SequenceLen","fStats","Bootstrap")
	}
	
	# Test significance using jackknife.
	if(test=="j")
	{
		max.rep<-AfterFilter-block.size+1
		drop.pos<-as.integer(seq(from=1,to=max.rep,length.out=repeats))
		jackknife<-sapply(1:repeats,function(x)
			fun(a[,-c(drop.pos[x]:(drop.pos[x]+block.size-1))])$f
			)
		SD<-sd(jackknife)
		Z<-abs(fStats$f/SD)
		P.val<-2*(1-pnorm(Z))	#######!!!
		Jackknife<-data.frame(n=repeats,BlockSize=block.size,SD,Z,P.val)
		b<-list(SequenceLen,fStats,Jackknife)
		names(b)<-c("SequenceLen","fStats","Jackknife")
	}
	
	# Output a list.
	return(b)
}

#######=======================================================================#######
# Some useful wrapper functions.
# slidingWindow runs ABBA BABA tests in sliding windows.
# popBatch runs a batch of ABBA BABA tests for population data.

# a: BStringSet object of sequence alignment;
#	all sequences MUST have the same length;
#	the nucleotides in the sequences must be ordered according to position.
# fun: a function to be implemented on each window.
# position: integral vector; positions of the nucleotide sites;
#	when NA, the increment of each nucleotide is one and the start is one.
slidingWindow<-function(a,fun,position=NA,windowSize=2e5,stepSize=1e5)
{
	if(is.na(position))
	{
		L<-width(a)[1]
		position<-1:L
	} else
	{
		L<-max(position)
	}
	result<-vector("list",0)
	i<-1
	while(TRUE)	# Sliding window
	{
		ST<-(i-1)*stepSize+1
		EN<-ST-1+windowSize
		index1<-which(position>=ST)[1]
		index2<-rev(which(position<=EN))[1]
		if(index2>=index1)
		{
			aSub<-subseq(a,index1,index2)
			result<-c(result,list(fun(aSub)))
		} else
		{
			result<-c(result,list(NA))
		}
		if(EN>=L)	# The window reaches the end of the chromosome
		{
			break
		}
		i<-i+1
	}
	return(result)
}

# a: sequence alignment, BStringSet object;
# popmap: two column data frame;
#	1st column sequence names;
#	2nd column population ID's.
# abPop: 4 element vector, the population ID's of P1, P2, P3, O.
popBatch<-function(a,popmap,abPop,n.random=10)
{
	pops<-as.character(sort(unique(popmap[,2])))
	indexList<-lapply(pops,function(x)which(popmap[,2]==x))
	names(indexList)<-pops
	abPop<-as.character(abPop)
	result<-integer(0)
	
	# Sequence-based.
	for(i in 1:n.random)
	{
		index<-sapply(abPop,function(x)sample(indexList[[x]],1))
		temp<-abba.baba2(a[index],missing.char=".",test="b")	##
		temp<-as.data.frame(c(P=popmap[index,1],temp))
		result<-rbind(result,temp)
	}
	
	# Genotype-based.
	index<-unlist(indexList[abPop])
	ns.seq<-sapply(abPop,function(x)length(indexList[[x]]))
	temp<-abba.baba2(a[index],genotype=TRUE,ns.seq,missing.char=".",test="b")	##
	result<-list(SequenceBased=result,GenotypeBased=temp)
	
	return(result)
}
