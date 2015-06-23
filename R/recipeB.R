`rbreak` <- function(filename, re=0.5, be=0.9, fe=0.5, studies = c(1,2), samplesizes=c(845,40), sensitivities=c(0.8, 0.6)) {
	library.dynam("recipeB", "recipeB", paste(system.file(package="recipeB"), "/../", sep=""))
	outfile = paste(substr(filename, 1, nchar(filename)-4), ".out", sep="")
	logfile = paste(substr(filename, 1, nchar(filename)-4), ".log", sep="")
	d = read.table(filename); types = unique(d[,6]); types=types[order(types)]; stu = unique(d[,7]);
	if( length(stu)!= length(studies) ) { warnings("You are a silly person, more studies are in the input file than specified to this function!! - my results will now be rubbish!") }
	system(paste("rm", outfile, logfile, sep=" ")) 
	.C("rbreak",
		"f1" = as.character(filename) 
		,"f2" = as.character(outfile)
		,"f3" = as.character(logfile)
		,"re" = as.double(re) 
		,"be" = as.double(be) 
		,"fe" = as.double(fe)
		,"types" = as.double(types) 
		,"studies" = as.double(studies)
		,"samplesizes" = as.double(samplesizes) 
		,"sensitivities" = as.double(sensitivities) 
		,"size1" = as.integer(length(types)) 
		,"size2" = as.integer(length(studies))
		,"PACKAGE" = "recipeB")
	library.dynam.unload("recipeB", system.file(package="recipeB"))	
}

`uber.overlap` <- function(set1=NULL, set2=NULL, re1=0, re2=0) {
	counts = coverlap(set1, set2, re1, re2)
	ind = length(counts[1,]); set1 = counts
	pin = 1; forward_overs = vector(); backwards_overs = vector(); overlap_indexes = vector()
	
	u = unique(set1[,1]); idss = 1:length(set2[,1])
	
	for(x in 1:length(u)) {	
		s1 = set1[as.character(set1[,1])==as.character(u[x]),]
		s2 = set2[as.character(set2[,1])==as.character(u[x]),]
		ids = idss[as.character(set2[,1])==as.character(u[x])]
		ss1 = s1[,2]; ss2 = s1[,3]; ss3 = s2[,2]; ss4 = s2[,3];
		for(y in 1:length(ss1)) {
			r1 = vector(length=set1[pin,ind])
			r2 = vector(length=set1[pin,ind])
			indexes = vector(length=set1[pin,ind])
			if(length(r1)==0) { r1 = 0; r2 = 0; indexes=0;}	
			r = .C( "uoverlap", 
				   "rep1" = as.double(re1),
				   "rep2" = as.double(re2),
				   "a" = as.integer(ss1[y]), 
				   "b" = as.integer(ss2[y]), 
				   "c" = as.integer(ss3), 
				   "d" = as.integer(ss4),
				   "ss" = as.integer(ids),
				   "r1" = as.double(r1),
				   "r2" = as.double(r2),
				   "indexes" = as.integer(indexes),
				   "size1" = as.integer(length(ss3)),
				   PACKAGE = "recipeB" )
			forward_overs[pin]=paste(r$r1, collapse=":")
			backwards_overs[pin]=paste(r$r2, collapse=":")
			overlap_indexes[pin]=paste(r$indexes, collapse=":")
			pin = pin+1
		}
	}
	dat = cbind(set1, forward_overs, backwards_overs, overlap_indexes)
invisible(dat)
}


`ooverlap` <- function(set1, set2) {
	r = NULL;
	si1 = length(set1[,1]); si2 = length(set2[,1]); 
	set1=set1[set1[,2]<=set1[,3],]; set2=set2[set2[,2]<=set2[,3],];
	if( si1>length(set1[,1]) | si2>length(set2[,1]) ) { warning("i removed some rubbish rows - where start was larger then stop!") }
	if(length(set1)>0 & length(set2)>0) {
	set1 = set1[order(set1[,1], set1[,2],set1[,3]),];
	set2 = set2[order(set2[,1], set2[,2],set2[,3]),];
	u = unique(set1[,1])
	for(x in 1:length(u)) {
		s1 = set1[as.character(set1[,1])==u[x],]
		s2 = set2[as.character(set2[,1])==u[x],]
		r1 = vector(length=length(s1[,1]));
		r2 = vector(length=length(s1[,1]));
		if(length(s1)>0 & length(s2)>0) {
		if(is.null(dim(s1))) { s1 = rbind(s1,s1); }
		if(is.null(dim(s2))) { s2 = rbind(s2,s2); }
			res <- .C("overlap",
				"a" = as.integer(s1[,2]) 
				,"b" = as.integer(s1[,3])
				,"c" = as.integer(s2[,2]) 
				,"d" = as.integer(s2[,3]) 
				,"r1" = as.double(r1) 
				,"r2" = as.double(r2) 
				,"size1" = as.integer(length(s1[,1]))
				,"size2" = as.integer(length(s2[,2])) 
				,"PACKAGE" = "recipeB")
			r = rbind(r, cbind(res$r1, res$r2))
		} else {
			r = rbind(r, cbind(r1, r2))
		}
	}
	rset = cbind(set1,r)
	return(rset)
	} else {
		 warning("everything was rubbish rows - where start was larger then stop!")
	}
}


`toverlap` <- function(set1, set2=NULL, rarecommon=0.01, sep=T) {
	if(is.null(set2)) { data(cnv2.1_separate_plusX); set2=cnv2.1_separate_plusX; } 
	classover <- function(set1, set2) {
		r = NULL;
		si1 = length(set1[,1]); si2 = length(set2[,1]); 
		set1=set1[set1[,2]<=set1[,3],]; set2=set2[set2[,2]<=set2[,3],];
		if( si1>length(set1[,1]) | si2>length(set2[,1]) ) { warning("i removed some rubbish rows - where start was larger then stop!") }
		if(length(set1)>0 & length(set2)>0) {
			set1 = set1[order(set1[,1], set1[,2],set1[,3]),];
			set2 = set2[order(set2[,1], set2[,2],set2[,3]),];
			u = unique(set1[,1])
			for(x in 1:length(u)) {
				s1 = set1[as.character(set1[,1])==u[x],]
				s2 = set2[as.character(set2[,1])==u[x],]
				r1 = vector(length=length(s1[,1]));
				r2 = vector(length=length(s1[,1]));
				if(length(s1)>0 & length(s2)>0) {
					if(is.null(dim(s1))) { s1 = rbind(s1,s1); }
					if(is.null(dim(s2))) { s2 = rbind(s2,s2); }
					res <- .C("toverlap",
							  "a" = as.integer(s1[,2]) 
							  ,"b" = as.integer(s1[,3])
							  ,"c" = as.integer(s2[,2]) 
							  ,"d" = as.integer(s2[,3]) 
							  ,"r1" = as.double(r1) 
							  ,"r2" = as.double(r2) 
							  ,"size1" = as.integer(length(s1[,1]))
						  	,"size2" = as.integer(length(s2[,2])) 
						  	,"PACKAGE" = "recipeB")
					r = rbind(r, cbind(res$r1, res$r2))
				} else {
					r = rbind(r, cbind(r1, r2))
				}
			}
			rset = cbind(set1,r)
		return(rset)
		} else {
			warning("everything was rubbish rows - where start was larger then stop!")
		}
	}
	
	rare_del = set2[set2[,11]<rarecommon & set2[,13]<=0,]
	common_del = set2[set2[,11]>=rarecommon & set2[,13]<=0,]
	rare_dup = set2[set2[,11]<rarecommon & set2[,13]>=0,]
	common_dup = set2[set2[,11]>=rarecommon & set2[,13]>=0,]
	
	dels = NULL
	if(length(set1[set1[,4]<0,1])>0 ) {
		rdels = classover(set1[set1[,4]<0,], rare_del)
		cdels = classover(set1[set1[,4]<0,], common_del)
		len = length(cdels[1,])-1
		dels=cbind(rdels, cdels[,len:length(cdels[1,])])
	}
	
	dups=NULL
	if(length(set1[set1[,4]>0,1])>0 ) {
		rdups = classover(set1[set1[,4]>0,], rare_dup)
		cdups = classover(set1[set1[,4]>0,], common_dup)
		len = length(cdups[1,])-1
		dups=cbind(rdups, cdups[,len:length(cdups[1,])])
	}
	
	dat = rbind(dels,dups)
	dat = dat[order(dat[,1], dat[,2], dat[,3]),]
	
	return(dat)
}

`roverlap` <- function(set1, set2) {
	r = NULL;
	si1 = length(set1[,1]); si2 = length(set2[,1]); 
	set1=set1[set1[,2]<=set1[,3],]; set2=set2[set2[,2]<=set2[,3],];
	if( si1>length(set1[,1]) | si2>length(set2[,1]) ) { warning("i removed some rubbish rows - where start was larger then stop!") }
	if(length(set1)>0 & length(set2)>0) {
	set1 = set1[order(set1[,1], set1[,2],set1[,3]),];
	set2 = set2[order(set2[,1], set2[,2],set2[,3]),];
	u = unique(set1[,1])
	for(x in 1:length(u)) {
		s1 = set1[as.character(set1[,1])==u[x],]
		s2 = set2[as.character(set2[,1])==u[x],]
		r1 = vector(length=length(s1[,1]));
		r2 = vector(length=length(s1[,1]));
		if(length(s1)>0 & length(s2)>0) {
		if(is.null(dim(s1))) { s1 = rbind(s1,s1); }
		if(is.null(dim(s2))) { s2 = rbind(s2,s2); }
			res <- .C("roverlap",
				"a" = as.integer(s1[,2]) 
				,"b" = as.integer(s1[,3])
				,"c" = as.integer(s2[,2]) 
				,"d" = as.integer(s2[,3]) 
				,"r1" = as.double(r1) 
				,"r2" = as.double(r2) 
				,"size1" = as.integer(length(s1[,1]))
				,"size2" = as.integer(length(s2[,2])) 
				,"PACKAGE" = "recipeB")
			r = rbind(r, cbind(res$r1, res$r2))
		} else {
			r = rbind(r, cbind(r1, r2))
		}
	}
	rset = cbind(set1,r)
	return(rset)
	} else {
		 warning("everything was rubbish rows - where start was larger then stop!")
	}
}

`foverlap` <- function(set1, set2=NULL) {
	r = NULL;
	if(is.null(set2)) { data(cnv2.1); set2=cnv2.1; } 
	si1 = length(set1[,1]); si2 = length(set2[,1]); 
	set1=set1[set1[,2]<=set1[,3],]; set2=set2[set2[,2]<=set2[,3],];
	if( si1>length(set1[,1]) | si2>length(set2[,1]) ) { warning("i removed some rubbish rows - where start was larger then stop!") }
	if(length(set1)>0 & length(set2)>0) {
	set1 = set1[order(set1[,1], set1[,2],set1[,3]),];
	set2 = set2[order(set2[,1], set2[,2],set2[,3]),];
	u = unique(set1[,1])
	for(x in 1:length(u)) {
		s1 = set1[as.character(set1[,1])==u[x],]
		s2 = set2[as.character(set2[,1])==u[x],]
		r1 = vector(length=length(s1[,1]));
		r2 = vector(length=length(s1[,1]));
		rf1 = vector(length=length(s1[,1]));
		rf2 = vector(length=length(s1[,1]));
		rf3 = vector(length=length(s1[,1]));
		rty = vector(length=length(s1[,1]));
		if(length(s1)>0 & length(s2)>0) {
		if(is.null(dim(s1))) { s1 = rbind(s1,s1); }
		if(is.null(dim(s2))) { s2 = rbind(s2,s2); }
			res <- .C("foverlap",
				"a" = as.integer(s1[,2]) 
				,"b" = as.integer(s1[,3])
				,"c" = as.integer(s2[,2]) 
				,"d" = as.integer(s2[,3])
				,"fre1" = as.double(s2[,5]) 
				,"fre2" = as.double(s2[,8]) 
				,"fre3" = as.double(s2[,11]) 
				,"r1" = as.double(r1) 
				,"r2" = as.double(r2)
				,"rf1" = as.double(rf1) 
				,"rf2" = as.double(rf2) 
				,"rf3" = as.double(rf3)
				,"ty" = as.double(s2[,13]) 
				,"rty" = as.double(rty)
				,"size1" = as.integer(length(s1[,1]))
				,"size2" = as.integer(length(s2[,2])) 
				,"PACKAGE" = "recipeB")
			r = rbind(r, cbind(res$r1, res$r2, res$rf1, res$rf2, res$rf3, res$rty))
		} else {
			r = rbind(r, cbind(r1, r2, -1, -1, -1, -1))
		}
	}
	rset = cbind(set1,r)
	return(rset)
	} else {
		 warning("everything was rubbish rows - where start was larger then stop!")
	}
}


`soverlap` <- function(set1, set2=NULL) {
	r = NULL;
	if(is.null(set2)) { data(cnv2.1_separate); set2=cnv2.1_separate; } 
	si1 = length(set1[,1]); si2 = length(set2[,1]); 
	set1=set1[set1[,2]<=set1[,3],]; set2=set2[set2[,2]<=set2[,3],];
	if( si1>length(set1[,1]) | si2>length(set2[,1]) ) { warning("i removed some rubbish rows - where start was larger then stop!") }
	
	if(length(set1)>0 & length(set2)>0) {
	set1 = set1[order(set1[,1], set1[,2],set1[,3]),];
	set2 = set2[order(set2[,1], set2[,2],set2[,3]),];
	u = unique(set1[,1])
	for(x in 1:length(u)) {
		
		s1 = set1[as.character(set1[,1])==u[x],]
		rr1 = 0; rr2 = 0; rrf1 = 0; rrf2 = 0; rrf3 = 0; rrty = 0;
		
		studies = unique(set2[,15]); first = 1; fin2 = NULL;
		for( y in 1:length(studies) ) {
		s2 = set2[as.character(set2[,1])==u[x] & set2[,15]==studies[y],]
		s2 = s2[order(s2[,1],s2[,2],s2[,3]),]
		
		r1 = vector(length=length(s1[,1]));
		r2 = vector(length=length(s1[,1]));
		rf1 = vector(length=length(s1[,1]));
		rf2 = vector(length=length(s1[,1]));
		rf3 = vector(length=length(s1[,1]));
		rty = vector(length=length(s1[,1]));
		if(length(s1)>0 & length(s2)>0) {
		if(is.null(dim(s1))) { s1 = rbind(s1,s1); }
		if(is.null(dim(s2))) { s2 = rbind(s2,s2); }
			
		res <- .C("foverlap",
				"a" = as.integer(s1[,2]) 
				,"b" = as.integer(s1[,3])
				,"c" = as.integer(s2[,2]) 
				,"d" = as.integer(s2[,3])
				,"fre1" = as.double(s2[,5]) 
				,"fre2" = as.double(s2[,8]) 
				,"fre3" = as.double(s2[,11]) 
				,"r1" = as.double(r1) 
				,"r2" = as.double(r2)
				,"rf1" = as.double(rf1) 
				,"rf2" = as.double(rf2) 
				,"rf3" = as.double(rf3)
				,"ty" = as.double(s2[,13]) 
				,"rty" = as.double(rty)
				,"size1" = as.integer(length(s1[,1]))
				,"size2" = as.integer(length(s2[,2])) 
				,"PACKAGE" = "recipeB")
			
				rr1 = res$r1; rr2 = res$r2; rrf1 = res$rf1; rrf2 = res$rf2; rrf3 = res$rf3; rrty = res$rty;
				fin1 = cbind(rr1, rr2, rrf1, rrf2, rrf3, rrty);
				
				if(first) {
					fin2 = fin1;
					first=0;
				} else {
					for(z in 1:length(fin1[1,])) {
						if(fin1[z,1] >= fin2[z,1]) {
							fin2[z,] = fin1[z,]
						}
						if(fin1[z,1]>0) {
							if(fin1[z,6] != fin2[z,6] | fin2[z,6]==0) {
								fin2[z,6]=0;
							}	
						}
					}
				}
				
			} else {
				r = rbind(r, cbind(r1, r2, -1, -1, -1, -1))
			}
	}
			r = rbind(r, fin2)
		
	}
	rset = cbind(set1,r)
	return(rset)
	} else {
		 warning("everything was rubbish rows - where start was larger then stop!")
	}
}

`coverlap` <- function(set1, set2=NULL, rep1=0.5, rep2=0) {
	r = NULL;
	set1 = set1[order(set1[,1], set1[,2],set1[,3]),];
	set2 = set2[order(set2[,1], set2[,2],set2[,3]),];
	u = unique(set1[,1])
	for(x in 1:length(u)) {
		
		s1 = set1[as.character(set1[,1])==u[x],]
		s2 = set2[as.character(set2[,1])==u[x],]
		s1 = s1[order(s1[,1],s1[,2],s1[,3]),]
		s2 = s2[order(s2[,1],s2[,2],s2[,3]),]	
		rr1 = 0; rr2 = 0; ccon =0;
		
		r1 = vector(length=length(s1[,1]));
		r2 = vector(length=length(s1[,1]));
		con = vector(length=length(s1[,1]));
			
		res <- .C("coverlap",
				"a" = as.integer(s1[,2]) 
				,"b" = as.integer(s1[,3])
				,"c" = as.integer(s2[,2]) 
				,"d" = as.integer(s2[,3])
				,"con" = as.double(con)
				,"r1" = as.double(r1)
				,"r2" = as.double(r2)
				,"rep1" = as.double(rep1)
				,"rep2" = as.double(rep2)
				,"size1" = as.integer(length(s1[,1]))
				,"size2" = as.integer(length(s2[,2])) 
				,"PACKAGE" = "recipeB")
			
				rr1 = res$r1; rr2 = res$r2; ccon=res$con;
				fin1 = cbind(rr1, rr2, ccon);
				

			r = rbind(r, cbind(s1, fin1))
	}

	return(r)
}

`soverlap2` <- function(set1, set2=NULL) {
	r = NULL;
	if(is.null(set2)) { data(cnv2.1_separate); set2=cnv2.1_separate; } 
	si1 = length(set1[,1]); si2 = length(set2[,1]); 
	set1=set1[set1[,2]<=set1[,3],]; set2=set2[set2[,2]<=set2[,3],];
	if( si1>length(set1[,1]) | si2>length(set2[,1]) ) { warning("i removed some rubbish rows - where start was larger then stop!") }
	if(length(set1)>0 & length(set2)>0) {
	set1 = set1[order(set1[,1], set1[,2],set1[,3]),];
	set2 = set2[order(set2[,1], set2[,2],set2[,3]),];
	u = unique(set1[,1])
	for(x in 1:length(u)) {
		s1 = set1[as.character(set1[,1])==u[x],]
		s2 = set2[as.character(set2[,1])==u[x],]
		r1 = vector(length=length(s1[,1]));
		r2 = vector(length=length(s1[,1]));
		rf1 = vector(length=length(s1[,1]));
		rf2 = vector(length=length(s1[,1]));
		rf3 = vector(length=length(s1[,1]));
		rty = vector(length=length(s1[,1]));
		if(length(s1)>0 & length(s2)>0) {
		if(is.null(dim(s1))) { s1 = rbind(s1,s1); }
		if(is.null(dim(s2))) { s2 = rbind(s2,s2); }
			res <- .C("foverlap",
				"a" = as.integer(s1[,2]) 
				,"b" = as.integer(s1[,3])
				,"c" = as.integer(s2[,2]) 
				,"d" = as.integer(s2[,3])
				,"fre1" = as.double(s2[,5]) 
				,"fre2" = as.double(s2[,8]) 
				,"fre3" = as.double(s2[,11]) 
				,"r1" = as.double(r1) 
				,"r2" = as.double(r2)
				,"rf1" = as.double(rf1) 
				,"rf2" = as.double(rf2) 
				,"rf3" = as.double(rf3)
				,"ty" = as.double(s2[,13]) 
				,"rty" = as.double(rty)
				,"size1" = as.integer(length(s1[,1]))
				,"size2" = as.integer(length(s2[,2])) 
				,"PACKAGE" = "recipeB")
			r = rbind(r, cbind(res$r1, res$r2, res$rf1, res$rf2, res$rf3, res$rty))
		} else {
			r = rbind(r, cbind(r1, r2, -1, -1, -1, -1))
		}
	}
	rset = cbind(set1,r)
	return(rset)
	} else {
		 warning("everything was rubbish rows - where start was larger then stop!")
	}
}


set.params <- function(plist = NULL) {
	if(plist==NULL) {
		warning("i now expect that all study field in your input file is set to zero")
		plist$studies = 0;
		plist$samplesizes=1;
		plist$sensitivities=1;
	}
		.C("setallparams"
			,"studies" = as.integer(plist$studies)
			,"samplesizes" = as.integer(plist$samplesizes)
			,"sensitivities" = as.integer(plist$sensitivities)
			,"size" = as.integer(length(plist$studies))
			,"PACKAGE" = "recipeB"
		)
	return(plist)	
}