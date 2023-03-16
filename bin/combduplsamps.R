## Initial setup:
options(warn=-1)
packs <- c("dplyr", "tidyr", "stringr")
invisible(suppressMessages(lapply(packs, require, character.only = TRUE)))

args <- commandArgs(T)
duptbl <- args[1]
headers <- args[2]
out <- args[3]
storestats <- TRUE

## Load table & set duplicated names as "name.1":
if(headers == 0){h=F} else {h=T}
gd <- read.table(duptbl, sep = " ", stringsAsFactors = F, header = h)
gdids <- gd[,1]
gdids <- make.unique(gdids)

## Transpose table & replace -1 by NA:
gdt <- setNames(data.table::transpose(gd[,-1]), gdids)
row.names(gdt) <- gsub("\\.", "-", colnames(gd)[-1])
gdt[gdt == -1] <- NA

## Create combined table:
getmode <- function(vec){
	veccl <- vec[!is.na(vec)]
	uniqv <- unique(veccl)
	freqs <- tabulate(match(veccl, uniqv))
	if(length(vec[which(!is.na(vec))])==0){return(NA)}
	else if(length(freqs)==1){return(uniqv)}
	else if(length(freqs)>1 & var(freqs)==0){return(veccl[1])}
	else{return(uniqv[which.max(freqs)])}
}
all_cols_eq <- function(dat, getalleq=T, na.rm=F){
	if(na.rm){dat <- dat[complete.cases(dat),]}
	if(getalleq){aeqdat <- dat[apply(dat, 1, function(x) length(unique(x))) == 1,]}
	else{aeqdat <- dat[apply(dat, 1, function(x) length(unique(x))) > 1,]}
	return(aeqdat)
}

gdcomb <- data.frame("snps"=row.names(gdt))
badcorresp <- NULL
if(storestats){cstattbl <- data.frame()}
for(n in unique(gsub("\\.1", "", colnames(gdt)))){
	message("COMBINING SAMPLES WITH ID: ", n)
	ctbl <- select(gdt, starts_with(n))
	naperc <- apply(ctbl, 2, function(x){sum(is.na(x))/nrow(ctbl)})
	message(" NA proportion at each repetition:")
	print(naperc)
	message(" Most informative sample: ", names(which.min(naperc)))
	if(ncol(ctbl)==1){
		message("> WARNING!!! Sample ", n, " has no repetitions.\n\n")
		next
	} else if(ncol(ctbl) == 2){
		#Sort NA vector and table columns by increasing NA%
		naperc <- sort(naperc, decreasing = F)
		ctbl <- ctbl[,names(naperc)]
		colnames(ctbl) <- paste0("rep", 1:ncol(ctbl))
		names(naperc) <- colnames(ctbl)
		#Combine info: Get the most frequent non-NA value; if all non-NA are equally frequent, get the first one (rep1 is always the most informative).
		ctbl$comb <- apply(ctbl, 1, function(x){getmode(x)})
		#Combination statistics:
		combstats <- data.frame("sample" = n,
			"nrep" = ncol(ctbl)-1,
			"corresp_total" = round((sum(ctbl$rep1==ctbl$rep2, na.rm=T) + sum(is.na(ctbl$rep1) & is.na(ctbl$rep2)))/nrow(ctbl), 4),
			"neqG" = sum(ctbl$rep1 == ctbl$rep2, na.rm = T),
			"ndiffG" = sum(ctbl$rep1 != ctbl$rep2, na.rm = T),
			"corresp_ident" = round(sum(ctbl$rep1 == ctbl$rep2, na.rm = T)/(sum(ctbl$rep1 == ctbl$rep2, na.rm = T)+sum(ctbl$rep1 != ctbl$rep2, na.rm = T)), 4),
			"nfilled" = sum((is.na(ctbl$rep1) & !is.na(ctbl$rep2)) | (!is.na(ctbl$rep1) & is.na(ctbl$rep2))),
			"finalNA" = sum(is.na(ctbl$rep1) & is.na(ctbl$rep2)),
			"finalNAperc" = round((sum(is.na(ctbl$comb))/nrow(ctbl))*100, 4))
	} else if(ncol(ctbl) >= 3){
		#Sort NA vector and table columns by increasing NA%
		naperc <- sort(naperc, decreasing = F)
		ctbl <- ctbl[,names(naperc)]
		colnames(ctbl) <- paste0("rep", 1:ncol(ctbl))
		names(naperc) <- colnames(ctbl)
		#Combine info:
		ctbl$comb <- apply(ctbl, 1, function(x){getmode(x)})		
		message("> WARNING!!! Sample ", n, " has more than 2 repetitions. Combination statistics are done comparing all repetitions at a time.\n\n")
		#Combination statistics:
		reptbl <- select(ctbl, starts_with("rep"))
		combstats <- data.frame("sample" = n,
			"nrep" = ncol(reptbl),
			"corresp_total" = round(nrow(all_cols_eq(reptbl,na.rm=F))/nrow(ctbl), 4),
			"neqG" = nrow(all_cols_eq(reptbl,na.rm=T)),
			"ndiffG" = nrow(all_cols_eq(reptbl,getalleq=F,na.rm=T)),
			"corresp_ident" = round(nrow(all_cols_eq(reptbl,na.rm=T))/(nrow(all_cols_eq(reptbl,na.rm=T)) + nrow(all_cols_eq(reptbl,getalleq=F,na.rm=T))), 4),
			"nfilled" = nrow(reptbl[rowSums(!is.na(reptbl)) > 0 & rowSums(!is.na(reptbl)) < ncol(reptbl),]), 
			"finalNA" = nrow(reptbl[rowSums(!is.na(reptbl)) == 0,]), 
			"finalNAperc" = round((sum(is.na(ctbl$comb))/nrow(ctbl))*100, 4))
		rm(reptbl)
	}
	if(storestats){
		cstattbl <- rbind(cstattbl, combstats)
		message("\n")
	} else {
		message(" Total correspondence:                ", combstats$corresp_total, "\n",
				" Number of equal genotypes:           ", combstats$neqG, "\n",
				" Number of different genotypes:       ", combstats$ndiffG, "\n",
				" No-NA correspondence:                ", combstats$corresp_ident, "\n",
				" Number of genotypes filled from NA:  ", combstats$nfilled, "\n",
				" Final number of NAs:                 ", combstats$finalNA, "\n",
				" Final percentage of NAs:             ", combstats$finalNAperc, "\n")
	}
	if(combstats$corresp_ident >= 0.95){
		gdcomb <- cbind(gdcomb,ctbl$comb)
		names(gdcomb)[ncol(gdcomb)] <- n
	} else {
		badcorresp <- c(badcorresp, n)
	}
}

## Transpose, format & save combined table:
gdcomb_fin <- setNames(data.table::transpose(gdcomb[,-1]), gdcomb$snps)
gdcomb_fin <- cbind("fish-id"=colnames(gdcomb)[-1],gdcomb_fin)

write.table(gdcomb_fin, out, quote = F, sep=" ", row.names = F)
write.table(badcorresp, "lowcorrsamples.txt", quote = F, row.names = F, , col.names = F)
if(storestats){write.table(cstattbl, "duplstats.txt", quote = F, sep=" ", row.names = F)}
