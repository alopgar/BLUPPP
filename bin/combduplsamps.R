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
	} else if(ncol(ctbl)>2){
		message("> WARNING!!! Sample ", n, " has more than 2 repetitions. Check your file.\n\n")
		next
	}
	colnames(ctbl) <- paste0("rep", 1:ncol(ctbl))
	names(naperc) <- colnames(ctbl)
	ctbl$comb <- ifelse(is.na(ctbl$rep1), ctbl$rep2,
						  ifelse(is.na(ctbl$rep2), ctbl$rep1,
								   ifelse(ctbl$rep1 != ctbl$rep2, ctbl[,names(which.min(naperc))], ctbl$rep1)))
	combstats <- data.frame("sample" = n,
	  "corresp_total" = round((sum(ctbl$rep1==ctbl$rep2, na.rm=T) + sum(is.na(ctbl$rep1) & is.na(ctbl$rep2)))/nrow(ctbl), 4),
	  "neqG" = sum(ctbl$rep1 == ctbl$rep2, na.rm = T),
	  "ndiffG" = sum(ctbl$rep1 != ctbl$rep2, na.rm = T),
	  "corresp_ident" = sum(ctbl$rep1 == ctbl$rep2, na.rm = T)/(sum(ctbl$rep1 == ctbl$rep2, na.rm = T)+sum(ctbl$rep1 != ctbl$rep2, na.rm = T)),
	  "nfilled" = sum((is.na(ctbl$rep1) & !is.na(ctbl$rep2)) | (!is.na(ctbl$rep1) & is.na(ctbl$rep2))),
	  "finalNA" = sum(is.na(ctbl$rep1) & is.na(ctbl$rep2)),
	  "finalNAperc" = round((sum(is.na(ctbl$comb))/nrow(ctbl))*100, 4))
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
