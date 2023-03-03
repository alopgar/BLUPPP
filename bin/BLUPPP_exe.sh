#!/bin/bash

#############################################################################################################
##### 											BLUPPP													#####
#####  BLUP Pre-Processing is a program that prepares genotype files for their use in BLUPf90 software  #####
#############################################################################################################

usage(){
cat << _EUs_
$(basename "$0") [PARAMETER FILE] [OPTIONS]... -- Initiates the pre-processing of files for BLUPF90 pipeline.\n
DESCRIPTION:\n
	\tThis program prepares the files for BLUPF90 analysis, including the formatting of the genotypes file,\n
	\ta pipeline of data pruning and the preparation of the BLUPF90 directory.\n
	\tIt also includes an additional pipeline for BLUP using ASRgenomics R package.\n
	\tIt requires a parameter file, which is further described.\n
	\tIMPORTANT!! The installation path of BLUPPP must be checked inside this script.\n\n
PARAMETER FILE:\n
	\t-i, --input\t\tMANDATORY. Path and name of the parameter file (.par). It must contain the following variables.\n
				\t\t\t\t  > WDIR: Working directory; where the outputs will be stored.\n
				\t\t\t\t  > genotypes: Genotype file which will be used and pre-processed to follow BLUPF90 requirements (if needed).\n
				\t\t\t\t  > phenot: File including phenotypic values of each individual, in the format required by BLUPF90.\n
				\t\t\t\t  > parent: File including parental relationships of each individual, in the format required by BLUPF90.\n
				\t\t\t\t  > BF90D: Directory containing data and results from BLUPF90 analysis. Will be created if missing.\n
				\t\t\t\t  > bf90param: Parameter file for BLUPF90 analysis.\n
				\t\t\t\t  > Other processing parameters: map, appendg, headers, NAto5, modif, EXTRA_PROC() (see details in example file).\n\n
OPTIONS:\n
	\t-p, --processing\tValues: blupf90, asrgenomics, custom. If specified, the genotypes file will be pre-processed and converted\n
			\t\t\t\t   to a format compatible either with BLUPF90 or with ASRGenomics.\n
	\t-l, --plink\t\tIf flag is specified, PLINK filters from plink_filt.sh will be executed. A .map file must be included \n
			\t\t\t\t  in the parameter file for this step.\n
	\t-R, --Rmatrix\t\tIf flag is specified, the H matrix will be created through ASRgenomics in R.\n
	\t-h, --help\t\tPrints this screen.\n
_EUs_
}

JOINORSEP(){
	modnamenm="${1%.*}"
	newname=$modnamenm'_'$2
	# Length params calculation:
	if [ "$2" = 'join' ]; then
		#Extra modification: When the ids have variable lengths, the id column must be trailed with spaces to reach the same length:
		ids_len=$(awk 'BEGIN{FS=" ";OFS = "\t"} {print $1}' $1 | awk '{print length}' | sort -nrk1,1 | head -1)
		awk -F " " -v num=$ids_len '{printf "%-*.*s\n", num, num, $1}' $1 > ids.tmp
		#Number of SNPs: n=(nchar/2)+1
		cut -f2- -d ' ' $1 > snps.tmp
		snps_len_sp=$( awk '{print length}' snps.tmp | sort -nrk1,1 | head -1)
		snps_len=$(( ( $snps_len_sp / 2 ) + 1 ))
		#joinorsep.f90 input file:
		paste -d ' ' ids.tmp snps.tmp > $modnamenm$3
		rm ids.tmp snps.tmp
	elif [ "$2" = 'separate' ]; then
		#Number of SNPs: n=nchar (in column 2)
		snps_len=$(awk 'BEGIN{FS=" ";OFS = "\t"} {print $2}' $1 | awk '{print length}' | sort -nrk1,1 | head -1)
		cp $1 $modnamenm$3
	fi
	## Fortran param file:
	printf "%s\n" "INPUT_FILE: '$WDIR/$modnamenm$3'" \
		"OUTPUT_FILE: '$WDIR/$newname$3'" \
		"OPERATION: '$2'" \
		"LENGTH_IDS: $ids_len" \
		"NUMBER_SNPS: $snps_len" > $WDIR/filter_param.txt
	gfortran $BINPATH/joinorsep.f90 -o $BINPATH/out_filter.x
	$BINPATH/out_filter.x $WDIR/filter_param.txt
	rm $modnamenm$3 $WDIR/filter_param.txt
}

## INITIAL VARIABLES
#############################################################################################################
## Paths & parameter file:
export BINPATH=/home/adrian.lopez@cai.cookeaqua.com/bin/BLUPPP
parfile=
## Options:
processing=
plink=false
Rmatrix=false

## RUNNING OPTIONS
#############################################################################################################
OPTS=`getopt -o i:p:l::R::h --long input:,processing:,plink::,Rmatrix::,help -- "$@"`
eval set -- "$OPTS"

while true; do
	case $1 in
		-i | --input)
			export parfile=$2; shift 2 ;;	
		-p | --processing)
			export processing=$2; shift 2 ;;
		-l | --plink)
			export plink=true; shift 2 ;;
		-R | --Rmatrix)
			export Rmatrix=true; shift 2 ;;	
		-h | --help)
			echo -e $(usage) | less ; exit ;;
		--) shift ; break ;;
        *) echo "Script definition error! Seems that one or more parameters are not valid..." ; exit 1 ;;
	esac
done

### START RUNNING
#############################################################################################################
printf "\n%s" "
##################################################
#####                 BLUPPP                 #####
#####  Running BLUP Pre-Processing pipeline  #####
##################################################
"

# PARAMETER FILE
## Read parameter file & convert to unix if necessary:
if [ -z "$parfile" ]; then
	printf "\n%s\n" "ERROR: Parameter file not specified."
	exit 1
elif [ ! -f "$parfile" ]; then
	printf "\n%s\n" "ERROR: Parameter file does not exist."
	exit 1
else
	printf "\n%s\n" "PARAMETER FILE: $parfile"
fi
if awk '/\r$/{exit 0}; 1{exit 1}' $parfile; then
	dos2unix $parfile
fi
source $parfile

## Assign variables for file names and extensions:
filename_ph=$(basename -- "$phenot")
filename_pa=$(basename -- "$parent")
filename_ge=$(basename -- "$genotypes")
extension_i=$([[ "$filename_ge" = *.* ]] && echo ".${filename_ge##*.}" || echo '')
extension_o='.txt'
export gname="${filename_ge%.*}"

printf "%s\n" "
WORKING DIRECTORY: $WDIR
GENOTYPES: $filename_ge
EXTENSION: $gname EXTENSION > $extension_i WILL BE CONVERTED TO > $extension_o
ANCESTRY FILE: $filename_pa
PHENOTYPES: $filename_ph

PROCESSING?? $processing
GENOTYPES TO APPEND: $appendg
PLINK MAP FILE: $smap
RMATRIX?? $Rmatrix"

cd $WDIR

# BLUPF90 file structure:
if [ ! -z "$BF90D" ]; then
	if [ "$Rmatrix" = false ]; then
		export BDIR=$WDIR/$BF90D
		# Initial controls:
		if [ -d $BDIR ]; then
			printf "\nWARNING: $BDIR already exists in the specified path and won't be replaced."
		elif [ ! -d $BDIR ]; then
			printf "\n> Creating BLUPF90 directory: $BDIR"
			mkdir -p $BDIR
		fi
		
		# Check file location:
		if [ ! -f $BDIR/$filename_ph ]; then
			cp $phenot $BDIR/$filename_ph
		fi
		if [ ! -f $BDIR/$filename_pa ]; then
			cp $parent $BDIR/$filename_pa
		fi
	fi
else
	printf "\n> BLUPF90 directory not specified. Files won't be placed for BLUPF90 preparation."
fi

### PROCESSING OF GENOTYPES
#############################################################################################################
if [ ! -z "$processing" ]; then
	# Activate processing if required:
	printf "\n\n%s\n" "--------- PROCESSING OF GENOTYPES REQUESTED ---------"
	if awk '/\r$/{exit 0}; 1{exit 1}' $genotypes; then
		dos2unix $genotypes
	fi
	
	## a) Add new genotypes (we need to keep in mind that the SNPs might not be the same in multiple datasets):
	if [ ! -z "$appendg" ]; then
		if [ "$headers" = 1 ] || [ "$headers" = 2 ]; then
			gname=$gname"_add"
			printf "  > Adding new genotypes from $appendg. NEW FILE: $gname.tmp \n"
			echo "options(warn=-1)
			g1 <- read.table('$genotypes', sep = ' ', stringsAsFactors = F, header = T)
			g2 <- read.table('$appendg', sep = ' ', stringsAsFactors = F, header = T)
			gu <- plyr::rbind.fill(g1,g2)
			snpcols <- colnames(Filter(is.numeric, gu))
			idcols <- colnames(gu)[!(colnames(gu) %in% snpcols)]
			gu <- gu[,c(idcols,sort(snpcols))]
			write.table(gu, '$gname.tmp', quote = F, sep=' ', row.names = F)" > R_add.R
			Rscript R_add.R
			rm R_add.R
		elif [ "$headers" = 'n' ]; then
			printf "\nERROR: It is not possible to merge genotype files if SNP names (i.e., headers) are not present."
			exit 1
		fi
	elif [ -z "$appendg" ]; then
		printf "  > Append new genotypes not requested; just copying the original file. NEW FILE: $gname.tmp \n"
		cp $genotypes $gname.tmp
	fi
	
	## b) Additional processing to get the required format:
	if [ ! -z $(declare -F "EXTRA_PROC") ]; then
		printf "  > Custom additional processing requested. Applying EXTRA_PROC function. UPDATING FILE: $gname.tmp \n"
		EXTRA_PROC "$gname.tmp"
	else
		printf "  > Custom additional processing not requested (EXTRA_PROC function not defined) \n"
	fi
	
	## c) Combine duplicated rows (genotyped individuals):
	if [ "$duplicates" = 'y' ]; then
		printf "  > Combining info from duplicated individuals. UPDATING FILE: $gname.tmp \n"
		awk 'seen[$1]++ {print $1}' $gname.tmp > duplicatedsamps.tmp
		awk -F" " 'NR==FNR{a[$1]; next} FNR==1 || $1 in a' duplicatedsamps.tmp $gname.tmp > $gname.dupl.tmp
		### c.1) R script for combining duplicates:
		Rscript $BINPATH/combduplsamps.R $gname.dupl.tmp $headers $gname.dupl.comb.tmp > combduplsamps.out
		#Delete original entries of combined duplicates (i.e., >95% correspondence) & add them after combination:
		awk 'NR>1 && NR==FNR {a[$1]} FNR!=NR {printf("%s",!($1 in a)? $0"\n": "")}' $gname.dupl.comb.tmp $gname.tmp > rmgooddups.tmp
		tail -n +2 $gname.dupl.comb.tmp | cat rmgooddups.tmp - > combdups.tmp
		#Remove duplicated samples which were not combined (i.e., <95% correspondence):
		awk 'NR==FNR {a[$0]} FNR!=NR {printf("%s",!($1 in a)? $0"\n": "")}' lowcorrsamples.txt combdups.tmp > nodups.tmp
		cp nodups.tmp $gname.tmp
		rm duplicatedsamps.tmp $gname.dupl.tmp $gname.dupl.comb.tmp rmgooddups.tmp nodups.tmp
	elif [ "$duplicates" = 'n' ] || [ -z "$duplicates" ]; then
		printf "  > There are no duplicated samples; skipping duplicates processing. NEW FILE: $gname.tmp \n"
	fi
	
	## d) Plink filters:
	if [ "$plink" = true ] && [ -z "$smap" ]; then
		printf "\n%s\n" "ERROR: PLINK filters requested but no .map file was specified."
		exit 1
	fi
	
	if [ "$plink" = true ]; then
		printf "  ## PLINK FILTERS REQUESTED.\n"
		### d.1) Create .fam file from ancestry file and sort it according to the order from genotypes:
		printf "    > Creating .fam file from ancestry ($filename_pa) and sorting by fish ids from genotypes file ($gname.tmp). NEW FILE: plink_filt.fam \n"
		cp $parent plink_filt.fam
		dos2unix plink_filt.fam
		awk -F' ' '{print $1}' $gname.tmp | tail -n +2 > plink_filt_iids.tmp
		awk -v OFS=' ' 'FNR==NR{arr[$1]=$2;arr2[$1]=$3;next} {print $0,($1 in arr?arr[$1]:0),($1 in arr2?arr2[$1]:0)}' plink_filt.fam plink_filt_iids.tmp > plink_filt.tmp
		rm plink_filt_iids.tmp
		
		### d.2) Use sorted .fam file (plink_filt.tmp) and genotypes to create a .ped file:
		printf "    > Creating .ped file from plink_filt.fam and $gname.tmp. NEW FILE: plink_filt.ped \n"
		printf "      Conversion: -1 => 'NA'; 2 => '2 2'; 1 => '1 2'; 0 => '1 1'; NA => '0 0' GENOTYPES FILE: plink_filt_genots.tmp \n"
		cut -d' ' -f2- $gname.tmp | tail -n +2 | sed 's/-1/NA/g' | sed 's/2/2 2/g; s/1/1 2/g; s/0/1 1/g; s/NA/0 0/g' > plink_filt_genots.tmp
		paste -d' ' plink_filt.tmp plink_filt_genots.tmp > plink_filt.ped
		rm plink_filt.tmp
		
		### d.3) Sort .map file according to the SNP order in genotypes file:
		printf "    > Creating sorted .map file. NEW FILE: plink_filt.map \n"
		#First: Take SNP names (headers from genotypes) and save as rows in temporary file.
		#Second: Sort .map file according to former file:
		cut -d' ' -f2- $gname.tmp | head -n1 | sed 's/\s/\n/g' > snp_order.tmp
		awk -v OFS=' ' 'FNR==NR {arr[$2]=$0;next} $1 in arr {print arr[$1]}' $smap snp_order.tmp > plink_filt.map
		rm snp_order.tmp
		
		### d.4) Execute plink filters:
		filname=$gname"_filt"
		printf "  > PLINK filters running... \n"
		sh $BINPATH/plink_filt.sh $WDIR plink_filt $gname.tmp $filname.tmp
		rm $gname.tmp
		printf "  > PLINK filters finished. NEW FILE: $filname.tmp \n"
	
	elif [ "$plink" = false ]; then
		filname=$gname
		printf "  > PLINK FILTERS NOT REQUESTED. \n"
	fi
	
	## e) Convert genotypes to the appropriate format (blupf90/asrgenomics):
	### e.0) Reassign parameters according to -p flag values:
	if [ "$processing" = 'blupf90' ]; then
		printf "  ## PROCESSING MODE: BLUPF90 \n"
		headers=1
		NAto5='y'
		modif='join'
	elif [ "$processing" = 'asrgenomics' ]; then
		printf "  ## PROCESSING MODE: ASRGenomics. NOTE: this mode is still under development.\n"
		headers=2
		NAto5='y'
	elif [ "$processing" = 'custom' ]; then
		printf "  ## PROCESSING MODE: custom. Please, check carefully the parameter options for recalling possible errors. \n"
	else
		printf "\n%s\n" "ERROR: -p value not valid. Please, use 'blupf90', 'asrgenomics' or 'custom' (default)."
		exit 1
	fi
	
	### e.1) Remove unnecessary rows (usually, column names):
	if [ "$headers" = 1 ]; then
		modname=$filname"_mod"
		printf "    > Removing headers. NEW FILE: $modname.tmp \n"
		tail -n +2 $filname.tmp > $modname.tmp
		rm $filname.tmp
	else
		modname=$filname
		printf "    > Headers not present or removing not requested.\n"
	fi
	
	### e.2) Replace NAs by 5:
	if [ "$NAto5" = 'y' ]; then
		printf "    > Replacing NA/-1 by 5. UPDATING FILE: $modname.tmp \n"
		if [ "$headers" = 0 ] || [ "$headers" = 1 ]; then
			printf "      HEADERS NO OR REMOVE=FALSE\n"
			paste -d' ' <(cut -d' ' -f1 $modname.tmp) <(cut -d' ' -f2- $modname.tmp | sed 's/NA/5/g; s/-1/5/g') > $modname'_5'.tmp
		elif [ "$headers" = 2 ]; then
			printf "      HEADERS YES AND REMOVE=FALSE\n"
			paste -d' ' <(cut -d' ' -f1 $filname.tmp) <(cut -d' ' -f2- $filname.tmp | sed '1 ! s/NA/5/g; 1 ! s/-1/5/g') > $filname'_5'.tmp
		fi
		mv -f $modname'_5'.tmp $modname.tmp
	fi
	
	### e.3) Join SNPs in 1 string (BLUPF90 format):
	if [ ! -z "$modif" ]; then
		printf "    > JOINORSEP function requested: You requested to $modif SNPs from genotypes file.\n"
		if [ "$headers" = 2 ] && [ "$modif" = 'join' ]; then
			printf "    ERROR: You requested JOINORSEP to join SNPs but headers are present. This operation cannot be done.\n"
			exit 1
		elif [ "$processing" = 'asrgenomics' ] && [ "$modif" = 'join' ]; then
			printf "    WARNING: You requested JOINORSEP to join SNPs in ASRGenomics mode. This mode requires SNPs to be separated, so JOINORSEP will not be executed.\n"
		else
			if [ "$processing" = 'asrgenomics' ] && [ "$modif" = 'separate' ]; then
				printf "    WARNING: You requested JOINORSEP to separate SNPs in ASRGenomics mode. The function will be executed, however:
					- If you are using a joined SNP file, you probably lack headers and ASRGenomics won't work properly.
					- If this request is an error, please re-run the script leaving modif blank.\n"	
			fi
			JOINORSEP "$modname.tmp" "$modif" "$extension_o"
			rm variables.mod $modname.tmp
		fi
	else
		mv -f $modname.tmp $modname$extension_o
	fi
	
	genotypes_mod=$(ls -t | grep $gname | head -n1)
	printf "  > FINAL FILE: $genotypes_mod \n\n"

elif [ -z "$processing" ]; then
	# Copy original file to working directory if processing is not required:
	printf "\n\n%s\n" "-------- GENOTYPES WILL NOT BE PROCESSED --------"
	genotypes_mod=$filename_ge
	if [ ! -f $genotypes_mod ]; then
		cp $genotypes $genotypes_mod
	fi
	printf "  > FINAL FILE: $genotypes_mod \n\n"
fi

## H matrix creation: ASRgenomics:
if [ "$Rmatrix" = true ]; then
	echo "----- ASRGenomics ACTIVATED -----"
	#Rscript $BINPATH/ASRgenomics.R $WDIR $genotypes_mod $parent
fi

if [ ! -z "$BF90D" ]; then
	cd $BDIR
	if [ ! -f $genotypes_mod ]; then ln -s $WDIR/$genotypes_mod $genotypes_mod; fi
fi

