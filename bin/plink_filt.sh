#!/bin/bash

WD=$1
input=$2
oldg=$3
newg=$4

if [ ! -d $WD/PLINK ]; then mkdir -p $WD/PLINK; fi
cd $WD/PLINK

ln -s ../$input.ped $input.ped
ln -s ../$input.map $input.map
ln -s ../$input.fam $input.fam

# 1) Critic filter: Call rates + MAF + HW
printf "\n%s\n" "## FILTER 1: CRITIC FILTER (CR + MAF + HWE)."
filtnm1=$input"_f1"
plink --file $input --no-fid --no-sex --no-pheno --allow-no-sex --nonfounders --chr-set 24 --mind 0.20 --maf 0.05 --geno 0.10 --hwe 0.001 \
      --recode --out $filtnm1

# 2) Outlier filter: 5-nearest neighbours
printf "\n%s\n" "## FILTER 2: OUTLIER SAMPLES (5-NN)."
filtnm2=$input"_f2"
plink --file $filtnm1 --allow-no-sex --nonfounders --chr-set 24 --cluster --neighbour 1 5 --recode --out $filtnm2
mv $filtnm2.log $filtnm2"_1".log

echo "options(warn=-1)
library(tidyverse)
data_5nn <- read.table('$filtnm2.nearest', header = T)
#data_5nn <- read.table('plink_filt_f2.nearest', header = T)
data_5nn.grp <- data_5nn %>% group_by(IID) %>% summarise(avgZ := mean(Z))
data_5nn <- merge(data_5nn, data_5nn.grp, by = 'IID', all.x = T, sort = F)
data_5nn <- data_5nn[order(data_5nn$IID, data_5nn$avgZ),]
neighbour <- -1.5
data_5nn$group <- ifelse(data_5nn$avgZ <= neighbour, 'out', 'keep') %>% factor()
data_5nn_filt <- data_5nn[which(data_5nn$group == 'out'),]
out <- data_5nn_filt[which(!duplicated(data_5nn_filt[,c(2,1)])),c(2,1)]
write.table(out, 'neighbour.nearest.out.txt', col.names = F, row.names = F, quote = F)" > R_nnout.R
Rscript R_nnout.R
rm R_nnout.R

plink --file $filtnm2 --allow-no-sex --chr-set 24 --remove neighbour.nearest.out.txt --recode --out $filtnm2
mv $filtnm2.log $filtnm2"_2".log

# 3) LD filter: linkage desequilibrium
printf "\n%s\n" "## FILTER 3: LINKAGE DESEQUILIBRIUM FILTER."
filtnm3=$input"_f3"
plink --file $filtnm2 --allow-no-sex --make-founders --chr-set 24 --r2 --ld-window-r2 0 --ld-window 100 --ld-window-kb 100 --recode --out $filtnm3
mv $filtnm3.log $filtnm3"_1".log
plink --file $filtnm3 --allow-no-sex --make-founders --chr-set 24 --indep-pairwise 50 5 0.25 --recode --out $filtnm3
mv $filtnm3.log $filtnm3"_2".log
plink --file $filtnm3 --allow-no-sex --chr-set 24 --exclude $filtnm3.prune.out --recode --out $filtnm3
mv $filtnm3.log $filtnm3"_3".log

# 4) EII filter: Identity by descendence
printf "\n%s\n" "## FILTER 4: EII FILTER (IDENTITY BY DESCENDENCE)."
filtnm4=$input"_f4"
plink --file $filtnm3 --allow-no-sex --nonfounders --chr-set 24 --genome --recode --out $filtnm4
mv $filtnm4.log $filtnm4"_1".log
awk -F' ' '$10>0.7{print $1,$2}' $filtnm4.genome | tail -n +2 > $filtnm4"_rm".genome
plink --file $filtnm4 --allow-no-sex --chr-set 24 --remove $filtnm4"_rm".genome --recode --out $filtnm4
mv $filtnm4.log $filtnm4"_2".log
rm $filtnm4"_rm".genome

# 5) Heterozygosis filter: F coefficient estimates
printf "\n%s\n" "## FILTER 5: HETEROZYGOSIS FILTER (F COEFFICIENT)."
filtnm5=$input"_f5"
plink --file $filtnm4 --allow-no-sex --nonfounders --chr-set 24 --het --recode --out $filtnm5
mv $filtnm5.log $filtnm5"_1".log
echo "options(warn=-1)
fname <- '$filtnm5'
hetdat <- read.table(paste0(fname,'.het'), stringsAsFactors = F, header = T)
colnames(hetdat) <- c('FID', 'IID', 'O_HOM', 'E_HOM', 'N_NM', 'F')
hetdat\$rate <- (hetdat\$N_NM-hetdat\$O_HOM)/hetdat\$N_NM
max <- mean(hetdat\$rate) + (3*sd(hetdat\$rate))
min <- mean(hetdat\$rate) - (3*sd(hetdat\$rate))
hetexcl <- hetdat[hetdat\$rate<min|hetdat\$rate>max,c(1,2)]
write.table(hetexcl, paste0(fname,'_rm.het'),col.names = F, row.names = F, quote = F)" > R_hetfilt.R
Rscript R_hetfilt.R
rm R_hetfilt.R
plink --file $filtnm5 --allow-no-sex --chr-set 24 --remove $filtnm5"_rm".het --recode --out $filtnm5
mv $filtnm5.log $filtnm5"_2".log
rm $filtnm5"_rm".het

# 6) FINAL GENOTYPES:
printf "\n%s\n" "## CREATING FINAL GENOTYPES..."
awk -F' ' '{print $2}' $filtnm5.ped > $filtnm5.ids.in
echo "options(warn=-1)
fname <- paste0('$filtnm5')
oldgen <- read.table(paste0('../','$oldg'), stringsAsFactors = F, header = T)
filtids <- read.table(paste0(fname, '.ids.in'), stringsAsFactors = F)[,1]
filtsnps <- read.table(paste0(fname, '.map'), stringsAsFactors = F)[,2]
filtsnps <- gsub('-', '.', filtsnps)
newgen <- oldgen[oldgen[,1] %in% filtids, c(colnames(oldgen)[1], filtsnps)]
write.table(newgen, paste0('../','$newg'), col.names = T, row.names = F, quote = F)" > R_finalgen.R
Rscript R_finalgen.R
rm R_finalgen.R
