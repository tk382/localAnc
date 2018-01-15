#those who have too many missing values

setwd('/Volumes/im-lab/nas40t2/tae/differentialeQTLs/multitissue/localAnc/400genes_result1/')
L = read.table('../../../data/finalL_smaller.txt', header = TRUE)
total = colnames(L)
names = colnames(L)
names = list.files()
tmp = strsplit(names, '[.]')
part1 = sapply(tmp, function(x) x[1])
part2 = sapply(tmp, function(x) x[2])
combined = paste0(part1, '.',part2)

ind = which(!total%in%combined)
residual = total[ind]
write.table(residual,'genelist_for_400genes_try2.txt',col.names=FALSE,row.names=FALSE,quote=FALSE)


