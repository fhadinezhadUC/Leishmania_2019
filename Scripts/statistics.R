bothdf <- genefile[genefile$foundby=="both",]
tsehasintron <- bothdf$tseintronbegin!=0
table(hasintron)
89
arahasintron <- bothdf$araintronbegin!="nointron"
106
both_hasintondf <- bothdf[tsehasintron & arahasintron,]

# 89 genes with intron by both gene finders
table(paste(both_hasintondf$arafunc,both_hasintondf$tsefunc,sep = "|"))
# N|Y Y|Y 
# 11  78 
# (manualy looked at the alignment) These 11 genes excluding their introns look exactly like other "Y" genes marked by both tse and ara 

table(paste(bothdf[!tsehasintron & arahasintron,]$arafunc,bothdf[!tsehasintron & arahasintron,]$tsefunc,sep = "|"))
# L|E L|M O|M V|V 
#  1   9   2   5 
# for the "V" tse reports as insertion but ara reports as intron 
# for L|M they look like other "M" genes. tse reports two insertions at 37 but ara reports intron size 4 before anticodon ((((( iiii AAA  )))))
#GCGCAGTTGGTCCAGTGGGAGAATTCCCGCATCATTaaAGCGGGGGGCCCGGGTTCGATTCCCGGACTGCGCA
#>>>>>>>..>>.>.......<.<<.>>>>>.........<<<<<....>>>>>.......<<<<<<<<<<<<.
#gcgcagttggtccagtgggagaattcccgcatcattaaagcggggggcccgggttcgattcccggactgcgca
#(((((((  ((.(ddddddd).))((((( iiii AAA  )))))   (((((ttttttt))))))))))))

# for O|M tse reports insertion after position 37 but ara reports intron at 34 (((((  AiiAA  ))))) these two sequences do not look like other sequences marked as M! 
# L|E this one gene looks more like M than other genes E or L!!! the tse score is 48! ara 113!

table(paste(bothdf[bothdf$tsefunc!=bothdf$arafunc,]$tsescore,bothdf[bothdf$tsefunc!=bothdf$arafunc,]$arafunc,bothdf[bothdf$tsefunc!=bothdf$arafunc,]$tsefunc,bothdf[bothdf$tsefunc!=bothdf$arafunc,]$arascore,sep = "|"))

#43.5|L|?|113.6 43.5|S|R|113.6 47.8|O|M|112.4 48.6|L|E|113.1 49.9|D|I|115.2 50.0|L|M|112.7 50.5|L|M|112.7 54.4|W|G|113.7 60.0|N|Y|101.9 72.3|N|Y|108.4 
#3              1              2              1              3              8              1              3              1             10 
table(paste(bothdf[bothdf$arafunc!=bothdf$tsefunc,]$arafunc,bothdf[bothdf$arafunc!=bothdf$tsefunc,]$tsefunc,sep = "|"))
# D|I L|? L|E L|M N|Y O|M S|R W|G 
#  3   3   1   9  11   2   1   3 

table(both_hasintondf$araintronbegin)
# 36 37 38 
# 10  1 78 

both_hasintondf[both_hasintondf$araintronbegin!=38,]
# Intersection set: 89 genes have intron and they are all type Tyr(Y) in TSE but 11 of them are marked as N in ara!

aradf <- genefile[genefile$foundby=="ara",]
aradf[aradf$araintronbegin=="nointron",]$arageness
table(aradf$araintronbegin)
# 27       28       29       30       31       32       33       34       35       36       37       38       39       40       41       42       43       44 nointron 
# 3       15       12       24       27       40       40       64       58       34       62       35       86       18        8        1        4        1      209 
hist(as.numeric(aradf$arascore))
# 
