### Fixing the indel immediately following a boundary bug

### Version 2 -- skipping over the INDEL
## The column 2 and 3 values placed right after the boundary should actually be the values AFTER the intron

# accomodae multiple boundaryIndel values -- future-proof vs insertions after more than a single boundary

boundaryIndel <- w[which((w+1) %in% v)] # Identifies the position of the ExonBoundary immediately after which an indel occurs

for(ff in 1:length(boundaryIndel)) { # boundary Indel Loop
 if(length(boundaryIndel[ff]) != 0) { # then we have an indel immediately following a boundary
      vsplitSub <- NA
        for(o in 1:length(vsplit)) {
          if((boundaryIndel[ff]+1) %in% vsplit[[o]]) {vsplitSub <- o}
              toReplace <- vsplit[[vsplitSub]] ### The actual corr_table row 1 position ids that need to be replaced
              numToReplace <- length(vsplit[[vsplitSub]]) ## The number of positions that need to be renamed
            }
        } # close length(boundaryIndel) loop

 ## The values in corr_table rows 2 and 3 need to be replaced with (boundaryIndel-1) . 1:length(toReplace)

row2val <- as.numeric(corr_table[[loci[[i]]]][2,boundaryIndel[ff]-1]) ## the value to extend as indels
row3val <- as.numeric(corr_table[[loci[[i]]]][3,boundaryIndel[ff]-1]) ## the value to exendd as indels

dotDecimal <- 1:numToReplace

corr_table[[loci[[i]]]][2,toReplace] <- paste(rep(row2val,numToReplace),dotDecimal,sep=".")
corr_table[[loci[[i]]]][3,toReplace] <- paste(rep(row3val,numToReplace),dotDecimal,sep=".")

corr_table[[loci[[i]]]][,toReplace]

} # close boundaryIndel loop
