######################################################################
#   Version  03.03.2012 * based on Hauke's adjusted expression dataset
######################################################################


#setwd("")


# load dataset

map = c(1,2,3,4,5,6,8,16,17,19,7,9,20,21,10,13,31,33,34,23,24,25,26,27,18,22,11,12,14,32,35,36,37,28,29,30,38,39,40,41,42,43,44,45,46,47,48,49) # map from Hauke to Sebastian W.

data = read.table("../Normalization_20120329/Adjusted_expression.xls", header=TRUE, sep="\t")
data = data[c(-117,-11671),] # these rows left out for now since problematic in Hauke's annotation
eset = data[,map]
rownames(eset) = data[,51]


# average multiple measurements (arithmetic)

duplicate_avg = log2( (2^eset[,c(1,25,13,14)] + 2^eset[,c(45,46,47,48)]) * 0.5 ) # 0h, 6h, 12h, 24h (Marisa, Barbara)
triplicate_avg = matrix(log2( (( (2^eset[,27] + 2^eset[,15]) / 2.0) + 2^eset[,28]) / 2.0)) # 12h ((technical, technical) biological)


# define experimental series

# NGF
eset.control = cbind(duplicate_avg[,1], eset[,c(8:9)], duplicate_avg[,2], eset[,10], duplicate_avg[,c(3:4)], eset[,26]) # 0h, 2h, 4h, 6h, 8h, 12h, 24h, 48h
eset.ngf = cbind(eset.control[,1], eset[,c(2:6,11,7,12)], triplicate_avg, eset[,c(16,29)]) # 0h, 0.5h, 1h, 2h, 3h, 4h, 5h, 6h, 8h, 12h, 24h, 48h
eset.inhibitor = cbind(eset.control[,1],eset[,c(17,30,18:19,31:33)]) # 0h, 1h, 2h, 4h, 8h, 12h, 24h, 48h
eset.ngf_inhibitor = cbind(eset.control[,1],eset[,c(20:24,34:36)]) # 0h, 1h, 2h, 4h, 6h, 8h, 12h, 24h, 48h

# EGF
eset.egf = cbind(eset.control[,1],eset[,c(37:44)]) # 0h, 1h, 2h, 3h, 4h, 6h, 8h, 12h, 24h


# define timepoints

timepoints.control = c(0.0,2.0,4.0,6.0,8.0,12.0,24.0,48.0)
timepoints.ngf = c(0.0,0.5,1.0,2.0,3.0,4.0,5.0,6.0,8.0,12.0,24.0,48.0)
timepoints.inhibitor = c(0.0,1.0,2.0,4.0,8.0,12.0,24.0,48.0)
timepoints.ngf_inhibitor = c(0.0,1.0,2.0,4.0,6.0,8.0,12.0,24.0,48.0)
timepoints.egf = c(0.0,1.0,2.0,3.0,4.0,6.0,8.0,12.0,24.0)

all_timepoints_48h = c(0.0,0.5,1.0,2.0,3.0,4.0,5.0,6.0,8.0,12.0,24.0,48.0)
all_timepoints_24h = c(0.0,0.5,1.0,2.0,3.0,4.0,5.0,6.0,8.0,12.0,24.0)


# interpolate

# NGF
eset.control_48h.interp = as.data.frame(t(apply((eset.control),1,function(x){i = approx(timepoints.control,x,xout=all_timepoints_48h); return(i$y)})))
eset.inhibitor_48h.interp = as.data.frame(t(apply((eset.inhibitor),1,function(x){i = approx(timepoints.inhibitor,x,xout=all_timepoints_48h); return(i$y)})))
eset.ngf_inhibitor_48h.interp = as.data.frame(t(apply((eset.ngf_inhibitor),1,function(x){i = approx(timepoints.ngf_inhibitor,x,xout=all_timepoints_48h); return(i$y)})))
eset.ngf_48h.interp = as.data.frame(t(apply((eset.ngf),1,function(x){i = approx(timepoints.ngf,x,xout=all_timepoints_48h); return(i$y)})))

# EGF
eset.egf_24h.interp = as.data.frame(t(apply((eset.egf),1,function(x){i = approx(timepoints.egf,x,xout=all_timepoints_24h); return(i$y)})))
eset.ngf_24h.interp = as.data.frame(t(apply((eset.ngf),1,function(x){i = approx(timepoints.ngf,x,xout=all_timepoints_24h); return(i$y)})))
eset.control_24h.interp = as.data.frame(t(apply((eset.control),1,function(x){i = approx(timepoints.control,x,xout=all_timepoints_24h); return(i$y)})))


# normalize

# NGF
eset.normed.ngf_48h = eset.ngf_48h.interp - eset.control_48h.interp
eset.normed.ngf_inhibitor_48h = eset.ngf_inhibitor_48h.interp - eset.inhibitor_48h.interp

# EGF
eset.normed.egf_24h = eset.egf_24h.interp - eset.control_24h.interp
eset.normed.ngf_24h = eset.ngf_24h.interp - eset.control_24h.interp 



mysym <- as.vector(as.character(mget(rownames(a),org.Rn.egSYMBOL,ifnotfound=NA)))






# fix headers

colnames(eset.normed.ngf_48h) = all_timepoints_48h
colnames(eset.normed.ngf_inhibitor_48h) = all_timepoints_48h
colnames(eset.normed.egf_24h) = all_timepoints_24h
colnames(eset.normed.ngf_24h) = all_timepoints_24h


# export normalized sets for MDS analysis

#write.table(eset.normed.ngf_48h, file="eset.normed.ngf_48h.csv", sep=",", quote=FALSE)
#write.table(eset.normed.ngf_inhibitor_48h, file="eset.normed.ngf_inhibitor_48h.csv", sep=",", quote=FALSE)
#write.table(eset.normed.egf_24h, file="eset.normed.egf_24h.csv", sep=",", quote=FALSE)
#write.table(eset.normed.ngf_24h, file="eset.normed.ngf_24h.csv", sep=",", quote=FALSE)



#######################
#   Version  20.03.2012
#######################


# normalize to time point zero only (for qPCR comparisons)

# NGF
eset.normed_t0.ngf_48h = eset.ngf_48h.interp - eset.ngf_48h.interp[,1]
eset.normed_t0.ngf_inhibitor_48h = eset.ngf_inhibitor_48h.interp - eset.ngf_inhibitor_48h.interp[,1]

# EGF
eset.normed_t0.egf_24h = eset.egf_24h.interp - eset.egf_24h.interp[,1]
eset.normed_t0.ngf_24h = eset.ngf_24h.interp - eset.ngf_24h.interp[,1]


# fix headers

colnames(eset.normed_t0.ngf_48h) = all_timepoints_48h
colnames(eset.normed_t0.ngf_inhibitor_48h) = all_timepoints_48h
colnames(eset.normed_t0.egf_24h) = all_timepoints_24h
colnames(eset.normed_t0.ngf_24h) = all_timepoints_24h


# export normalized sets

write.table(eset.normed_t0.ngf_48h, file="eset.normed_t0.ngf_48h.csv", sep=",", quote=FALSE)
write.table(eset.normed_t0.ngf_inhibitor_48h, file="eset.normed_t0.ngf_inhibitor_48h.csv", sep=",", quote=FALSE)
write.table(eset.normed_t0.egf_24h, file="eset.normed_t0.egf_24h.csv", sep=",", quote=FALSE)
write.table(eset.normed_t0.ngf_24h, file="eset.normed_t0.ngf_24h.csv", sep=",", quote=FALSE)