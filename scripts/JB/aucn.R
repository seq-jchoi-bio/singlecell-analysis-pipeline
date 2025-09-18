args <- commandArgs(trailingOnly=TRUE) 

########## calculate auprc #########
auPRC <- function (perf) {
	rec <- perf@x.values
	prec <- perf@y.values
	result <- list()
	for (i in 1:length(perf@x.values)) {
	result[i] <- list(sum((rec[[i]][2:length(rec[[i]])] - rec[[i]][2:length(rec[[i]])-1])*prec[[i]][-1]))
	}
	return(result)
}

########## plot ROC and PR-Curve #########
rocprc <- function(x, output="rocprc.pdf", plotting=T) {
	sink(NULL,type="message")
	options(warn=-1)
	suppressMessages(suppressWarnings(library('ROCR')))
	svmresult <- data.frame(x)
    if (length(x)==5) {
	colnames(svmresult) <- c("Seqid","Label", "PLabel", "Pred", "CV")
    }
    if (length(x)==4) {
	colnames(svmresult) <- c("Seqid","Pred", "Label", "CV")
    }

	linewd <- 1
	wd <- 4
	ht <- 4
	fig.nrows <- 1 
	fig.ncols <- 2
	pt <- 10
	cex.general <- 1 
	cex.lab <- 0.9
	cex.axis <- 0.9
	cex.main <- 1.2
	cex.legend <- 0.8
    
    preds<-svmresult$Pred
    labs<-svmresult$Label

	pred <- prediction(preds, labs)
	perf_roc <- performance(pred, 'tpr', 'fpr')
	perf_prc <- performance(pred, 'prec', 'rec')

	perf_auc <- performance(pred, 'auc')
	prcs <- auPRC(perf_prc)
    avgauc <- perf_auc@y.values[[1]]
    avgprc <- prcs[[1]]

    #cat(paste("auROC:", avgauc, "\n", sep=""))
    #cat(paste("auPRC:", avgprc, "\n", sep=""))

    if (plotting) {
	pdf(output, width=wd*fig.ncols, height=ht*fig.nrows)
	#png(output, width=wd*fig.ncols, height=ht*fig.nrows, unit="in", res=150)

	par(xaxs="i", yaxs="i", mar=c(3.5,3.5,2,2)+0.1, mgp=c(2,0.8,0), mfrow=c(fig.nrows, fig.ncols))

	plot(perf_roc, colorize=T, main="ROC curve", spread.estimate="stderror",
	xlab="1-Specificity", ylab="Sensitivity", cex.lab=1.2)
	text(0.2, 0.1, paste("AUC=", format(avgauc, digits=3, nsmall=3)))
	cat(paste("auROC=", format(avgauc, digits=3, nsmall=3), '\n'))

	plot(perf_prc, colorize=T, main="P-R curve", spread.estimate="stderror",
	xlab="Recall", ylab="Precision", cex.lab=1.2, xlim=c(0,1), ylim=c(0,1))
	text(0.2, 0.1, paste("AUC=", format(avgprc, digits=3, nsmall=3)))
	cat(paste("auPRC=", format(avgprc, digits=3, nsmall=3), '\n'))

	dev.off()
    }

    return(list(auroc=avgauc, auprc=avgprc))
}

############## main function #################
for (fn in args) {
    d <- read.table(fn)
    res<-rocprc(d, plotting=F)
    cat(paste(fn, nrow(d)/2, paste(res$auroc,"\n", sep=""), sep="\t"))
}
