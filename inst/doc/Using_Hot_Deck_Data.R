## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, echo=FALSE, include=FALSE-------------------------------------------------------------
options(useFancyQuotes=FALSE, width=100)

## ----echo=T, include=T----------------------------------------------------------------------------
library(hot.deck)
data(isq99)
out <- hot.deck(isq99, sdCutoff=3, IDvars = c("IDORIGIN", "YEAR"))

## ----numdonors, echo=T, include=T-----------------------------------------------------------------
numdonors <- sapply(out$donors, length)
numdonors <- sapply(out$donors, length)
numdonors <- ifelse(numdonors > 5, 6, numdonors)
numdonors <- factor(numdonors, levels=1:6, labels=c(1:5, ">5"))
table(numdonors)

## ----tscslag, echo=T, include=T-------------------------------------------------------------------
tscslag <- function(dat, x, id, time){
  obs <- apply(dat[, c(id, time)], 1, paste, collapse=".")
  tm1 <- dat[[time]] - 1
  lagobs <- apply(cbind(dat[[id]], tm1), 1, paste, collapse=".")
  lagx <- dat[match(lagobs, obs), x]
}
for(i in 1:length(out$data)){
  out$data[[i]]$lagAI <- tscslag(out$data[[i]], "AI", "IDORIGIN", "YEAR")
  out$data[[i]]$lagPCGNP <- tscslag(out$data[[i]], "PCGNP", "IDORIGIN", "YEAR")
  out$data[[i]]$lagLPOP <- tscslag(out$data[[i]], "LPOP", "IDORIGIN", "YEAR")
}

## ----pcgchange, echo=T, include=T-----------------------------------------------------------------
for(i in 1:length(out$data)){
  out$data[[i]]$pctchgPCGNP <- with(out$data[[i]], c(PCGNP-lagPCGNP)/lagPCGNP)
  out$data[[i]]$pctchgLPOP <- with(out$data[[i]], c(LPOP-lagLPOP)/lagLPOP)
}

## ----mods, echo=T, include=T----------------------------------------------------------------------
# initialize list
out <- hd2amelia(out)
results <- list()
# loop over imputed datasets
for(i in 1:length(out$imputations)){
    results[[i]] <- lm(AI ~ lagAI + pctchgPCGNP + PCGNP + pctchgLPOP + LPOP + MIL2 + LEFT +
    BRIT + POLRT + CWARCOW + IWARCOW2, data=out$imputations[[i]])
}
summary(mitools::MIcombine(results))

## ----conv, echo=T, include=T----------------------------------------------------------------------
out.mids <- miceadds::datalist2mids(out$imputations)
s <- summary(mice::pool(mice::lm.mids(AI ~ lagAI + pctchgPCGNP + PCGNP + pctchgLPOP + LPOP + MIL2 + LEFT +
BRIT + POLRT + CWARCOW + IWARCOW2, data=out.mids)))
print(s, digits=4)

