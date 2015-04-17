library(xlsx)
library(dplyr)
dat0 <- read.xlsx("Schimmack-JPSP.xlsx", 1), colClasses="numeric")
save(dat0, file="Schimmack.RData")

dat <- dat0 %>% filter(Title != "", Type=="FHT")
head(dat)

res <- c()
for (i in 1:nrow(dat)) {
	header <- paste0(dat$Title[i], " (", dat$Year[i], ") ", dat$Article[i], "/", dat$Study[i])
	
	if (!is.na(dat$F[i])) {
		res <- c(res, paste0(header, ": ", "F(", dat$df1[i], ", ", dat$df2[i], ")=", dat$F[i], "\n"))
	}
	
	if (!is.na(dat$chi.square[i])) {
		res <- c(res, paste0(header, ": ", "chi2(", dat$chi.df[i], ")=", dat$chi.square[i], "\n"))
	}
	
	if (!is.na(dat$r.beta[i])) {
		res <- c(res, paste0(header, ": ", "r(", dat$DF[i], ")=", dat$r.beta[i], "\n"))
	}
}

res <- gsub("#", "", res)
cat(res)

clipboard <- pipe("pbcopy", open="w")
write(res, clipboard)
close(clipboard)