source("p-curve-fun.R")

# H0:100 Papers with 5 studies each
res <- data.frame()
for (i in 1:100) {
	studies.H0 <- data.frame(sim_studies(k=5, n=50, sampleFrom="H0"))
	studies.H0$paper_id <- i
	studies.H0$study_id <- 1:nrow(studies.H0)
	studies.H0$true.d <- 0
	res <- rbind(res, data.frame(studies.H0))
	print(i)
}
res.txt <- paste0("H0_", res$paper_id, " (2015) S", res$study_id, ": t(", res$dfs, ")=", round(res$t.values, 3))
writeLines(res.txt, con=con <- file("H0_100x5.txt"))
close(con)


# H1 (d=0.5): 100 Papers with 5 studies each
res <- data.frame()
for (i in 1:100) {
	studies.H1 <- data.frame(sim_studies(k=5, n=50, sampleFrom="H1", d=0.5))
	studies.H1$paper_id <- i
	studies.H1$study_id <- 1:nrow(studies.H1)
	studies.H1$true.d <- 0
	res <- rbind(res, data.frame(studies.H1))
	print(i)
}
res.txt <- paste0("H1_d0.5_", res$paper_id, " (2015) S", res$study_id, ": t(", res$dfs, ")=", round(res$t.values, 3))
writeLines(res.txt, con=con <- file("H1_100x5.txt"))
close(con)


# H0 hack: 100 Papers with 5 studies each
res <- data.frame()
for (i in 1:100) {
	studies.hack <- data.frame(sim_studies(k=5, n=50, sampleFrom="hack", d=0))
	studies.hack$paper_id <- i
	studies.hack$study_id <- 1:nrow(studies.hack)
	studies.hack$true.d <- 0
	res <- rbind(res, data.frame(studies.hack))
	print(i)
}
res.txt <- paste0("hack_", res$paper_id, " (2015) S", res$study_id, ": t(", res$dfs, ")=", round(res$t.values, 3))
writeLines(res.txt, con=con <- file("H0_hack_100x5.txt"))
close(con)