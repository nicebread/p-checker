library(rio)
dat <- import("VadilloHardwickeShanks_effect_sizes edited.xlsx")


dat$id <- 1:nrow(dat)

plot(dat[, "d"], dat[, "p-checker d"])

# relative error:
dat$relErr <- dat[, "p-checker d"]/dat[, "d"]

dat[order(dat$relErr), c("Paper title", "# study", "id", "relErr")]

dat[47, ]