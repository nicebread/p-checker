library(effsize)

x <- c(27, 25, 30, 29, 30, 33, 31, 35)
y <- c(21, 25, 23, 26, 27, 26, 29, 31)
D <- x-y

t.test(D, mu=0)
t.test(x, y)
t.test(x, y, paired=TRUE)

cohen.d(x,y)
tes(t.test(x, y)$statistic, length(x), length(y))

cohen.d(x,y, paired=TRUE)

# t(14)=2.530
# t(7) = 4.513
# t(7) = 4.513