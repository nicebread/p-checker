x = "1: t(88)=2.1; crit = .10"
x = "r(147)=.246"
x = "F(1,100)=9.1"
x = "f(2,210)=4.45"
x = "Z=3.45"
x = "chi2(1)=9.1"
x = "r(77)=.47"

parse_ES1("3: t(88)=2.1; crit = .10; p < .04", round_up=TRUE)
parse_ES1("r(147)=.246")
parse_ES1("test: F(1,100)=9.1")
parse_ES1("f(2,210)=4.45")
parse_ES1("Z=3.45")
parse_ES1("chi2(1)=9.1")
parse_ES1("r(77)=.47")


x <- c("
a1: t(88)=2.1; crit=.10
a1: r(147)=.246
a1: F(1,100)=9.1; crit = .08
f(2,210)=4.45
Z=3.45
chi2(1)=9.1
r(77)=.47")


x <- c("
t(88)=2.1
r(147)=.246
F(1,100)=9.1
f(2,210)=4.45
Z=3.45
chi2(1)=9.1
r(77)=.47
chi2(2)=8.74
")

tbl2 <- parse_ES(txt)


# Test wrong syntax

parse_ES1("3: t(88, 2)=2.1; crit = .10; p < .04", round_up=TRUE)
r <- parse_ES("3: t(88, 2)=2.1; crit = .10; p < .04")
attr(r, "warnings")
