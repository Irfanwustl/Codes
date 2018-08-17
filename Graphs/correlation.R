> rsq(DATA[,78],DATA[,79])
[1] 0.04511633
> scatterplot(mpg ~ wt | cyl, data=mtcars,
+   xlab="Weight of Car", ylab="Miles Per Gallon",
+    main="Enhanced Scatter Plot",
+    labels=row.names(mtcars))
Error in scatterplot(mpg ~ wt | cyl, data = mtcars, xlab = "Weight of Car",  :
  could not find function "scatterplot"
> scatterplot(DATA[,78],DATA[,79])
Error in scatterplot(DATA[, 78], DATA[, 79]) :
  could not find function "scatterplot"
> plot(DATA[,78],DATA[,79])
> plot(log(DATA[,78)],log(DATA[,79]))
Error: unexpected ')' in "plot(log(DATA[,78)"
> plot(log(DATA[,78)+1],log(DATA[,79]+1))
Error: unexpected ')' in "plot(log(DATA[,78)"
> plot(log(DATA[,78]+1),log(DATA[,79]+1))
> abline(lm(x ~ y)
+
+ )
Error in eval(predvars, data, env) : object 'x' not found
> abline(lm(log(DATA[,78]+1) ~ log(DATA[,79]+1)))
Error in int_abline(a = a, b = b, h = h, v = v, untf = untf, ...) :
  plot.new has not been called yet
> abline(x=y)
Error in int_abline(a = a, b = b, h = h, v = v, untf = untf, ...) :
  object 'y' not found
> abline(1)
Error in int_abline(a = a, b = b, h = h, v = v, untf = untf, ...) :
  plot.new has not been called yet
> abline(x=x)
Error in int_abline(a = a, b = b, h = h, v = v, untf = untf, ...) :
  object 'x' not found
> abline(0,1)
Error in int_abline(a = a, b = b, h = h, v = v, untf = untf, ...) :
  plot.new has not been called yet
> plot(log(DATA[,78]+1),log(DATA[,79]+1))
> abline(0,1)
> cor.test(DATA[,78],DATA[,79])

	Pearson's product-moment correlation

data:  DATA[, 78] and DATA[, 79]
t = 4.008, df = 340, p-value = 7.528e-05
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 0.1088066 0.3114410
sample estimates:
     cor
0.212406

> t.test(DATA[,78],DATA[,79])

	Welch Two Sample t-test

data:  DATA[, 78] and DATA[, 79]
t = 0.029793, df = 675.15, p-value = 0.9762
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -1010.972  1042.125
sample estimates:
mean of x mean of y
 3024.796  3009.220

> pdf("H_N.pdf")
> plot(log(DATA[,78]+1),log(DATA[,79]+1))
> abline(0,1)
> dev.off()
quartz
     2
> pdf("H_N.pdf")
> plot(log(DATA[,78]+1),log(DATA[,79]+1),xlab="HTG",ylab="Nanostring",main="HTG vs Nanostring")
> abline(0,1)
> dev.off()
quartz
     2
> ?plot
> pdf("H_N.pdf")
> plot(log(DATA[,78]+1),log(DATA[,79]+1),xlab="HTG",ylab="Nanostring",main="HTG vs Nanostring",xlim=c(0,12),ylim=c(0,12))
> abline(0,1)
> dev.off()
quartz
     2
> pdf("H_N.pdf")
> plot(log(DATA[,78]+1),log(DATA[,79]+1),xlab="HTG",ylab="Nanostring",main="HTG vs Nanostring",xlim=c(3.5,12),ylim=c(3.5,12))
> abline(0,1)
> dev.off()
quartz
