panel.Ct.cor <-
function(x, y, text.cex=text.cex, Ct.max=Ct.max, ...)
{
	usr <- par("usr"); on.exit(par(usr))
	par(usr = c(0, 1, 0, 1))
	index	<- apply(cbind(x,y), 1, function(z) any(z>Ct.max))
	r <- abs(cor(x[!index], y[!index]))
	txt <- format(c(r, 0.123456789), digits=2)[1]
	text(0.5, 0.5, txt, cex=text.cex, ...)
}

