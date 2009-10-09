panel.Ct.scatter <-
function (x, y, point.col=point.col, point.pch=point.pch, draw.diag=draw.diag, point.cex=point.cex, ...) 
{
    points(x, y, col=point.col, pch=point.pch, cex=point.cex, ...)
	if (draw.diag)
		abline(0, 1, col="grey", lwd=2)	
}

