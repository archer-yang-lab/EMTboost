############################################################################
#      F(x) generator, Based on Friedman(1999)
############################################################################

"Fx.generator"<-function(X,N,P)
{
		a <- 0.1
		b <- 2.0
		r <- rexp(20, rate = 1/2)
		n_l <- floor(2.5+r)
		g <- matrix(0,N,20)
		n_l[n_l > P] <- P
		a_l <- runif(20,-1,1)

		for (i in 1:20)
		{
			tstMat <- array(rnorm(n_l[i]), dim=c(n_l[i],n_l[i]))
			U <- qr.Q(qr(tstMat))                       # generate random orthogonal matrix
			d<-runif(n_l[i],a,b)                        # eigenvalues
			D<-diag(d^2,n_l[i],n_l[i])
		    V<-U%*%D%*%t(U)                                # get V_l

			mu <- rnorm(n_l[i], 0, 1)                  # get mu_l
		    oo <- sample(1:P, n_l[i])
			z <- as.matrix(X[ ,oo])       # get z_l
			
			for (j in 1:N){
				c <- (z[j, ]-mu)%*%V%*%(z[j, ]-mu)
				g[j,i] <- exp(-0.5*c)
			}
		}
		F_x <- g %*% a_l
		F_x
}

# g1 = function(x) {
# 	2 * sin(pi*x)
# }
# 
# g2 = function(x) {
# 	exp(2 * x)
# }
# 
# g3 = function(x) {
# 	x^11 * (10*(1-x))^6 / 5 + 10^4*x^3*(1-x)^10
# }
# 
# g4 = function(x){
# 	0*x
# }
# 
# "Gx.generator"<-function(X)
# {
# 	y1 = g1(X[,1])
# 	y2 = g2(X[,2])
# 	y3 = g3(X[,3])
# 	y4 = g4(X[,4])
# 	Gx = 0.2 * y1 + 0.1 * y2 + 0.05 * y3 + y4
# 	Gx
# }

j1 = function(x) {
	2*x-ifelse(x>.5,1,0)
}

j2 = function(x) {
	ifelse((x>.75), (16/3)*x*(x-1)^2, 0) + ifelse((x<=.5), 4*x^2*(3-4*x), 0) + ifelse((x<=.75)&&(x>.5), (4/3)*x*(4*x^2-10*x+7)-3/2, 0)
}
	
j3 = function(x) {
	2-2*abs(x-.26)^.2*ifelse(x<=.26,1,0)-2*abs(x-.26)^.6*ifelse(x>.26,1,0)+ifelse(x>=.78,1,0)
}	
	

j4 = function(x) {
	2*sin(4*pi*x)-6*abs(x-.4)^.3-.5*sign(.7-x)
}	

# X = rep(seq(0,1,0.01),4)
# X = matrix(X,101,4)
"Jx.generator"<-function(X)
{
	y1 = sapply(X[,1],j1)
	y2 = sapply(X[,2],j2)
	y3 = sapply(X[,3],j3)
	y4 = sapply(X[,4],j4)
	Jx = y1 + y2 + y3 + y4
	Jx
}

h1 = function(x) {
	x[,1] * sin(4*pi*x[,2])
}
h2 = function(x) {
	sin(5*x[,1])*cos(5*(x[,2]+1)^2)
}

"hx.generator"<-function(X)
{
	y1 = h1(X[,1:2])
	y2 = h2(X[,3:4])
	hx = y1 + y2
	hx
}