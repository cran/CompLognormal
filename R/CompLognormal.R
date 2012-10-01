

dcomplnorm<-function(x, spec, sigma=1, theta=1, ...)
{	
	add.pars	<- match.call(expand.dots=FALSE)$...
	add.pars.val<- as.numeric(paste(match.call(expand.dots=FALSE)$...))

	fun1.args	<- function(y) list(y, meanlog=mu, sdlog=sigma)
	fun1		<- function(d) {function(z) do.call(paste(d,"lnorm",sep=""), fun1.args(z))}
	fun2.args	<- function(y)
	{
		fun.args						<- list()
		fun.args						<- setNames(fun.args[1:(length(add.pars)+1)],c("", names(add.pars)))
		fun.args[[1]]					<- y
		fun.args[2:(length(add.pars)+1)]		<- add.pars.val
		return(fun.args)
	}
	fun2		<- function(d) {function(z) do.call(paste(d,spec,sep=""), fun2.args(z))}

	mu		<- log(theta) + sigma^2 + theta*sigma^2*grad(fun2("d"),theta) / fun2("d")(theta)
	phi		<- ( fun1("d")(theta)*(1-fun2("p")(theta) ) / ( fun2("d")(theta)*fun1("p")(theta) ) )

	pdf		<- ifelse( x<0|x==0, 0,
				ifelse( x<theta | x==theta,
				(1/(1+phi)) / fun1("p")(theta)*fun1("d")(x),
				(phi/(1+phi)) / (1-fun2("p")(theta))*fun2("d")(x)
				)
			   )
	return(pdf)
}


pcomplnorm<-function(x, spec, sigma=1, theta=1, ...)
{	

	add.pars	<- match.call(expand.dots=FALSE)$...
	add.pars.val<- as.numeric(paste(match.call(expand.dots=FALSE)$...))

	fun1.args	<- function(y) list(y, meanlog=mu, sdlog=sigma)
	fun1		<- function(d) {function(z) do.call(paste(d,"lnorm",sep=""), fun1.args(z))}
	fun2.args	<- function(y)
	{
		fun.args						<- list()
		fun.args						<- setNames(fun.args[1:(length(add.pars)+1)],c("", names(add.pars)))
		fun.args[[1]]					<- y
		fun.args[2:(length(add.pars)+1)]		<- add.pars.val
		return(fun.args)
	}
	fun2		<- function(d) {function(z) do.call(paste(d,spec,sep=""), fun2.args(z))}

	mu		<- log(theta) + sigma^2 + theta*sigma^2*grad(fun2("d"),theta) / fun2("d")(theta)
	phi		<- ( fun1("d")(theta)*(1-fun2("p")(theta) ) / ( fun2("d")(theta)*fun1("p")(theta) ) )

	cdf		<- ifelse( x<0 | x==0, 0,
				ifelse( x<theta | x==theta,
				(1/(1+phi))*fun1("p")(x) / fun1("p")(theta),
				(1/(1+phi)) + (phi/(1+phi))*(fun2("p")(x)-fun2("p")(theta))/(1-fun2("p")(theta))
				)
			   )
	return(cdf)
}



qcomplnorm<-function(p, spec, sigma=1, theta=1, ...)
{	

	add.pars	<- match.call(expand.dots=FALSE)$...
	add.pars.val<- as.numeric(paste(match.call(expand.dots=FALSE)$...))

	fun1.args	<- function(y) list(y, meanlog=mu, sdlog=sigma)
	fun1		<- function(d) {function(z) do.call(paste(d,"lnorm",sep=""), fun1.args(z))}
	fun2.args	<- function(y)
	{
		fun.args						<- list()
		fun.args						<- setNames(fun.args[1:(length(add.pars)+1)],c("", names(add.pars)))
		fun.args[[1]]					<- y
		fun.args[2:(length(add.pars)+1)]		<- add.pars.val
		return(fun.args)
	}
	fun2		<- function(d) {function(z) do.call(paste(d,spec,sep=""), fun2.args(z))}

	mu		<- log(theta) + sigma^2 + theta*sigma^2*grad(fun2("d"),theta) / fun2("d")(theta)
	phi		<- ( fun1("d")(theta)*(1-fun2("p")(theta) ) / ( fun2("d")(theta)*fun1("p")(theta) ) )

	qf		<- ifelse( p<0 | p>1, NaN,
				c( fun1("q")( p[p<(1/(1+phi)) | p==(1/(1+phi))]*(1+phi) * fun1("p")(theta) ),
				fun2("q")( (p[p>(1/(1+phi))]*(1+phi)-1)/phi * (1-fun2("p")(theta))
                + fun2("p")(theta) ) )
			   )
	return(qf)
}



rcomplnorm<-function(n, spec, sigma=1, theta=1, ...)
{	
	add.pars	<- match.call(expand.dots=FALSE)$...
	add.pars.val<- as.numeric(paste(match.call(expand.dots=FALSE)$...))


	fun1.args	<- function(y) list(y, meanlog=mu, sdlog=sigma)
	fun1		<- function(d) {function(z) do.call(paste(d,"lnorm",sep=""), fun1.args(z))}
	fun2.args	<- function(y)
	{
		fun.args						<- list()
		fun.args						<- setNames(fun.args[1:(length(add.pars)+1)],c("", names(add.pars)))
		fun.args[[1]]					<- y
		fun.args[2:(length(add.pars)+1)]		<- add.pars.val
		return(fun.args)
	}
	fun2		<- function(d) {function(z) do.call(paste(d,spec,sep=""), fun2.args(z))}

	mu		<- log(theta) + sigma^2 + theta*sigma^2*grad(fun2("d"),theta) / fun2("d")(theta)
	phi		<- ( fun1("d")(theta)*(1-fun2("p")(theta) ) / ( fun2("d")(theta)*fun1("p")(theta) ) )

	rand		<- runif(n)

	rf		<- ifelse( rand<0 | rand>1, 0,
				c( fun1("q")( rand[rand<(1/(1+phi)) | rand==(1/(1+phi))] * (1+phi)
                     * fun1("p")(theta) ),
				fun2("q")( (rand[rand>(1/(1+phi))]*(1+phi)-1)/phi * (1-fun2("p")(theta))
                    + fun2("p")(theta) ) )
			   )
	return(rf)
}



nlogl.complnorm<-function(p,x,spec,...)
{	
	sigma		<- p[1]
	theta		<- p[2]
	add.pars	<- paste(match.call(expand.dots=FALSE)$...)

	tt						<- 1.0e20

if(all(p>0)){

	fun1.args	<- function(y) list(y, meanlog=mu, sdlog=sigma)
	fun1		<- function(d) {function(z) do.call(paste(d,"lnorm",sep=""), fun1.args(z))}
	fun2.args	<- function(y)
	{
		fun.args						<- list()
		fun.args						<- setNames(fun.args[1:(length(add.pars)+1)],c("", add.pars))
		fun.args[[1]]					<- y
		fun.args[2:(length(add.pars)+1)]		<- p[3:(length(add.pars)+2)]
		return(fun.args)
	}
	fun2		<- function(d) {function(z) do.call(paste(d,spec,sep=""), fun2.args(z))}

	mu		<- log(theta) + sigma^2 + theta*sigma^2*grad(fun2("d"),theta) / fun2("d")(theta)
	phi		<- fun1("d")(theta)*(1-fun2("p")(theta)) / ( fun2("d")(theta)*fun1("p")(theta) )

	tt		<- -sum(log( (1/(1+phi))/fun1("p")(theta)*fun1("d")(x[x>0&x<theta])))
        - sum(log( (phi/(1+phi))/(1-fun2("p")(theta))*fun2("d")(x[x>theta]) ))

}

return(tt)
}
