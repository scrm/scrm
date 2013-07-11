library(pracma)

ee_tmrca=function(n){
	return (1-1/n);
}

sd_tmrca=function(n){
	i=seq(2,n,1);
	return(sqrt( sum(1/i^2/(i-1)^2)));
}

ee_seg=function(n, rate){
	i=seq(1,n-1,1);
	return(rate*sum(1/i));
}

sd_seg_norecomb=function(n, theta){
	i=seq(1,n-1,1);
	return(sqrt( ee_seg(n,theta) + theta^2*sum(1/i^2)));
}

f2x=function(x){
	return ( (x+18)/(x^2+13*x+18) );
}

fnx=function(x,n){
#	return (n/(2*x*(n-1)));
	i=seq(1,n-1,1);
	return(f2x(x)*sum(1/i^2));
}


f2x_integrand=function(x,rho){
	return ( (rho-x)*f2x(x));
}

fnx_integrand=function(x,rho,n){
	return ( (rho-x)*fnx(x,n=n));
}

sd_seg_recomb=function(n, theta, rho){
#	if (n==2){
		return ( sqrt( ee_seg(n, theta) + theta^2 * 2 /rho^2 *quad(f2x_integrand, xa=0, xb=rho, rho=rho) ) ); # Hein 2005 5.25
#	#	return ( sqrt( theta + theta^2 * 2 /rho^2 *quad(seg_integrand, xa=0, xb=rho, rho=rho) ) ); #Wakeley 2008 7.20, maybe wrong
#	}
#	else{
#		return ( sqrt( ee_seg(n, theta) + theta^2 * 2 /rho^2 *quad(fnx_integrand, xa=0, xb=rho, rho=rho, n=n) ) );
#	}
}




sd_recomb=function(rho,n){
#	if (n==2){
		return( sqrt(ee_seg(n, rho) + 2 * quad(f2x_integrand, xa=0, xb=rho, rho=rho) ) );
#	}
#	else{
#		return( sqrt( ee_seg(n, rho) + 2 * quad(fnx_integrand, xa=0, xb=rho, rho=rho, n=n) ) );
#	}
}


PSk=function(n, k, theta){
	i=seq(2,n,1);
	return ( sum( (-1)^i * choose(n-1, i-1) * (i-1)/(theta+i-1) * (theta/(theta+i-1)^k) ));
}
