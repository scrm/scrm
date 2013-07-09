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

seg_integrand=function(x,rho){
	return ( (rho-x)*f2x(x));
}

sd_seg_recomb=function(n, theta, rho){
	return ( sqrt( theta + theta^2 * 2 /rho^2 *quad(seg_integrand, xa=0, xb=rho, rho=rho) ) );

}

fnx=function(x,n){
	return (n/2/x/(n-1));
#	i=seq(1,n-1,1);
#	return(f2x(x)*sum(1/i^2));
}

recomb_integrand=function(x,rho,n){
	return ( (rho-x)*fnx(x,n=n));
}

sd_recomb=function(rho,n){
	if (n==2){
		return( sqrt(ee_seg(n, rho) + 2 * quad(seg_integrand, 0, rho, rho=rho) ) );
	}
	else{
		return( sqrt( ee_seg(n, rho) + 2 * quad(recomb_integrand, 0, rho, rho=rho, n=n) ) );
	}
}
