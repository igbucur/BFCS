/* R callable functions for trait analysis */

/************************   NOTE   ************************

   All these routines expect to be called with .C(..., DUP=FALSE)
   The result argument MUST be specifically allocated for this purpose
   by the calling routine.  If not, it will interfere with R's memory
   management, which doesn't always make copies of data until needed.

************************************************************/

/* -I$RHOME when compiling */
#include <math.h>
#include <assert.h>
#include "R.h"
#include "Rdefines.h"
#include "Rinternals.h"

/* assume that arguments have the right type already!!! */

/* indicate the order (position) of elements in two sorted lists,
   as if they were concatenated.  */

void mergeorder(int *nn1, double *x1, int *nn2, double *x2, int *out) {
int i1,i2,i3,n1,n2,n3;
i1=i2=i3=0;
n1=*nn1;
n2=*nn2;
n3=n1+n2;;

while(i1<n1 && i2<n2) {
   if(x1[i1]<=x2[i2]) out[i3++]=1+i1++;
   else out[i3++]=1+n1+i2++;
   }
while(i1<n1) out[i3++]=1+i1++;
while(i2<n2) out[i3++]=1+n1+i2++;
}

/* comparison function for order_c */

int pdcmp(const void *x, const void *y) {
   double xx, yy;
   xx=**(double**)x;
   yy=**(double**)y;
   if(xx>yy) return 1;
   else if(xx<yy) return -1;
   else return 0;
   }

void order_c(double *in, int *n, int *out) {
   int i;
   double **p;
// create array of pointers
   p = Calloc(*n, double*);
   for(i=0;i<*n;i++) p[i]=in+i;
// sort indirect through pointers
   qsort(p,*n,sizeof(*p),pdcmp);
// convert pointers to positions
   for(i=0;i<*n;i++) out[i]=p[i]-in+1;
   Free(p);
   }

/******************  NOTE  *****************************

The following routines often calculate statistics like

  sum((x[i]-mean(x))^2)   

  I do this as sum(x^2)-mean(x)^2*length(x)

  which is mathematically equivalent, but not always computationally.
  If x values vary widely, especially of different sign, the results
  can round differently.   I don't believe that is a problem here, but
  am noting it just in case.

*********************************************************/

void sec_link_stat(int *n, double *trait, double *cistrait, 
         int *genotype, int *ng, double *lik0, double *lik1) {

int ntr, nct, i, j;
int *ll;
double *mx, *mc, *xc, *xx, *cc;
double  mt;
double Scc, Sct, Stt0, Stt1, Sctx;

nct=ntr=*n;
//allocation
ll = Calloc(*ng, int);
mx = Calloc(*ng, double);
mc = Calloc(*ng, double);
xc = Calloc(*ng, double);
xx = Calloc(*ng, double);
cc = Calloc(*ng, double);


//initialize
for(i = 0; i < *ng; i++){
	ll[i] = 0;
	mx[i] = mc[i] = xc[i] = xx[i] = cc[i] = 0.;
}
Scc=Sct=Stt0=Stt1=0.;

/* sum up the appropriate entries from each vector,
  counting them to calculate a mean */

for(i=0;i<ntr;i++) {
	j = genotype[i] - 1;
	ll[j]++;;
	mx[j] += trait[i];
	mc[j] += cistrait[i];
	xc[j] += trait[i]*cistrait[i];
	xx[j] += trait[i]*trait[i];
	cc[j] += cistrait[i]*cistrait[i];
	if (genotype[i] > *ng) assert("bad genotype");
   }

/* convert sums to means */
mt = 0.;
for(i = 0; i < *ng; i ++){
	mt = mt + mx[i];
	mx[i] = mx[i] / ll[i];
	mc[i] = mc[i] / ll[i];
	
}

/* generate the result statistic */
Stt1 = Stt0 = 0.;
Scc = Sct = 0.;
int L = 0;

for(i = 0; i < *ng; i ++){
	Stt1 = Stt1 + (xx[i] - (mx[i] * mx[i] * ll[i]));
	Stt0 = Stt0 + xx[i];
	Scc = Scc + (cc[i] - (mc[i] * mc[i] * ll[i]));
	Sct = Sct + (xc[i] - (mc[i] * mx[i] * ll[i]));
	L = L + ll[i];
}
mt = mt/L;
Stt0 = Stt0 - mt*mt*L;
Sctx = Sct*Sct/Scc;

*lik0 = Stt0-Sctx;
*lik1 = Stt1-Sctx;

Free(ll); Free(mx); Free(mc); Free(xc); Free(xx); Free(cc);
}

/* like sec.link.stat, but for an matrix of traits */
void sec_link_stat_x(int *n, int *kk, double *trait, double *cistrait, 
         int *genotype, int *ng, double *lik0, double* lik1) {

int ntr, nct, i, j, k, nk;
int *ll;
double *mx, *mc, *xc, *xx, *cc;
double mt, tr, ct;
double Scc, Sct, Stt0, Stt1, Sctx;

//allocation
ll = Calloc(*ng, int);
mx = Calloc(*ng, double);
mc = Calloc(*ng, double);
xc = Calloc(*ng, double);
xx = Calloc(*ng, double);
cc = Calloc(*ng, double);


nct = ntr = *n;
nk = *kk;

/* sum up the appropriate entries from each vector,
  counting them to calculate a mean */

for(k=0;k<nk;k++) {
	//initialize
	for(i = 0; i < *ng; i++){
		ll[i] = 0;
		mx[i] = mc[i] = xc[i] = xx[i] = cc[i] = 0.;
	}
	Scc=Sct=Stt0=Stt1=0.;
	
	for(i=0;i<ntr;i++) {
		j = genotype[i] - 1;
		tr = trait[k+nk*i];
		ct = cistrait[i];
		ll[j]++;
		mx[j] += tr;
		mc[j] += ct;
		xc[j] += tr*ct;
		xx[j] += tr*tr;
		cc[j] += ct*ct;
		if (genotype[i] > *ng) assert("bad genotype");
	}


	/* convert sums to means */
	mt = 0.;
	for(i = 0; i < *ng; i ++){
		mt = mt + mx[i];
		mx[i] = mx[i] / ll[i];
		mc[i] = mc[i] / ll[i];
	}


	/* generate the result statistic */
	Stt1 = Stt0 = 0.;
	Scc = Sct = 0.;
	int L = 0;

	for(i = 0; i < *ng; i ++){
		Stt1 = Stt1 + (xx[i] - (mx[i] * mx[i] * ll[i]));
		Stt0 = Stt0 + xx[i];
		Scc = Scc + (cc[i] - (mc[i] * mc[i] * ll[i]));
		Sct = Sct + (xc[i] - (mc[i] * mx[i] * ll[i]));
		L = L + ll[i];
	}
	mt = mt/L;
	Stt0 = Stt0 - mt*mt*L;
	Sctx = Sct*Sct/Scc;

	lik0[k]=Stt0-Sctx;
	lik1[k]=Stt1-Sctx;
	}
Free(ll); Free(mx); Free(mc); Free(xc); Free(xx); Free(cc);
}




void condi_indep_stat(int *n, double *trait, double *cistrait, 
         int *genotype, int *ng, double *stat) {

int ntr, nct, i, j;

/* variables to sum over combinations of the inputs */
double *mx, *mc, *xc, *xx, *cc;
double mxa, mca, xca, xxa, cca;
double cor, corn, cord, xry, llik0, llik1, m;
int *ll;

nct = ntr = *n;

//allocation
ll = Calloc(*ng, int);
mx = Calloc(*ng, double);
mc = Calloc(*ng, double);
xc = Calloc(*ng, double);
xx = Calloc(*ng, double);
cc = Calloc(*ng, double);


//initialize
for(i = 0; i < *ng; i++){
	ll[i] = 0;
	mx[i] = mc[i] = xc[i] = cc[i] = xx[i] = 0.;
}


/* sum up the appropriate entries from each vector,
  counting them to calculate a mean */

for(i=0;i<ntr;i++) {
	j = genotype[i]-1;
	ll[j]++;
	mx[j] += trait[i];
	mc[j] += cistrait[i];
	xc[j] += trait[i]*cistrait[i];
	xx[j] += trait[i]*trait[i];
	cc[j] += cistrait[i]*cistrait[i];
	if (genotype[i] > *ng) assert("bad genotype");
	}

// compute cor before dividing means

mxa = mca = xca = xxa = cca = 0.;
for(i = 0; i < *ng; i ++){
	mxa = mxa + mx[i];
	mca = mca + mc[i];
	xca = xca + xc[i];
	xxa = xxa + xx[i];
	cca = cca + cc[i];
}


/* calculate correlation coefficient */

corn = ntr*xca-mxa*mca;
cord = (ntr*xxa-mxa*mxa)*(ntr*cca-mca*mca);
cor = corn/sqrt(cord);

/* generate the intermediate results */
xry = (mxa-cor*mca)/ntr;
m = (xxa-2*cor*xca+cor*cor*cca)/ntr-xry*xry;
llik0 = ntr*log(m);

llik1 = 0.;
for(i = 0; i < *ng; i++) {
	xry = (mx[i] - cor*mc[i])/ll[i];
	m = (xx[i] - 2*cor*xc[i] + cor*cor*cc[i])/ll[i]-xry*xry;
	llik1 = llik1 + ll[i]*log(m);
}

*stat = llik0-llik1;
Free(ll); Free(mx); Free(mc); Free(xc); Free(xx); Free(cc);
}

void condi_indep_stat_x(int *n, int *kk, double *trait, double *cistrait, 
         int *genotype, int *ng, double *stat) {

int ntr, nct, i, j, k, nk;

/* variables to sum over combinations of the inputs */
double ct, tr;
double *mx, *mc, *xc, *xx, *cc;
double mxa, mca, xca, xxa, cca;
double cor, corn, cord, xry, llik0, llik1, m;
int *ll;

nct = ntr = *n;
nk = *kk;

//allocation
ll = Calloc(*ng, int);
mx = Calloc(*ng, double);
mc = Calloc(*ng, double);
xc = Calloc(*ng, double);
xx = Calloc(*ng, double);
cc = Calloc(*ng, double);


/* sum up the appropriate entries from each vector,
  counting them to calculate a mean */

for(k=0;k<nk;k++) {
   //initialize
	for(i = 0; i < *ng; i++){
		ll[i] = 0;
		mx[i] = mc[i] = xc[i] = cc[i] = xx[i] = 0.;
	}

	for(i = 0; i < ntr; i++) {
		tr = trait[k+nk*i];
		ct = cistrait[i];
		j = genotype[i] - 1;
		ll[j]++;
		mx[j] += tr;
		mc[j] += ct;
		xc[j] += tr*ct;
		xx[j] += tr*tr;
		cc[j] += ct*ct;
		if (genotype[i] > *ng) assert("bad genotype");
	}

	// compute cor before dividing means

	mxa = mca = xca = xxa = cca = 0.;
	for(i = 0; i < *ng; i ++){
		mxa = mxa + mx[i];
		mca = mca + mc[i];
		xca = xca + xc[i];
		xxa = xxa + xx[i];
		cca = cca + cc[i];
	}

	/* calculate correlation coefficient */

	corn=ntr*xca-mxa*mca;
	cord=(ntr*xxa-mxa*mxa)*(ntr*cca-mca*mca);
	cor=corn/sqrt(cord);

	/* generate the intermediate results */

	xry = (mxa-cor*mca)/ntr;
	m = (xxa - 2*cor*xca+cor*cor*cca)/ntr-xry*xry;
	llik0 = ntr*log(m);

	llik1 = 0.;
	for(i=0; i < *ng; i++) {
		xry = (mx[i] - cor*mc[i])/ll[i];
		m = (xx[i] - 2*cor*xc[i] + cor*cor*cc[i])/ll[i]-xry*xry;
		llik1 = llik1 + ll[i]*log(m);
	}
	
	stat[k] = llik0 - llik1;
	}
Free(ll); Free(mx); Free(mc); Free(xc); Free(xx); Free(cc);
}


void condi_indep_stat_rx(int *n, int *kk, double *trait, double *cistrait, 
         int *genotype, int *ng, double *stat) {

int ntr, nct, i, j, k, nk;

/* variables to sum over combinations of the inputs */
double ct, tr;
double *mx, *mc, *xc, *xx, *cc;
double mxa, mca, xca, xxa, cca;
double cor, corn, cord, xry, llik0, llik1, m;
int *ll;

nct = ntr = *n;
nk = *kk;

//allocation
ll = Calloc(*ng, int);
mx = Calloc(*ng, double);
mc = Calloc(*ng, double);
xc = Calloc(*ng, double);
xx = Calloc(*ng, double);
cc = Calloc(*ng, double);


/* sum up the appropriate entries from each vector,
  counting them to calculate a mean */

for(k=0;k<nk;k++) {
   //initialize
	for(i = 0; i < *ng; i++){
		ll[i] = 0;
		mx[i] = mc[i] = xc[i] = cc[i] = xx[i] = 0.;
	}

	for(i = 0; i < ntr; i++) {
		tr = trait[i];
		ct = cistrait[k+nk*i];
		j = genotype[i] - 1;
		ll[j]++;
		mx[j] += tr;
		mc[j] += ct;
		xc[j] += tr*ct;
		xx[j] += tr*tr;
		cc[j] += ct*ct;
		if (genotype[i] > *ng) assert("bad genotype");
	}

	// compute cor before dividing means

	mxa = mca = xca = xxa = cca = 0.;
	for(i = 0; i < *ng; i ++){
		mxa = mxa + mx[i];
		mca = mca + mc[i];
		xca = xca + xc[i];
		xxa = xxa + xx[i];
		cca = cca + cc[i];
	}

	/* calculate correlation coefficient */

	corn=ntr*xca-mxa*mca;
	cord=(ntr*xxa-mxa*mxa)*(ntr*cca-mca*mca);
	cor=corn/sqrt(cord);

	/* generate the intermediate results */

	xry = (mxa-cor*mca)/ntr;
	m = (xxa - 2*cor*xca+cor*cor*cca)/ntr-xry*xry;
	llik0 = ntr*log(m);

	llik1 = 0.;
	for(i=0; i < *ng; i++) {
		xry = (mx[i] - cor*mc[i])/ll[i];
		m = (xx[i] - 2*cor*xc[i] + cor*cor*cc[i])/ll[i]-xry*xry;
		llik1 = llik1 + ll[i]*log(m);
	}
	
	stat[k] = llik0 - llik1;
	}
Free(ll); Free(mx); Free(mc); Free(xc); Free(xx); Free(cc);
}



void link_stat(int* n, double *trait, int *genotype, int* ng, double *lik0, double* lik1 ) {
//ng is the number of group to compare, must be user input
int ntr, i, j;
int* ll;
double mt;
double* mx, *xs;
double Stt0, Stt1;

//allocation
ll = Calloc(*ng, int);
mx = Calloc(*ng, double);
xs = Calloc(*ng, double);

//initialize
for(i = 0; i < *ng; i++){
	ll[i] = 0;
	mx[i] = 0.;
	xs[i] = 0.;
}

ntr=*n;
/* sum up the appropriate entries from each vector,
  counting them to calculate a mean */

for(i=0;i<ntr;i++) {
   j= genotype[i] - 1;
   ll[j]++;
   mx[j] += trait[i];
   xs[j] += trait[i]*trait[i];
   if (j > (*ng-1)) assert("bad genotype");
   }

/* convert sums to means */
mt = 0.;
for(i = 0; i < *ng; i ++){
	mt = mt + mx[i];
	mx[i] = mx[i] / ll[i];
}

/* generate the result statistic */
Stt1 = 0.;
Stt0 = 0.;
int L = 0;

for(i = 0; i < *ng; i ++){
	Stt1 = Stt1 + (xs[i] - (mx[i] * mx[i] * ll[i]));
	Stt0 = Stt0 + xs[i];
	L = L + ll[i];

}

mt = mt/L;
*lik0 = Stt0 - mt*mt*L;
*lik1 = Stt1;
Free(ll); Free(mx); Free(xs);
}




void link_stat_xx(int *n, int *kk, double *trait, 
     int *mm, int *genotype, int *ng, double *lik0, double * lik1) {

int ntr, i, j, k, m, nk, nm;
int *ll;
double mt, tr;
double *mx, *xs;
double Stt0, Stt1;

ntr=*n;
nk=*kk;
nm=*mm;

//allocation
ll = Calloc(*ng, int);
mx = Calloc(*ng, double);
xs = Calloc(*ng, double);


/* sum up the appropriate entries from each vector,
  counting them to calculate a mean */

for(m=0;m<nm;m++) for(k=0;k<nk;k++) {
	for(i = 0; i < *ng; i++){
		ll[i] = 0;
		mx[i] = 0.;
		xs[i] = 0.;
	}

   for(i = 0; i < ntr; i++) {
      j = genotype[m + nm * i] - 1;
      ll[j]++;
      tr = trait[k + nk * i];
      mx[j] += tr;
      xs[j] += tr * tr;
	  if (j > (*ng-1)) assert("bad genotype");
      }

	/* convert sums to means */
	mt = 0.;
	for(i = 0; i < *ng; i ++){
		mt = mt + mx[i];
		mx[i] = mx[i] / ll[i];
		
	}

	/* generate the result statistic */
	Stt1 = 0.;
	Stt0 = 0.;
	int L = 0;

	for(i = 0; i < *ng; i ++){
		Stt1 = Stt1 + (xs[i] - (mx[i] * mx[i] * ll[i]));
		Stt0 = Stt0 + xs[i];
		L = L + ll[i];
	}
	
	mt  = mt/L;
	lik0[k+nk*m] = Stt0 - mt*mt*L;
	lik1[k+nk*m] = Stt1;
	}
Free(ll); Free(mx); Free(xs);
}
