
namespace DockTK { namespace nr {

template<typename T_num>
void eigsrt(T_num d[], T_num **v, int n,bool ascend)
{
	int k,j,i;
	T_num p;

	for (i=1;i<n;i++) {
		p=d[k=i];
		for (j=i+1;j<=n;j++)
		  if(ascend) {
		    if (d[j] <= p) p=d[k=j];
		  }
		  else {
		    if (d[j] >= p) p=d[k=j];
		  }
		if (k != i) {
			d[k]=d[i];
			d[i]=p;
			for (j=1;j<=n;j++) {
				p=v[j][i];
				v[j][i]=v[j][k];
				v[j][k]=p;
			}
		}
	}
}

} } // namespace DockTK::nr


// instantiate

namespace DockTK { namespace nr {

template void eigsrt(float d[], float **v, int n, bool ascend);
template void eigsrt(double d[], double **v, int n, bool ascend);

} } // namespace DockTK::nr

/* (C) Copr. 1986-92 Numerical Recipes Software =$jX16#S'j#(+1. */
