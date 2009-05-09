// ary3d.h ///////////////////////////////////////////////////////////
// 2D/3D array classes 
//  double3d , float3d , int3d , short3d , char3d
//  double2d , float2d , int2d , short2d , char2d
//
//  Usage:
//   int x=5,y=4,z=6,i,j,k;                int x=6,y=5,i,j;
//
//   // declaration                        // declaration
//   double3d a(x,y,z);                    float2d a(x,y);
//
//   // assignment                         // assignment
//   for(i=0;i<x;i++)                      for(i=0;i<x;i++) 
//     for(j=0;j<y;j++)                      for(j=0;j<y;j++)
//       for(k=0;i<z;k++)                      a(i,j) = i*2+j;
//         a(i,j,i) = i*j+k;
//
//   // reference                          // reference
//   for(i=0;i<x;i++)                      for(i=0;i<x;i++)
//     for(j=0;j<y;j++)                      for(j=0;j<y;j++)
//       for(k=0;k<z;k++)                      cout << a(i,j) << endl;
//         cout << a(i,j,k) << endl;
//
// NOTE1:
//   It is recommended to precompile this header file. Precompiled version
//  is faster and more reliable.
//     $  makecint -mk Makefile -dl ary3d.dll -H ary3d.h
//     $  make -f Makefile
//
// NOTE2:
//   Only simple run-time array index check is done. Run-time array index
//  check can be turned off by undefining INDEXCHECK macro.
//
// NOTE3:
//   a[i][j][k] notation is available if INDEXOPR macro is defined.
//  This notation is less efficient compared to a(i,j,k). Not recommended.
//
/////////////////////////////////////////////////////////////////////
#include <iostream>

//#define INDEXOPR
#define INDEXCHECK
/////////////////////////////////////////////////////////////////////
// 2 dimentional array class template
/////////////////////////////////////////////////////////////////////
template<class T> class array2d {
	public:
		// constructors
		array2d(int jx, int kx) 
		{ 
			maxj = jx; 
			maxk = kx;
			isnew = 1;
			p = new T[maxj*maxk];
		}
		array2d(array2d& obj) { copy(obj); }

#ifdef INDEXOPR
		array2d() { p=NULL; }
		//array2d(T* pin,int jx,int kx) {init(pin,jx,kx);}
		void init(T* pin,int jx,int kx) 
		{ 
			maxj = jx;  
			maxk = kx;
			isnew=0;
			p = pin;
		} 
#endif

		// destructor
		~array2d() {if(isnew&&p) delete[] p;}

		// operator overloading
		array2d& operator=(array2d& obj) 
		{
			if(p==obj.p) return(*this);
			if(isnew&&p) delete[] p;
			copy(obj);
			return(*this);
			}
			T& operator()(int j,int k) {
#ifdef INDEXCHECK
			if(j < 0 || maxj <= j|| k < 0|| maxk <= k) 
			{
				cerr << "Bad index ("<<j<<","<<k<<") > ("<<maxj<<","<<maxk<<")"<< endl;
				return(*p);
			}
#endif
			return( *(p + maxk*j + k) );
		}
#ifdef INDEXOPR
		T* operator[](int j) {return(p+maxk*j);}
#endif
		friend bool operator==(array2d<T>& a,array2d<T>& b);
		friend bool operator!=(array2d<T>& a,array2d<T>& b);
		friend ostream& operator<<(ostream& ost,array2d<T>& a);

	private:
		T* p;
		int maxj,maxk;
		int isnew;

		// utility function
		void copy(array2d& obj) 
		{
			maxj=obj.maxj; maxk=obj.maxk;
			if(isnew) {
				isnew = 1;
				p = new T[maxj*maxk];
				memcpy(p,obj.p,sizeof(T)*maxj*maxk);
			}
			else {
				isnew = 0;
				p = obj.p;
			}
		}
};

template<class T>
bool operator==(array2d<T>& a,array2d<T>& b) 
{
	if(a.maxj != b.maxj || a.maxk != b.maxk) 
		return false;

	int i, max = (a.maxj * a.maxk);
	for( i = 0; i < max; i++ ) 
	{
		if(a.p[i] != b.p[i]) 
			return false;
	}
	return true;
}

template<class T>
bool operator!=(array2d<T>& a,array2d<T>& b) 
{
  return(!(a==b));
}

template<class T>
ostream& operator<<(ostream& ost,array2d<T>& a) 
{
	int j,k;
	ost << '(' ;
	for(j=0;j<a.maxj;j++) 
	{
		if(j) ost << ',';
		ost << '(' ;
		for(k=0;k<a.maxk;k++) 
		{
			if(k) ost << ',';
				ost << a(j,k);
		}
		ost << ')' ;
	}
	ost << ')' ;
	return(ost);
}

/////////////////////////////////////////////////////////////////////
// 3 dimentional array class template
/////////////////////////////////////////////////////////////////////
template<class T> class array3d {
public:
	// constructor
	array3d(int ix,int jx,int kx) 
	{ 
		maxi = ix;
		maxj = jx;
		maxk = kx;
		isnew = 1;
		p = new T[maxi*maxj*maxk];
	}
	array3d(array3d& obj) {copy(obj);}
#ifdef INDEXOPR
	array3d() { p = NULL; }
	//array3d(T* pin,int ix,int jx,int kx) {init(pin,ix,jx,kx);}
	void init(T* pin,int ix,int jx,int kx) 
	{ 
		maxi = ix; 
		maxj = jx; 
		maxk = kx;
		isnew = 0;
		p = pin;
	}
#endif

	// destructor
	~array3d() {if(isnew&&p) delete[] p;}

	// operator overloading
	array3d& operator=(array3d& obj) 
	{
		if(p==obj.p) 
			return(*this);
		if(isnew&&p) 
			delete[] p;
		copy(obj);
		return(*this);
	}
	T& operator()(int i,int j,int k) 
	{
#ifdef INDEXCHECK
		if(i<0||maxi<=i||j<0||maxj<=j||k<0||maxk<=k) 
		{
			cerr << "Bad index ("<<i<<","<<j<<","<<k<<") > (";
			cerr <<maxi<<","<<maxj<<","<<maxk<<")"<< endl;
			return(*p);
		}
#endif
		return( *(p + ((maxj*maxk)*i + maxk*j + k)) );
	}
#ifdef INDEXOPR
	array2d<T>& operator[](int i) 
	{
		subary.init(p+maxj*maxk*i,maxj,maxk);
		return(subary);
	}
#endif
	friend bool operator==(array3d<T>& a,array3d<T>& b);
	friend bool operator!=(array3d<T>& a,array3d<T>& b);
	friend ostream& operator<<(ostream& ost,array3d<T>& a);

private:
	T* p;
	int maxi, maxj, maxk;
	int isnew;
#ifdef INDEXOPR
	array2d<T> subary;
#endif

	// utility function
	void copy(array3d& obj) 
	{
		maxi = obj.maxi; 
		maxj = obj.maxj; 
		maxk = obj.maxk;
		if(isnew) 
		{
			isnew = 1;
			p = new T[maxi*maxj*maxk];
			memcpy(p,obj.p, sizeof(T)*maxi*maxj*maxk);
		}
		else 
		{
			isnew = 0;
			p = obj.p;
		}
	}
};

template<class T>
bool operator==(array3d<T>& a,array3d<T>& b) 
{
	if(a.maxi!=b.maxi || a.maxj!=b.maxj || a.maxk!=b.maxk) 
		return false;
	int i, max = (a.maxi * a.maxj * a.maxk);
	for(i = 0 ; i < max ; i++) 
	{
		if(a.p[i] != b.p[i]) 
			return false;
	}
	return true;
}

template<class T>
bool operator!=(array3d<T>& a,array3d<T>& b) 
{
	return(!(a==b));
}

template<class T>
ostream& operator<<(ostream& ost,array3d<T>& a) 
{
	int i,j,k;
	ost << '(' ;
	for(i = 0 ; i < a.maxi ; i++) 
	{
		if(i) 
			ost << ',';
		ost << '(' ;
		for(j = 0 ; j < a.maxj ; j++) 
		{
			if(j) ost << ',';
			ost << '(' ;
			for(k = 0 ; k < a.maxk ; k++) 
			{
				if(k) 
					ost << ',';
				ost << a(i,j,k);
			}
			ost << ')' ;
		}
		ost << ')' ;
	}
	ost << ')' ;
	return(ost);
}

/////////////////////////////////////////////////////////////////////
// instantiation of class template as typedef
/////////////////////////////////////////////////////////////////////
typedef array3d<double> double3d;
typedef array2d<double> double2d;
typedef array3d<float>  float3d;
typedef array2d<float>  float2d;
typedef array3d<int>    int3d;
typedef array2d<int>    int2d;
typedef array3d<short>  short3d;
typedef array2d<short>  short2d;
typedef array3d<char>   char3d;

