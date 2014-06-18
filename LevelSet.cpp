#include "LevelSet.h"
#include"ImageF.h"


LevelSet::LevelSet(void)
{
}


LevelSet::~LevelSet(void)
{
}
ImageF drlse_edge(ImageF phi_0,ImageF g,float lambda,float mu,float alfa,float epsilon,int timestep, int iter,string sada,char * potentialFunction)
{
	ImageF phi=phi_0;

	for(int i=0;i<iter;i++)
	{
		phi=NeumannBoundCond(phi);
		ImageF *vx=new ImageF();
		ImageF *vy=new ImageF();
		(*vy)*(vx);
		ImageF *phi_x=new ImageF();
		ImageF *phi_y=new ImageF();
		*vx=gradientx(g);
		*vy=gradienty(g);
		*phi_x=gradientx(g);
		*phi_y=gradienty(g);
		ImageF *s=new ImageF();
		*s=ImageFSqrt( *vx, *vy);
		float smallNumber=1e-10;
		ImageF *Nx,*Ny;
		*Nx=*phi_x/(*s+smallNumber);
		*Ny=*phi_y/(*s+smallNumber);
		ImageF * curvature=new ImageF();
		*curvature=div(*Nx,*Ny);
		ImageF distRegTerm;
		if (strcmp(potentialFunction,"single-well"))
			/*
			compute distance regularization term in equation (13) 
			with the single-well potential p1.
			*/
			distRegTerm= 4*del2(phi)-*curvature;
		//printf("");

		else if (strcmp(potentialFunction,"double-well"))
		{
			distRegTerm=distReg_p2(phi);  // compute the distance regularization term in eqaution (13) with the double-well potential p2.

			printf("asda");

		}
		else printf("EEROR");

		//float eplsion=0.1;
		ImageF diracPhi=Dirac(phi,epsilon);
		ImageF areaTerm=diracPhi*g; 
		ImageF  *edgeTerm=new ImageF();
		*edgeTerm=diracPhi*(*vx**Nx+*vy**Ny) + diracPhi*g*(*curvature);
		//vx*vy;
		 phi=phi + timestep*(mu*distRegTerm + lambda**edgeTerm + alfa*areaTerm);
		return phi_0; 
	}
}
ImageF Dirac(ImageF phi,float epsoilon)
{

}
ImageF gradientx(ImageF g )
{
	


}
ImageF gradienty(ImageF g )
{
	


}

ImageF NeumannBoundCond(ImageF img)
{

}

ImageF ImageFSqrt(ImageF x,ImageF y)
{


}
ImageF operator *(int p,ImageF x)
{
}
ImageF operator * (ImageF x,ImageF y)
{

}
ImageF operator/(ImageF &x,ImageF& y)
{
	/*
	实现数组除法，对应位置的元素相除
	*/
}
ImageF operator+(ImageF &x,float s)
{
	/*
	数组+标量s：每个元素都加s
	*/
}
ImageF div(ImageF x,ImageF y)
{
	return (x+y);
}
ImageF operator+(ImageF x,ImageF y)
{
}

ImageF operator -(ImageF x,ImageF y)
{
}
	
ImageF del2(ImageF phi)
{
	return phi;
}
ImageF regFunction(ImageF s)
{
	//if((s>=0) && (s<=1))
		return s;

}
ImageF distReg_p2(ImageF phi)
{
	ImageF *phi_x=new ImageF();
	ImageF *phi_y=new ImageF();
	*phi_x=gradientx(phi);
	*phi_y=gradienty(phi);
	//ImageF s=ImageFSqrt(((*phi_x)*(*phi_x)) + ((*phi_y)*(*phi_y)));
	ImageF  s=ImageFSqrt(*phi_x,*phi_y);
	ImageF a=regFunction(s);
	ImageF b=regFunction(s);//此处的函数需该
	ImageF ps=

//ImageF b=(s>1);
}
