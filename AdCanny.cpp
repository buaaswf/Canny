#include "AdCanny.h"
#include"math.h"
#include<afx.h>

AdCanny::AdCanny(void)
{
}

/*
this algorithm is changged canny,after compute the grad(x,y,z) then use the  directions(8*2 0π) and the grad to filter the image  
1、guass filter 1d
2、compute the x,y,z,diffrerention store in three images
3、compute the tan=deta y /detax;include 8 directions 0,45 90 135 180 225 270 315 360
4、vlue+= vr *grad + vl *grad;

input :double sigma, double **pdKernel, int *pnWidowSize
output:


*/


/*
2.高斯滤波器的kernel是可分离的(separable)，也就是说，可以将2D的高斯kernel分解为两个1D的kernel，
先沿x方向对图像进行1D高斯kernel的卷积，然后沿y方向对图像进行1D的高斯kernel卷积，
最后的结果和使用一个2D高斯kernel对图像卷积效果是一样的。这样一来，针对每个像素，滤波器的算法复杂度降为O(r)。
 一维高斯分布函数，用于平滑函数中生成的高斯滤波系数 ,也可以用而为高斯函数，甚至三维高斯函数，效果有待比较

 */ 


/*intput:sigma正态分布的方差，pdKernel 存放高斯系数的数组，

*/
void CreatGauss(double sigma, double **pdKernel, int *pnWidowSize)
{
 
 LONG i;
 
 //数组中心点
 int nCenter;
 
 //数组中一点到中心点距离
 double dDis;
 
 //中间变量
 double dValue;
 double dSum;
 dSum = 0;
 
 // [-3*sigma,3*sigma] 以内数据，会覆盖绝大部分滤波系数
 *pnWidowSize = 1+ 2*ceil(3*sigma);
 
 nCenter = (*pnWidowSize)/2;
 
 *pdKernel = new double[*pnWidowSize];
 
 //生成高斯数据
 for(i=0;i<(*pnWidowSize);i++)
 {
  dDis = double(i - nCenter);
  dValue = exp(-(1/2)*dDis*dDis/(sigma*sigma))/(sqrt(2*3.1415926)*sigma);
  (*pdKernel)[i] = dValue;
  dSum+=dValue;
 
 }
 //归一化
 for(i=0;i<(*pnWidowSize);i++)
 {
  (*pdKernel)[i]/=dSum;
 }
 
}
 
//用高斯滤波器平滑原图像
void GaussianSmooth(SIZE sz, LPBYTE pGray, LPBYTE pResult, double sigma)
{
 LONG x, y;
 LONG i;
 
 //高斯滤波器长度
 int nWindowSize;
 
 //窗口长度
 int nLen;
 
 //一维高斯滤波器
 double *pdKernel;
 
 //高斯系数与图像数据的点乘
 double dDotMul;
 
 //滤波系数总和 
 double dWeightSum;
 
 double *pdTemp;
 pdTemp = new double[sz.cx*sz.cy];
 
 //产生一维高斯数据
 CreatGauss(sigma, &pdKernel, &nWindowSize);
 
 nLen = nWindowSize/2;
 
 //x方向滤波
 for(y=0;y<sz.cy;y++)
 {
  for(x=0;x<sz.cx;x++)
  {
   dDotMul = 0;
   dWeightSum = 0;
   for(i=(-nLen);i<=nLen;i++)
   {
    //判断是否在图像内部
    if((i+x)>=0 && (i+x)<sz.cx)
    {
     dDotMul+=(double)pGray[y*sz.cx+(i+x)] * pdKernel[nLen+i];
     dWeightSum += pdKernel[nLen+i];
    }
   }
   pdTemp[y*sz.cx+x] = dDotMul/dWeightSum;
  }
 }
 
 //y方向滤波
 for(x=0; x<sz.cx;x++)
 {
  for(y=0; y<sz.cy; y++)
  {
   dDotMul = 0;
   dWeightSum = 0;
   for(i=(-nLen);i<=nLen;i++)
   {
    if((i+y)>=0 && (i+y)< sz.cy)
    {
     dDotMul += (double)pdTemp[(y+i)*sz.cx+x]*pdKernel[nLen+i];
     dWeightSum += pdKernel[nLen+i];
    }
   }
   pResult[y*sz.cx+x] = (unsigned char)dDotMul/dWeightSum;
  }
 }
 
 delete []pdKernel;
 pdKernel = NULL;
 
 delete []pdTemp;
 pdTemp = NULL;
 
}
 
// 方向导数,求梯度
/*
input:SIZE sz, LPBYTE pGray,int *pGradX, int *pGradY, int pGradz,int *pMag


*/
void Grad(SIZE sz, LPBYTE pGray,int *pGradX, int *pGradY,int *pGradZ, int *pMag)
{
 LONG y,x,z;
 
 //x方向的方向导数
 for(y=1;y<sz.cy-1;y++)
 {
  for(x=1;x<sz.cx-1;x++)
  {
   pGradX[y*sz.cx +x] = (int)( pGray[y*sz.cx+x+1]-pGray[y*sz.cx+ x-1]  );
  }
 }
 
 //y方向方向导数
 for(x=1;x<sz.cx-1;x++)
 {
  for(y=1;y<sz.cy-1;y++)
  {
   pGradY[y*sz.cx +x] = (int)(pGray[(y+1)*sz.cx +x] - pGray[(y-1)*sz.cx +x]);
  }
 }
 
 //求梯度
 //z方向 方向导数


 //中间变量
 double dSqt1;
 double dSqt2;
 double dSqt3;
 
 for(y=0; y<sz.cy; y++)
 {
  for(x=0; x<sz.cx; x++)
  {
   //二阶范数求梯度
   dSqt1 = pGradX[y*sz.cx + x]*pGradX[y*sz.cx + x];
   dSqt2 = pGradY[y*sz.cx + x]*pGradY[y*sz.cx + x];
   pMag[y*sz.cx+x] = (int)(sqrt(dSqt1+dSqt2)+0.5);
  }
 }
}
void Direction(SIZE sz, LPBYTE pGray,int *pGradX, int *pGradY,int *pMag,int *tan )
{
	float *tan;
	int i;
	for(i=1;i<sz.cx-1;i++)
	{
		tan[i]=pGradY[i]/pGradX[i];
	}
	switch (tan[i])
	{
	case 1
		tan[i]=1;break;
	case 2
		tan[i]=2;break;
	case 3
		tan[i]=3;break;
	default
	}
}
/*
input:
function:value+=
*/



void Value(SIZE sz,int *tan,int *pMag,LPBYTE* pGray)
{
	LPBYTE value=pGray;
	
	int i=1;
	for(i=0;i<sz.cx-1;i++)
	{
		value[i]+=
	}
	
}

AdCanny::~AdCanny(void)
{
}
