#include "AdCanny.h"
#include"math.h"
#include<afx.h>

AdCanny::AdCanny(void)
{
}

/*
this algorithm is changged canny,after compute the grad(x,y,z) then use the  directions(8*2 0��) and the grad to filter the image  
1��guass filter 
2��compute the x,y,z,diffrerention store in three images
3��compute the tan=deta y /detax;include 8 directions 0,45 90 135 180 225 270 315 360
4��vlue+= vr *grad + vl *grad;

input :double sigma, double **pdKernel, int *pnWidowSize
output:


*/
//  һά��˹�ֲ�����������ƽ�����������ɵĸ�˹�˲�ϵ��
void CreatGauss(double sigma, double **pdKernel, int *pnWidowSize)
{
 
 LONG i;
 
 //�������ĵ�
 int nCenter;
 
 //������һ�㵽���ĵ����
 double dDis;
 
 //�м����
 double dValue;
 double dSum;
 dSum = 0;
 
 // [-3*sigma,3*sigma] �������ݣ��Ḳ�Ǿ��󲿷��˲�ϵ��
 *pnWidowSize = 1+ 2*ceil(3*sigma);
 
 nCenter = (*pnWidowSize)/2;
 
 *pdKernel = new double[*pnWidowSize];
 
 //���ɸ�˹����
 for(i=0;i<(*pnWidowSize);i++)
 {
  dDis = double(i - nCenter);
  dValue = exp(-(1/2)*dDis*dDis/(sigma*sigma))/(sqrt(2*3.1415926)*sigma);
  (*pdKernel)[i] = dValue;
  dSum+=dValue;
 
 }
 //��һ��
 for(i=0;i<(*pnWidowSize);i++)
 {
  (*pdKernel)[i]/=dSum;
 }
 
}
 
//�ø�˹�˲���ƽ��ԭͼ��
void GaussianSmooth(SIZE sz, LPBYTE pGray, LPBYTE pResult, double sigma)
{
 LONG x, y;
 LONG i;
 
 //��˹�˲�������
 int nWindowSize;
 
 //���ڳ���
 int nLen;
 
 //һά��˹�˲���
 double *pdKernel;
 
 //��˹ϵ����ͼ�����ݵĵ��
 double dDotMul;
 
 //�˲�ϵ���ܺ� 
 double dWeightSum;
 
 double *pdTemp;
 pdTemp = new double[sz.cx*sz.cy];
 
 //����һά��˹����
 CreatGauss(sigma, &pdKernel, &nWindowSize);
 
 nLen = nWindowSize/2;
 
 //x�����˲�
 for(y=0;y<sz.cy;y++)
 {
  for(x=0;x<sz.cx;x++)
  {
   dDotMul = 0;
   dWeightSum = 0;
   for(i=(-nLen);i<=nLen;i++)
   {
    //�ж��Ƿ���ͼ���ڲ�
    if((i+x)>=0 && (i+x)<sz.cx)
    {
     dDotMul+=(double)pGray[y*sz.cx+(i+x)] * pdKernel[nLen+i];
     dWeightSum += pdKernel[nLen+i];
    }
   }
   pdTemp[y*sz.cx+x] = dDotMul/dWeightSum;
  }
 }
 
 //y�����˲�
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
 
// ������,���ݶ�
void Grad(SIZE sz, LPBYTE pGray,int *pGradX, int *pGradY, int *pMag)
{
 LONG y,x;
 
 //x����ķ�����
 for(y=1;y<sz.cy-1;y++)
 {
  for(x=1;x<sz.cx-1;x++)
  {
   pGradX[y*sz.cx +x] = (int)( pGray[y*sz.cx+x+1]-pGray[y*sz.cx+ x-1]  );
  }
 }
 
 //y��������
 for(x=1;x<sz.cx-1;x++)
 {
  for(y=1;y<sz.cy-1;y++)
  {
   pGradY[y*sz.cx +x] = (int)(pGray[(y+1)*sz.cx +x] - pGray[(y-1)*sz.cx +x]);
  }
 }
 
 //���ݶ�
 
 //�м����
 double dSqt1;
 double dSqt2;
 double dSqt3;
 
 for(y=0; y<sz.cy; y++)
 {
  for(x=0; x<sz.cx; x++)
  {
   //���׷������ݶ�
   dSqt1 = pGradX[y*sz.cx + x]*pGradX[y*sz.cx + x];
   dSqt2 = pGradY[y*sz.cx + x]*pGradY[y*sz.cx + x];
   pMag[y*sz.cx+x] = (int)(sqrt(dSqt1+dSqt2)+0.5);
  }
 }
}
AdCanny::~AdCanny(void)
{
}
