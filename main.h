
#pragma once
#include <afxwin.h>  // necessary for MFC to work properly
#include <math.h> 
#include "../../src/blepo.h"

using namespace blepo;
using namespace blepo_ex;
//using namespace blepo_ex;

int initial_filename_recognize( const CString path_input, const CString filename_input)
{
	if(   fopen(path_input+filename_input,"r")==0   )
	  { return 1; }
	
	else return 0;
				
}


//convolution_x   with kernel
void separable_convolution_x( const ImgFloat &img_in,  const int kernel_width, const float gaussian_kernel[], ImgFloat* img_out ) 
{

	(*img_out).Reset( img_in.Width(),img_in.Height() );
	Set(img_out, 0.0);

	
	for( int y=0; y< img_in.Height(); y++ ) 
	{
		for( int x=(kernel_width-1)/2; x< img_in.Width()-(kernel_width-1)/2; x++ ) 
		{
			float tmp_value=0;
			for(int ctn=0; ctn<kernel_width ; ctn++)
			{

				if(	  ( (x+ (kernel_width-1)/2 - ctn) < 0 ) 
					||( (x+ (kernel_width-1)/2 - ctn) >= img_in.Width() )   
				  )

				{
					tmp_value+=gaussian_kernel[ctn] * 0.0f;				
				}

				else 
				{
					tmp_value+=gaussian_kernel[ctn] * img_in(x+ (kernel_width-1)/2 - ctn,y);
				}

			}//for ctn

			(*img_out)(x,y)=tmp_value;
			

		}//for x
	}//for y

	//return 0;

}


//convolution_y   with kernel
void separable_convolution_y( const ImgFloat &img_in,  const int kernel_width, const float gaussian_kernel[], ImgFloat* img_out ) 
{
	(*img_out).Reset( img_in.Width() , img_in.Height() );
	Set(img_out, 0.0);

	for( int x=0; x< img_in.Width(); x++) 
	{
		for(  int y=(kernel_width-1)/2; y< img_in.Height()-(kernel_width-1)/2; y++ ) 
		{
			float tmp_value=0;

			for(int ctn=0; ctn<kernel_width ; ctn++)
			{

				if(	  ( (y+ (kernel_width-1)/2 - ctn) < 0 ) 
					||( (y+ (kernel_width-1)/2 - ctn) >= img_in.Height() )   
				  )

				{
					tmp_value+=gaussian_kernel[ctn] * 0.0f;				
				}

				else 
				{
					tmp_value+=gaussian_kernel[ctn] * img_in(x,y + (kernel_width-1)/2 - ctn );
				}

			}//for ctn
			
			(*img_out)(x,y)=tmp_value;
			

		}//for x
	}//for y

	//return 0;

}



//convolution x and y  with kernel
int separable_convolution_2d( const ImgGray &img_in_g,  const int kernel_width,const float gaussian_kernel[] ,const float gaussian_deri_kernel[], ImgFloat* img_out, const int mode=0 ) 
{
	ImgFloat img_in_f;
	ImgFloat img_mid;
	

	Convert(img_in_g, &img_in_f);


	if(mode == 0) // x  y'
	{
		separable_convolution_y( img_in_f,  kernel_width, gaussian_deri_kernel, &img_mid );
		separable_convolution_x( img_mid,  kernel_width, gaussian_kernel,  img_out);
		
		return 0;
		
	}

	else if(mode == 1)// y x'
	{
		separable_convolution_x( img_in_f,  kernel_width, gaussian_deri_kernel, &img_mid );
		separable_convolution_y( img_mid,  kernel_width, gaussian_kernel,  img_out);
		
		return 0;
	}

	else
	{
		return 1;
	}
		
}


//make magnitude and anlge figure
int grad_magnitude_angle( const ImgFloat &img_in_x,  const ImgFloat &img_in_y,  ImgFloat* img_out_magnitude,  ImgFloat* img_out_angle ) 
{

	(*img_out_magnitude).Reset( img_in_x.Width(),img_in_x.Height() );
	Set(img_out_magnitude, 0.0);

	(*img_out_angle).Reset( img_in_x.Width(),img_in_x.Height() );
	Set(img_out_angle, 0.0);


	for( int y=0; y< img_in_x.Height(); y++ ) 
	{
		for( int x=0; x< img_in_x.Width(); x++ ) 
		{
			//(*img_out_magnitude)(x,y)=sqrt( abs( img_in_x(x,y) )*abs( img_in_x(x,y) ) +abs( img_in_y(x,y) )*abs( img_in_y(x,y) ) );
			(*img_out_magnitude)(x,y)=max(abs( img_in_x(x,y) ), abs( img_in_y(x,y)) );
			(*img_out_angle)(x,y)=atan2(img_in_y(x,y),img_in_x(x,y) );
			

		}//for x
	}//for y

	return 0;

}


//non max suppression 
int nms_filter( const ImgFloat &img_in_magnitude, const ImgFloat &img_in_angle, ImgFloat* img_out) 
{
	float point_max_value=0;
	(*img_out).Reset( img_in_magnitude.Width() , img_in_magnitude.Height() );
	Set(img_out, 0.0);
	//*img_out=img_in_magnitude;

	//+-pi8
	const float _pi8=3.14f/8;
	const float _npi8=-3.14f/8;

	//+-3*pi/8
	const float _3pi8=3*3.14f/8;
	const float _n3pi8=-3*3.14f/8;

	//+-5*pi/8
	const float _5pi8=5*3.14f/8;
	const float _n5pi8=-5*3.14f/8;

	//+-7*pi/8
	const float _7pi8=7*3.14f/8;
	const float _n7pi8=-7*3.14f/8;

	for( int y=1; y< img_in_magnitude.Height()-1; y++ ) 
	{
		for( int x=1; x< img_in_magnitude.Width()-1; x++ ) 
		{
			float angle=0;
			if (img_in_angle(x,y)<0 )  angle=img_in_angle(x,y)+3.14f;
			else angle=img_in_angle(x,y);

			if (angle<_pi8)// 0-pi/8  <-7pi/8
			{
				if(img_in_magnitude(x,y) < img_in_magnitude(x-1,y) || img_in_magnitude(x,y)< img_in_magnitude(x+1,y))
					(*img_out)(x,y)=0 ;
				else (*img_out)(x,y)=img_in_magnitude(x,y);
			}
				

			else if (angle<_3pi8)// pi/8-3pi/8   <-5pi/8  >-7pi/8
			{
				if(img_in_magnitude(x,y) < img_in_magnitude(x-1,y-1) || img_in_magnitude(x,y)< img_in_magnitude(x+1,y+1))
					(*img_out)(x,y)=0 ;
				else (*img_out)(x,y)=img_in_magnitude(x,y);
				
			}

			else if (angle<_5pi8)// 3pi/8-5pi/8   <-3pi/8 >-5pi/8
			{
				if(img_in_magnitude(x,y) < img_in_magnitude(x,y-1) || img_in_magnitude(x,y)< img_in_magnitude(x,y+1))
					(*img_out)(x,y)=0 ;
				else (*img_out)(x,y)=img_in_magnitude(x,y);
			}


			else if (angle<_7pi8)// 5pi/8-7pi/8   <-pi/8 >-3pi/8
			{
				if(img_in_magnitude(x,y) < img_in_magnitude(x+1,y-1) || img_in_magnitude(x,y)< img_in_magnitude(x-1,y+1))
					(*img_out)(x,y)=0 ;
				else (*img_out)(x,y)=img_in_magnitude(x,y);
				
				
			}

			else// 0-  -pi/8
			{	
				if(img_in_magnitude(x,y) < img_in_magnitude(x-1,y) || img_in_magnitude(x,y)< img_in_magnitude(x+1,y))
					(*img_out)(x,y)=0 ;
				else (*img_out)(x,y)=img_in_magnitude(x,y);
			}
					
		}//for x
	}//for y

	return 0;
	
}


//use cross to floodfill
void connect_binary_cross( const ImgBinary &img_hi, const ImgBinary &img_lo, const int x, const int y, ImgBinary* img_out )
{
	ImgInt point_xy;
	point_xy.Reset(img_hi.Height() * img_hi.Width()+1,2);
	Set(&point_xy,0);

	ImgInt point_fuc;
	point_fuc.Reset(img_hi.Height() * img_hi.Width(),1);
	Set(&point_fuc,0);
	
	int loop=1;
	long int ctn=2;
	
	
	(*img_out)(x,y) = 1;
	point_xy(1,0)=x;
	point_xy(1,1)=y;
	point_fuc(1,0)=1;

	while ( loop )
	{

		int status=0;
		int current_x=0;
		int current_y=0;
		int current_search=0;


		for(int search_point=1; search_point<img_hi.Height() * img_hi.Width()+1; search_point++)
		{
			//find seed point
			if( point_fuc(search_point,0)==1)
			{
				current_x=point_xy(search_point,0);
				current_y=point_xy(search_point,1);
				current_search=search_point;
			}


		}
		
		//exit 
		if(current_search ==0)
		{
			break;	  
		}

				
		//cross add
		if (current_y > 0 && img_lo(current_x,current_y-1) == 1 && (*img_out)(current_x,current_y-1) != 1)  
		{	


			point_xy(ctn,0)=current_x;
			point_xy(ctn,1)=current_y-1;
			point_fuc(ctn)=1;
			(*img_out)(current_x,current_y-1) = 1;
			ctn++;

		}
		else status+=1;

		if (current_y < img_hi.Height() &&img_lo(current_x,current_y+1) == 1 && (*img_out)(current_x,current_y+1) != 1)  
		{	

			point_xy(ctn,0)=current_x;
			point_xy(ctn,1)=current_y+1;
			point_fuc(ctn)=1;
			(*img_out)(current_x,current_y+1) = 1;
			ctn++;

		}
		else status+=1;



		if (current_x > 0 && img_lo(current_x-1,current_y) == 1 && (*img_out)(current_x-1,current_y) != 1)  
		{	


			point_xy(ctn,0)=current_x-1;
			point_xy(ctn,1)=current_y;
			point_fuc(ctn)=1;
			(*img_out)(current_x-1,current_y) = 1;
			ctn++;
		}
		else status+=1;


		if (current_x < img_hi.Width() && img_lo(current_x+1,current_y) == 1 && (*img_out)(current_x+1,current_y) != 1)  
		{	

			point_xy(ctn,0)=current_x+1;
			point_xy(ctn,1)=current_y;
			point_fuc(ctn)=1;
			(*img_out)(current_x+1,current_y) = 1;
			ctn++;

		}
		else status+=1;
		
		point_fuc(current_search,0)=0;
		//printf("  status=%d  ",status);
		
	}//while
}//func





void connect_outline3x3( const ImgBinary &img_hi, const ImgBinary &img_lo, const int x, const int y, ImgBinary* img_out )
{

	for( int y=0; y< img_lo.Height(); y++ ) 
	{
		for( int x=0; x< img_lo.Width(); x++ ) 
		{

			if ((*img_out)(x,y)==1) continue;
			
			//operate --x
			
			int offset_x_left=0;
			

			while(1)
			{					
				if( x-offset_x_left>=0 )
				{
					
					if(img_lo(x-offset_x_left,y)!=0 && (*img_out)(x-offset_x_left,y)!=1)
					{
						(*img_out)(x-offset_x_left,y)=1;
						
						int offset_y_up_a=0;
						int offset_y_dw_a=0;


						//opt --y
						while(1)
						{				
							if( y-offset_y_up_a>=0)
							{

								if(img_lo(x-offset_x_left,y-offset_y_up_a)!=0 && (*img_out)(x-offset_x_left,y-offset_y_up_a)!=1 )
								{
									(*img_out)(x-offset_x_left,y-offset_y_up_a)=1;
								}
								else
								{
									break;					
								}
								offset_y_up_a++;
							}
							else
							{
									break;					
							}
						}//opt  --y



						//opt y++
						while(1)
						{				
							if( y+offset_y_dw_a<img_lo.Height())
							{

								if(img_lo(x-offset_x_left,y+offset_y_dw_a)!=0 && (*img_out)(x-offset_x_left,y+offset_y_dw_a)!=1)
								{
									(*img_out)(x-offset_x_left,y+offset_y_dw_a)=1;
								}
								else
								{
									break;					
								}

								offset_y_dw_a++;
							}
							else
							{
									break;					
							}
						}//opt  y++
				


					}
					else
					{
						break;					
					}
					offset_x_left++;
				}
				else
				{
					break;					
				}
			}
			
			
			
			
			int offset_x_right=0;
			


			//operate x++
			while(1)
			{				

				if( x+offset_x_right<img_lo.Width())
				{
					
					if(img_lo(x+offset_x_right,y)!=0 && (*img_out)(x+offset_x_right,y)!=1)
					{
						(*img_out)(x+offset_x_right,y)=1;

						int offset_y_up_b=0;
						int offset_y_dw_b=0;
						//opt --y
						while(1)
						{				
							if( y-offset_y_up_b>=0)
							{

								if(img_lo(x+offset_x_right,y-offset_y_up_b)!=0 && (*img_out)(x+offset_x_right,y-offset_y_up_b)!=1)
								{
									(*img_out)(x+offset_x_right,y-offset_y_up_b)=1;
								}
								else
								{
									break;					
								}
								offset_y_up_b++;
							}
							else
							{
									break;					
							}

						}//opt  --y



						//opt y++
						while(1)
						{				
							if( y+offset_y_dw_b<img_lo.Height())
							{

								if(img_lo(x+offset_x_right,y+offset_y_dw_b)!=0 &&(*img_out)(x+offset_x_right,y+offset_y_dw_b)!=1 )
								{
									(*img_out)(x+offset_x_right,y+offset_y_dw_b)=1;
								}
								else
								{
									break;					
								}

								offset_y_dw_b++;
								
							}
							else
							{
									break;					
							}

						}//opt  y++

					}
					else
					{
						break;					
					}
					offset_x_right++;
				}
				else
				{
					break;					
				}


			}
			
		}// for x
	}//for y

}


//bubble sort
int bub_sort(ImgFloat* sort_operating , const int s_length)
{
	if( (*sort_operating).Width()==0 &&(*sort_operating).Height()==0 )
		return 1;


	float swap_value=0;

	for(int i=0;i<s_length;i++)
	{	
		for(int j=0;j<s_length - i-1;j++)
		{
			if( (*sort_operating)(j,0)>(*sort_operating)(j+1,0)  )
			{
				swap_value=(*sort_operating)(j,0);
				(*sort_operating)(j,0)=(*sort_operating)(j+1,0);
				(*sort_operating)(j+1,0)=swap_value;
			}
		}// for j
	}//for i
	return 0;
}



//Edge linking with hysteresis
int edge_linking( const ImgFloat &img_in_magnitude, float* treshold_hi_val, float* treshold_lo_val, ImgBinary* img_treshold_mix) 
{
	/*int sum_pixel_img_1= (img_in_magnitude.Width() * img_in_magnitude.Height())/10;
	int sum_pixel_img_2= (img_in_magnitude.Width() * img_in_magnitude.Height())/9;*/
	//float cell=10;


	ImgFloat values_mag;
	int values_mag_length=0;
	
	
	float ppercent_val= 0.8f;
	float threshold_ratio=5;
	int ppercent_position=0;


	values_mag.Reset( img_in_magnitude.Width() * img_in_magnitude.Height() ,1);
	Set(&values_mag,0);
	


	for( int y=0; y< img_in_magnitude.Height(); y++ ) 
	{
		for( int x=0; x< img_in_magnitude.Width(); x++ ) 
		{				
			if(img_in_magnitude(x,y)>0 ) 
				values_mag(values_mag_length++, 0)=img_in_magnitude(x,y);

		}//for x
	}//for y
	if (bub_sort(&values_mag , values_mag_length)==1) return 1;



	ppercent_position = Round( ppercent_val * values_mag_length );

	*treshold_hi_val=values_mag(ppercent_position, 0);
	*treshold_lo_val=*treshold_hi_val/threshold_ratio;

	/*printf("\n treshold_hi_val= %f \n",*treshold_hi_val);
	printf("\n treshold_lo_val= %f \n",*treshold_lo_val);	*/


	////test hysteresis p=10, 10%
	//printf("sum_pixel_img_1= %d\n",sum_pixel_img_1);
	//printf("sum_pixel_img_2= %d\n",sum_pixel_img_2);

	//for( int ctn=0; ctn<50; ctn++)
	//{
	//	printf("sum_values_mag[%d]= %d\n", ctn,sum_values_mag[ctn]);

	//}

	ImgBinary img_treshold_hi;
	ImgBinary img_treshold_lo;
	//ImgBinary img_treshold_mix;
	
	img_treshold_hi.Reset(img_in_magnitude.Width(), img_in_magnitude.Height() );
	Set(&img_treshold_hi,0);

	img_treshold_lo.Reset(img_in_magnitude.Width(), img_in_magnitude.Height() );
	Set(&img_treshold_lo,0);
	
	(*img_treshold_mix).Reset(img_in_magnitude.Width(), img_in_magnitude.Height() );
	Set(img_treshold_mix,0);

	/**treshold_hi_val=10.0;
	*treshold_lo_val=*treshold_hi_val/5;*/

	for( int y=1; y< img_in_magnitude.Height()-1; y++ ) 
	{
		for( int x=1; x< img_in_magnitude.Width()-1; x++ ) 
		{		
			if( img_in_magnitude(x,y)>*treshold_lo_val)
			{	
				img_treshold_lo(x,y)=1;
			}
			if( img_in_magnitude(x,y)>*treshold_hi_val)
			{	
				
				img_treshold_hi(x,y)=1;
			}
		}
	}
	

	for( int y=1; y< img_in_magnitude.Height()-1; y++ ) 
	{
		for( int x=1; x< img_in_magnitude.Width()-1; x++ ) 
		{	
			if( img_treshold_hi(x,y)==1)
			{	
				//connect_outline3x3( img_treshold_hi, img_treshold_lo ,x,y, img_treshold_mix );
				//FloodFill8( img_treshold_lo, x,  y, 5.0, img_treshold_mix);
				connect_binary_cross( img_treshold_hi, img_treshold_lo ,x,y, img_treshold_mix );
			}
		}//for x
	}//for y

	//test
	/*Figure lo,hi;
	lo.SetTitle("lo");
	hi.SetTitle("hi");
	lo.Draw(img_treshold_lo);
	hi.Draw(img_treshold_hi);*/
	return 0;

}



int Chamfering_distance_compute(const ImgFloat &img_in, const int ratio_chamfer, ImgFloat* img_out)
{
	float offset_value=100;
	float maximun_pixel = offset_value * ratio_chamfer;

	(*img_out).Reset(img_in.Width(), img_in.Height());
	//(*img_out)=img_in;
	
	//float median_pixel_chamfer = (img_in.Width()+img_in.Height() )/2*0.1f/2;
	
	Set(img_out,maximun_pixel);


	//first filter
	for( int y=1; y< img_in.Height(); y++ ) 
	{
		for( int x=1; x< img_in.Width(); x++ ) 
		{	
			
			if (img_in(x, y)==0) 
			{
				//(*img_out)(x,y) = img_in(x,y)?0 :min(ininity, 1+d(x-1,y), 1+d(x,y-1));
				float comparation_pixel = maximun_pixel;
				comparation_pixel = min( comparation_pixel, (*img_out)(x+1, y-1) + offset_value);
				comparation_pixel = min( comparation_pixel, (*img_out)(x-1, y-1) + offset_value);
				comparation_pixel = min( comparation_pixel, (*img_out)(x, y-1) + offset_value);
				comparation_pixel = min( comparation_pixel, (*img_out)(x-1, y) + offset_value);
				//comparation_pixel = min( comparation_pixel, (*img_out)(x-1, y+1) + 0.1f);
				(*img_out)(x,y) = comparation_pixel;
			}

			else  (*img_out)(x, y) = 0;
      			
			
		}//for x
	}//for y





	//Second filter
	for( int y=img_in.Height()-2; y>0 ; y-- ) 
	{
		for( int x=img_in.Width()-2; x>0; x-- ) 
		{	
			
			if (img_in(x, y)==0) 
			{
				//(*img_out)= img_in(x,y) ? 0 : min( d(x,y), 1+d(x+1,y), 1+d(x,y+1) );
				float comparation_pixel = (*img_out)(x, y);
				comparation_pixel = min( comparation_pixel, (*img_out)(x-1, y+1) + offset_value);
				comparation_pixel = min( comparation_pixel, (*img_out)(x+1, y+1) + offset_value);
				comparation_pixel = min( comparation_pixel, (*img_out)(x, y+1) + offset_value);
				comparation_pixel = min( comparation_pixel, (*img_out)(x+1, y) + offset_value);
				//comparation_pixel = min( comparation_pixel, (*img_out)(x+1, y-1) + 0.1f);
				(*img_out)(x,y) = comparation_pixel;
			}
						
			else  (*img_out)(x, y) = 0;
			
		}//for x
	}//for y


	return 0;
}




int canny_edge_dect(const ImgBgr &img_in, const int kernel_width ,const float gaussian_kernel[] ,const float gaussian_deri_kernel[] , ImgFloat* img_out)
{
	ImgGray img_mid;
	ImgFloat img_convolution_xd;
	ImgFloat img_convolution_yd;

	ImgFloat img_magnitude;
	ImgFloat img_angle;
	ImgFloat img_nms;

	ImgBinary img_out_mid;

	float treshold_hi_val=0;
	float treshold_lo_val=0;
	Convert(img_in,&img_mid);

	
	int status_value=0;
	
	//Gd x
	if ( separable_convolution_2d( img_mid, kernel_width , gaussian_kernel, gaussian_deri_kernel , &img_convolution_xd , 1) ==0)   status_value+=1;
	
	//Gd y
	if ( separable_convolution_2d( img_mid, kernel_width , gaussian_kernel , gaussian_deri_kernel, &img_convolution_yd , 0) ==0) status_value+=2;




	if ( grad_magnitude_angle( img_convolution_xd,  img_convolution_yd,  &img_magnitude,  &img_angle ) ==0) status_value+=3;

	if ( nms_filter( img_magnitude, img_angle, &img_nms)  ==0)  status_value+=4;


	if (  edge_linking( img_nms, &treshold_hi_val,&treshold_lo_val, &img_out_mid) ==0)  status_value+=5;

	Convert(img_out_mid,img_out);


	return status_value;

}








int item_dect(const ImgBgr &img_origin, const ImgFloat &img_chamfer ,const ImgFloat &img_templete, ImgFloat* img_probability, ImgBgr* img_out_circle )
{
	if( img_templete.Width()>=img_chamfer.Width() || img_templete.Height() >= img_chamfer.Height() )
		return 1;


	(*img_probability).Reset( img_origin.Width()-img_templete.Width() , img_origin.Height()-img_templete.Height() );
	Set(img_probability,0);

	(*img_out_circle).Reset( img_origin.Width() , img_origin.Height() );
	(*img_out_circle)=img_origin;
	
	

	float min_possibility=0; 
	int refence_point[2]={0};
	
	for( int y=0 ; y < img_origin.Height()-img_templete.Height()-1 ; y++) 
	{
		for( int x=0 ; x < img_origin.Width()-img_templete.Width()-1 ; x++ ) 
		{


			for(int xt=0; xt< img_templete.Width(); xt++ ) 
			{
				for(  int yt=0; yt< img_templete.Height(); yt++ ) 
				{
					if( img_templete(xt,yt)>0  )
						(*img_probability)(x,y)+= img_chamfer(xt+x,yt+y);

				}//for xt
			}//for yt

			if( x==0 && y==0) 
			{	
				min_possibility=(*img_probability)(x,y);
				refence_point[0]=x;
				refence_point[1]=y;
			}
			min_possibility<(*img_probability)(x,y)? min_possibility:( min_possibility=(*img_probability)(x,y), refence_point[0]=x,	refence_point[1]=y );

		}//for x
	}//for y
	Rect can_rect;
	Point left_top( refence_point[0], refence_point[1] );
	Point Right_bottom( refence_point[0]+img_templete.Width(), refence_point[1]+img_templete.Height() );

	can_rect=CRect(left_top,Right_bottom );

	DrawRect(can_rect, img_out_circle, Bgr(255,0,0), 1);

	


	return 0;

}

