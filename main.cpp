


//Vs2012 update4 multi-bytes


#include <afxwin.h>  // necessary for MFC to work properly
#include <math.h> 
#include "main.h"
#include "../../src/blepo.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif





using namespace blepo;
using namespace blepo_ex;


int main(int argc, const char* argv[], const char* envp[])
{
	int third_commandpara_specified=0;
	if( argc > 4)	
	{
		printf("Error: more command parameters than use \n");
		third_commandpara_specified=1;		
	}
	
	
	else
	{
		printf("argc = %d\n", argc);
		for (int i=0; i<argc ; i++) 
		{
			printf("argv[%d]=%s\n\n", i, argv[i]);
		}
	}
	if( argc == 4 )
	{
		third_commandpara_specified=1;		
	}

	HMODULE hModule = ::GetModuleHandle(NULL);
	if (hModule == NULL || !AfxWinInit(hModule, NULL, ::GetCommandLine(), 0))
	{
		printf("Fatal Error: MFC initialization failed (hModule = %x)\n", hModule);
		return 1;
	}


	/*
	a.	Reads 1 command-line parameter, which we will call filename.  
	b.	Loads filename from the blepo/images directory into a Grayscale or BGR image and displays it in a figure window.
	*/
	CString path_blepo_images ="../../images/";
	CString filename_input1=argv[1];
	CString filename_input2=argv[2];
	CString filename_input3=argv[3];

	CString filename_forloading_sigma = argv[1];
	CString filename_forloading_origin = path_blepo_images + argv[2];
	CString filename_forloading_templete = path_blepo_images + argv[3];
	
	//Judge File name that is whether incomplete before loading
	int status_filename1=initial_filename_recognize(path_blepo_images, filename_input2);
	
	if( status_filename1==1 )
	{
		printf("Error: File cannot be found : %s \n",filename_forloading_origin );	
		printf("Act: Program halt, Please close \n");
		EventLoop();
		return 0;
	}
	
	printf("\n Start: \n");
	//Val initial 
	float sigma_value=_ttof(filename_input1);
	/*filename_input1 = _T("655350");
	_stscanf(filename_input1, _T("%f"), &sigma_value);*/


	printf(" Sigma = %f \n", sigma_value);
		
	ImgBgr img_loaded_origin_bgr;
	

	//Figure and title
	Figure fig_loaded_image;
	fig_loaded_image.SetTitle("Loaded Image Origin");
	
	//a.
	// Loading filename2 3 and show it 
	Load(filename_forloading_origin, &img_loaded_origin_bgr);
	fig_loaded_image.Draw(img_loaded_origin_bgr);


	/*b.	Construct a 1D Gaussian convolution kernel and a 1D Gaussian derivative convolution kernel, 
		both with the specified value for sigma.  
		The length of each kernel should be automatically selected as explained in the notes.*/

	float f=2.5;
	int kernel_halfwidth= Round( f * sigma_value -0.5); //
	int kernel_width=2*kernel_halfwidth +1;

	//input protection
	if(kernel_width>=100 || sigma_value==0)
	{
		sigma_value=2.0;
		kernel_halfwidth=5;
		kernel_width=11;
		printf("Error: Input Sigma is wrong or too large.\n");
		printf("Act: go back to the default value Sigma=2.0.\n");
		
	}


	float gaussian_kernel[100]={0};
	float gaussian_sum=0;
	float sum1=0;
	
	//Build Gaussian  kernel
	for(int ctn=0; ctn<kernel_width ; ctn++)
	{
		gaussian_kernel[ctn] = exp( 
						-(ctn-kernel_halfwidth) * (ctn-kernel_halfwidth)  
					   /(2 * sigma_value * sigma_value) );
		gaussian_sum+=gaussian_kernel[ctn];
	
	}

	for(int ctn=0; ctn<kernel_width ; ctn++)
	{
		gaussian_kernel[ctn]/= gaussian_sum;
		sum1+=gaussian_kernel[ctn];
	
	}
	printf("b.\n");
	printf("kernel_halfwidth = %d \n", kernel_halfwidth);
	printf("gaussian_sum = %f \n", gaussian_sum);
	printf("Sum1 = %f \n", sum1);
	for(int ctn=0; ctn<kernel_width ; ctn++)
	{
		printf("gaussian_kernel[%d] = %f \n",ctn , gaussian_kernel[ctn]);
	
	
	}


	//Build Gaussian derivative kernel
	float gaussian_deri_kernel[100]={0};
	float gaussian_deri_sum=0;
	float sum2=0;

	for(int ctn=0; ctn<kernel_width ; ctn++)
	{
		gaussian_deri_kernel[ctn] = (ctn-kernel_halfwidth)
				*exp(  -(ctn-kernel_halfwidth) * (ctn-kernel_halfwidth)  
					   /(2 * sigma_value * sigma_value) );

		gaussian_deri_sum+=gaussian_deri_kernel[ctn] * ctn;
	
	}

	for(int ctn=0; ctn<kernel_width ; ctn++)
	{
		gaussian_deri_kernel[ctn]= (-1)*gaussian_deri_kernel[ctn]/gaussian_deri_sum;
		sum2+=gaussian_deri_kernel[ctn];
	
	}

	
	printf("\n\n");
	printf("gaussian_deri_sum = %f \n", gaussian_deri_sum);
	printf("Sum2 = %f \n", sum2);
	printf("With flip \n");
	for(int ctn=0; ctn<kernel_width ; ctn++)
	{
		
		printf("gaussian_deri_kernel[%d] = %f \n",ctn , gaussian_deri_kernel[ctn]);
	}
	

	/*c.	Compute the gradient of the image.  
		Use the principle of separability to apply the 1D Gaussian and 1D Gaussian derivative kernels.  
		Do not worry about image borders; the simplest solution is to simply set the border pixels in the convolution result to zero 
		rather than extending the image, but extension is fine, too.
*/
	printf("\nc. \n");
	ImgGray img_loaded_origin_flo;
	ImgFloat img_convolution_xd;
	ImgFloat img_convolution_yd;

	ImgFloat img_magnitude;
	ImgFloat img_angle;
	ImgFloat img_nms;

	Convert(img_loaded_origin_bgr,&img_loaded_origin_flo);

	//Gd x
	if ( separable_convolution_2d( img_loaded_origin_flo, kernel_width , gaussian_kernel ,gaussian_deri_kernel , &img_convolution_xd , 1) ==0)   printf(" Gx \n");
	
	//Gd y
	if ( separable_convolution_2d( img_loaded_origin_flo, kernel_width , gaussian_kernel , gaussian_deri_kernel, &img_convolution_yd , 0) ==0) printf(" Gy \n");


	if ( grad_magnitude_angle( img_convolution_xd,  img_convolution_yd,  &img_magnitude,  &img_angle ) ==0) printf("\n mag and angle \n");


  //show fig
	Figure gradient_x,gradient_y,gradient_magnitude, gradient_angle;
	Figure nms,threshold_mix,chamgering_distance;
	Figure fig_templete_edge, fig_probability_detection, fig_item_detection;

	gradient_x.SetTitle("gradient_x");
	gradient_y.SetTitle("gradient_y");
	gradient_magnitude.SetTitle("gradient_magnitude");
	gradient_angle.SetTitle("gradient_angle");
	
	nms.SetTitle("non-maximum supression");
	threshold_mix.SetTitle("double threshold");
	chamgering_distance.SetTitle("chamgering_distance");
	fig_templete_edge.SetTitle("templete_edge");
	fig_probability_detection.SetTitle("probability_detection");
	fig_item_detection.SetTitle("item_detection");

	gradient_x.Draw(img_convolution_xd);
	gradient_y.Draw(img_convolution_yd);
	gradient_magnitude.Draw(img_magnitude);
	gradient_angle.Draw(img_angle);
	

	//d.	Perform non-maximum suppression using the gradient magnitude and phase.
	if ( nms_filter( img_magnitude, img_angle, &img_nms)  ==0)  printf("d. \n");

	nms.Draw(img_nms);



	//e.	Perform thresholding with hysteresis (i.e., double-thresholding).  Automatically compute the threshold values based upon image statistics.
	float treshold_hi_val=0;
	float treshold_lo_val=0;

	ImgBinary img_threshold_mix;

	if (  edge_linking( img_nms, &treshold_hi_val,&treshold_lo_val, &img_threshold_mix) ==0)  printf("e. \n");

	threshold_mix.Draw(img_threshold_mix);


	/*f.	Compute the chamfer distance on the edge image with the Manhattan distance.
*/

	
	/*ImgGray imgx;
	ImgInt imgy;
	Convert(img_threshold_mix, &imgx);

	Chamfer(imgx, &imgy);
	chamgering_distance.Draw(imgy);*/

	 ImgFloat img_threshold_mix_f;
	 ImgFloat img_chamgering_distance;

	 Convert(img_threshold_mix,&img_threshold_mix_f);

	 const int ratio_chamfer=5;

	 if (  Chamfering_distance_compute(img_threshold_mix_f, ratio_chamfer,&img_chamgering_distance) ==0 )  printf("f. \n");


	 chamgering_distance.Draw(img_chamgering_distance);


	//g
	if( 1==third_commandpara_specified )
	{
		//test template name is right to open
		int status_filename2=initial_filename_recognize(path_blepo_images, filename_input3);
		if( status_filename2==1 )
		{
			printf("Error: File cannot be found : %s \n",filename_forloading_templete );	
			printf("Act: Program halt, Please close \n");
			EventLoop();
			return 0;
		}

		//load templete
		ImgBgr img_loaded_template_bgr;
		Load(filename_forloading_templete, &img_loaded_template_bgr);
		
		Figure fig_loaded_template;
		fig_loaded_template.SetTitle("Loaded Image Template");
		fig_loaded_template.Draw(img_loaded_template_bgr);


		//canny edge processing
		ImgFloat img_templete_edge;

		if (  canny_edge_dect(img_loaded_template_bgr, kernel_width ,gaussian_kernel ,gaussian_deri_kernel , &img_templete_edge) ==15 )  printf("g. \n");


		fig_templete_edge.Draw(img_templete_edge);

		//Matching
		ImgFloat img_probability;
		ImgBgr img_out_circle;


		if (   item_dect(img_loaded_origin_bgr, img_chamgering_distance , img_templete_edge , &img_probability, &img_out_circle ) ==0 )  printf("h. \n");
		else printf("Error: too large size for templete image");


		fig_probability_detection.Draw(img_probability);
		fig_item_detection.Draw(img_out_circle);

	}//if 3 command parameter

	EventLoop();
	return 0;
}



