#include <algorithm>
using namespace std;
#define QX_DEF_CHAR_MAX 255

/***********************************************************************
\Author:	Qingxiong Yang
\Function:	Recursive Bilateral Filtering
\Reference:	Qingxiong Yang, Recursive Bilateral Filtering,
European Conference on Computer Vision (ECCV) 2012, 399-413.
************************************************************************/


//Large sigma_spatial/sigma_range parameter may results in visible artifact which can be removed by 
//an additional filter with small sigma_spatial/sigma_range parameter.

unsigned char * recursive_bf(
	unsigned char * texture,
	double sigma_spatial, double sigma_range, 
	int width, int height, int channel)
{
	double * out = new double[width * height * channel];
	double * temp = new double[width * height * channel];
	double * temp_factor = new double[width * (height * 2 + 2)];
	double * factor_a = new double[width * height];
	double * factor_b = new double[width];
	double * factor_c = new double[width];
	double * factor_cb = new double[width * channel];
	double * factor_cc = new double[width * channel];

	double * in = new double[width * height * channel];
	for (int i = 0; i < width * height * channel; ++i)
		in[i] = texture[i];
	
	//compute a lookup table
	double range_table[QX_DEF_CHAR_MAX + 1];
	double inv_sigma_range = 1.0 / (sigma_range * QX_DEF_CHAR_MAX);
	for (int i = 0; i <= QX_DEF_CHAR_MAX; i++) 
		range_table[i] = exp(-i * inv_sigma_range);

	double alpha = exp(-sqrt(2.0) / (sigma_spatial * width));
	double ypr, ypg, ypb, ycr, ycg, ycb;
	double fp, fc;
	double inv_alpha_ = 1 - alpha;
	for (int y = 0; y < height; y++)
	{
		double * temp_x = &temp[y * width * channel];
		double * in_x = &in[y * width * channel];
		unsigned char * texture_x = &texture[y * width * channel];
		*temp_x++ = ypr = *in_x++; 
		*temp_x++ = ypg = *in_x++; 
		*temp_x++ = ypb = *in_x++;
		unsigned char tpr = *texture_x++; 
		unsigned char tpg = *texture_x++;
		unsigned char tpb = *texture_x++;

		double * temp_factor_x = &temp_factor[y * width];
		*temp_factor_x++ = fp = 1;
		for (int x = 1; x < width; x++) // from left to right
		{
			unsigned char tcr = *texture_x++; 
			unsigned char tcg = *texture_x++; 
			unsigned char tcb = *texture_x++;
			unsigned char dr = abs(tcr - tpr);
			unsigned char dg = abs(tcg - tpg);
			unsigned char db = abs(tcb - tpb);
			int range_dist = (((dr << 1) + dg + db) >> 2);
			double weight = range_table[range_dist];
			double alpha_ = weight*alpha;
			*temp_x++ = ycr = inv_alpha_*(*in_x++) + alpha_*ypr; 
			*temp_x++ = ycg = inv_alpha_*(*in_x++) + alpha_*ypg; 
			*temp_x++ = ycb = inv_alpha_*(*in_x++) + alpha_*ypb;
			tpr = tcr; tpg = tcg; tpb = tcb;
			ypr = ycr; ypg = ycg; ypb = ycb;
			*temp_factor_x++ = fc = inv_alpha_ + alpha_*fp;
			fp = fc;
		}
		*--temp_x; *temp_x = 0.5*((*temp_x) + (*--in_x));
		*--temp_x; *temp_x = 0.5*((*temp_x) + (*--in_x));
		*--temp_x; *temp_x = 0.5*((*temp_x) + (*--in_x));
		tpr = *--texture_x; 
		tpg = *--texture_x; 
		tpb = *--texture_x;
		ypr = *in_x; ypg = *in_x; ypb = *in_x;

		*--temp_factor_x; *temp_factor_x = 0.5*((*temp_factor_x) + 1);
		fp = 1;

		for (int x = width - 2; x >= 0; x--) // from right to left
		{
			unsigned char tcr = *--texture_x; 
			unsigned char tcg = *--texture_x; 
			unsigned char tcb = *--texture_x;
			unsigned char dr = abs(tcr - tpr);
			unsigned char dg = abs(tcg - tpg);
			unsigned char db = abs(tcb - tpb);
			int range_dist = (((dr << 1) + dg + db) >> 2);
			double weight = range_table[range_dist];
			double alpha_ = weight * alpha;

			ycr = inv_alpha_ * (*--in_x) + alpha_ * ypr; 
			ycg = inv_alpha_ * (*--in_x) + alpha_ * ypg; 
			ycb = inv_alpha_ * (*--in_x) + alpha_ * ypb;
			*--temp_x; *temp_x = 0.5*((*temp_x) + ycr);
			*--temp_x; *temp_x = 0.5*((*temp_x) + ycg);
			*--temp_x; *temp_x = 0.5*((*temp_x) + ycb);
			tpr = tcr; tpg = tcg; tpb = tcb;
			ypr = ycr; ypg = ycg; ypb = ycb;

			fc = inv_alpha_ + alpha_*fp;
			*--temp_factor_x; 
			*temp_factor_x = 0.5*((*temp_factor_x) + fc);
			fp = fc;
		}
	}
	alpha = exp(-sqrt(2.0) / (sigma_spatial * height));
	inv_alpha_ = 1 - alpha;
	double * ycy, * ypy, * xcy;
	unsigned char * tcy, * tpy;
	memcpy(out, temp, sizeof(double)* width * channel);

	double * in_factor = temp_factor;
	double*ycf, *ypf, *xcf;
	memcpy(factor_a, in_factor, sizeof(double) * width);
	for (int y = 1; y < height; y++)
	{
		tpy = &texture[(y - 1) * width * channel];
		tcy = &texture[y * width * channel];
		xcy = &temp[y * width * channel];
		ypy = &out[(y - 1) * width * channel];
		ycy = &out[y * width * channel];

		xcf = &in_factor[y * width];
		ypf = &factor_a[(y - 1) * width];
		ycf = &factor_a[y * width];
		for (int x = 0; x < width; x++)
		{
			unsigned char dr = abs((*tcy++) - (*tpy++));
			unsigned char dg = abs((*tcy++) - (*tpy++));
			unsigned char db = abs((*tcy++) - (*tpy++));
			int range_dist = (((dr << 1) + dg + db) >> 2);
			double weight = range_table[range_dist];
			double alpha_ = weight*alpha;
			for (int c = 0; c < channel; c++) 
				*ycy++ = inv_alpha_*(*xcy++) + alpha_*(*ypy++);
			*ycf++ = inv_alpha_*(*xcf++) + alpha_*(*ypf++);
		}
	}
	int h1 = height - 1;
	ycf = factor_b;
	ypf = factor_c;
	memcpy(ypf, &in_factor[h1 * width], sizeof(double) * width);
	for (int x = 0; x < width; x++) 
		factor_a[h1 * width + x] = 0.5*(factor_a[h1 * width + x] + ypf[x]);

	ycy = factor_cb;
	ypy = factor_cc;
	memcpy(ypy, &temp[h1 * width * channel], sizeof(double)* width * channel);
	int k = 0; 
	for (int x = 0; x < width; x++) {
		for (int c = 0; c < channel; c++) {
			int idx = (h1 * width + x) * channel + c;
			out[idx] = 0.5*(out[idx] + ypy[k++]) / factor_a[h1 * width + x];
		}
	}

	for (int y = h1 - 1; y >= 0; y--)
	{
		tpy = &texture[(y + 1) * width * channel];
		tcy = &texture[y * width * channel];
		xcy = &temp[y * width * channel];
		double*ycy_ = ycy;
		double*ypy_ = ypy;
		double*out_ = &out[y * width * channel];

		xcf = &in_factor[y * width];
		double*ycf_ = ycf;
		double*ypf_ = ypf;
		double*factor_ = &factor_a[y * width];
		for (int x = 0; x < width; x++)
		{
			unsigned char dr = abs((*tcy++) - (*tpy++));
			unsigned char dg = abs((*tcy++) - (*tpy++));
			unsigned char db = abs((*tcy++) - (*tpy++));
			int range_dist = (((dr << 1) + dg + db) >> 2);
			double weight = range_table[range_dist];
			double alpha_ = weight*alpha;

			double fcc = inv_alpha_*(*xcf++) + alpha_*(*ypf_++);
			*ycf_++ = fcc;
			*factor_ = 0.5*(*factor_ + fcc);

			for (int c = 0; c < channel; c++)
			{
				double ycc = inv_alpha_*(*xcy++) + alpha_*(*ypy_++);
				*ycy_++ = ycc;
				*out_ = 0.5*(*out_ + ycc) / (*factor_);
				*out_++;
			}
			*factor_++;
		}
		memcpy(ypy, ycy, sizeof(double) * width * channel);
		memcpy(ypf, ycf, sizeof(double) * width);
	}


	unsigned char * out_img = new unsigned char[width * height * channel];
	for (int i = 0; i < width * height * channel; ++i)
		out_img[i] = static_cast<unsigned char>(out[i]);

	delete[] in;
	delete[] out;
	delete[] temp;
	delete[] temp_factor;
	delete[] factor_a;
	delete[] factor_b;
	delete[] factor_c;
	delete[] factor_cb;
	delete[] factor_cc;

	return out_img;
}