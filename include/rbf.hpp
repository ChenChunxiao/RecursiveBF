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
	float sigma_spatial, float sigma_range, 
	int width, int height, int channel)
{
	const int width_height = width * height;
	const int width_channel = width * channel;
	const int width_height_channel = width * height * channel;

	float * out = new float[width_height_channel];
	float * temp = new float[width_height_channel];
	float * temp_factor = new float[width * (height * 2 + 2)];
	float * factor_a = new float[width_height];
	float * factor_b = new float[width];
	float * factor_c = new float[width];
	float * factor_cb = new float[width_channel];
	float * factor_cc = new float[width_channel];
	
	//compute a lookup table
	float range_table[QX_DEF_CHAR_MAX + 1];
	float inv_sigma_range = 1.0f / (sigma_range * QX_DEF_CHAR_MAX);
	for (int i = 0; i <= QX_DEF_CHAR_MAX; i++) 
		range_table[i] = exp(-i * inv_sigma_range);

	float alpha = static_cast<float>(exp(-sqrt(2.0) / (sigma_spatial * width)));
	float ypr, ypg, ypb, ycr, ycg, ycb;
	float fp, fc;
	float inv_alpha_ = 1 - alpha;
	for (int y = 0; y < height; y++)
	{
		float * temp_x = &temp[y * width_channel];
		unsigned char * in_x = &texture[y * width_channel];
		unsigned char * texture_x = &texture[y * width_channel];
		*temp_x++ = ypr = *in_x++; 
		*temp_x++ = ypg = *in_x++; 
		*temp_x++ = ypb = *in_x++;
		unsigned char tpr = *texture_x++; 
		unsigned char tpg = *texture_x++;
		unsigned char tpb = *texture_x++;

		float * temp_factor_x = &temp_factor[y * width];
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
			float weight = range_table[range_dist];
			float alpha_ = weight*alpha;
			*temp_x++ = ycr = inv_alpha_*(*in_x++) + alpha_*ypr; 
			*temp_x++ = ycg = inv_alpha_*(*in_x++) + alpha_*ypg; 
			*temp_x++ = ycb = inv_alpha_*(*in_x++) + alpha_*ypb;
			tpr = tcr; tpg = tcg; tpb = tcb;
			ypr = ycr; ypg = ycg; ypb = ycb;
			*temp_factor_x++ = fc = inv_alpha_ + alpha_*fp;
			fp = fc;
		}
		*--temp_x; *temp_x = 0.5f*((*temp_x) + (*--in_x));
		*--temp_x; *temp_x = 0.5f*((*temp_x) + (*--in_x));
		*--temp_x; *temp_x = 0.5f*((*temp_x) + (*--in_x));
		tpr = *--texture_x; 
		tpg = *--texture_x; 
		tpb = *--texture_x;
		ypr = *in_x; ypg = *in_x; ypb = *in_x;

		*--temp_factor_x; *temp_factor_x = 0.5f*((*temp_factor_x) + 1);
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
			float weight = range_table[range_dist];
			float alpha_ = weight * alpha;

			ycr = inv_alpha_ * (*--in_x) + alpha_ * ypr; 
			ycg = inv_alpha_ * (*--in_x) + alpha_ * ypg; 
			ycb = inv_alpha_ * (*--in_x) + alpha_ * ypb;
			*--temp_x; *temp_x = 0.5f*((*temp_x) + ycr);
			*--temp_x; *temp_x = 0.5f*((*temp_x) + ycg);
			*--temp_x; *temp_x = 0.5f*((*temp_x) + ycb);
			tpr = tcr; tpg = tcg; tpb = tcb;
			ypr = ycr; ypg = ycg; ypb = ycb;

			fc = inv_alpha_ + alpha_*fp;
			*--temp_factor_x; 
			*temp_factor_x = 0.5f*((*temp_factor_x) + fc);
			fp = fc;
		}
	}
	alpha = static_cast<float>(exp(-sqrt(2.0) / (sigma_spatial * height)));
	inv_alpha_ = 1 - alpha;
	float * ycy, * ypy, * xcy;
	unsigned char * tcy, * tpy;
	memcpy(out, temp, sizeof(float)* width_channel);

	float * in_factor = temp_factor;
	float*ycf, *ypf, *xcf;
	memcpy(factor_a, in_factor, sizeof(float) * width);
	for (int y = 1; y < height; y++)
	{
		tpy = &texture[(y - 1) * width_channel];
		tcy = &texture[y * width_channel];
		xcy = &temp[y * width_channel];
		ypy = &out[(y - 1) * width_channel];
		ycy = &out[y * width_channel];

		xcf = &in_factor[y * width];
		ypf = &factor_a[(y - 1) * width];
		ycf = &factor_a[y * width];
		for (int x = 0; x < width; x++)
		{
			unsigned char dr = abs((*tcy++) - (*tpy++));
			unsigned char dg = abs((*tcy++) - (*tpy++));
			unsigned char db = abs((*tcy++) - (*tpy++));
			int range_dist = (((dr << 1) + dg + db) >> 2);
			float weight = range_table[range_dist];
			float alpha_ = weight*alpha;
			for (int c = 0; c < channel; c++) 
				*ycy++ = inv_alpha_*(*xcy++) + alpha_*(*ypy++);
			*ycf++ = inv_alpha_*(*xcf++) + alpha_*(*ypf++);
		}
	}
	int h1 = height - 1;
	ycf = factor_b;
	ypf = factor_c;
	memcpy(ypf, &in_factor[h1 * width], sizeof(float) * width);
	for (int x = 0; x < width; x++) 
		factor_a[h1 * width + x] = 0.5f*(factor_a[h1 * width + x] + ypf[x]);

	ycy = factor_cb;
	ypy = factor_cc;
	memcpy(ypy, &temp[h1 * width_channel], sizeof(float)* width_channel);
	int k = 0; 
	for (int x = 0; x < width; x++) {
		for (int c = 0; c < channel; c++) {
			int idx = (h1 * width + x) * channel + c;
			out[idx] = 0.5f*(out[idx] + ypy[k++]) / factor_a[h1 * width + x];
		}
	}

	for (int y = h1 - 1; y >= 0; y--)
	{
		tpy = &texture[(y + 1) * width_channel];
		tcy = &texture[y * width_channel];
		xcy = &temp[y * width_channel];
		float*ycy_ = ycy;
		float*ypy_ = ypy;
		float*out_ = &out[y * width_channel];

		xcf = &in_factor[y * width];
		float*ycf_ = ycf;
		float*ypf_ = ypf;
		float*factor_ = &factor_a[y * width];
		for (int x = 0; x < width; x++)
		{
			unsigned char dr = abs((*tcy++) - (*tpy++));
			unsigned char dg = abs((*tcy++) - (*tpy++));
			unsigned char db = abs((*tcy++) - (*tpy++));
			int range_dist = (((dr << 1) + dg + db) >> 2);
			float weight = range_table[range_dist];
			float alpha_ = weight*alpha;

			float fcc = inv_alpha_*(*xcf++) + alpha_*(*ypf_++);
			*ycf_++ = fcc;
			*factor_ = 0.5f * (*factor_ + fcc);

			for (int c = 0; c < channel; c++)
			{
				float ycc = inv_alpha_*(*xcy++) + alpha_*(*ypy_++);
				*ycy_++ = ycc;
				*out_ = 0.5f * (*out_ + ycc) / (*factor_);
				*out_++;
			}
			*factor_++;
		}
		memcpy(ypy, ycy, sizeof(float) * width_channel);
		memcpy(ypf, ycf, sizeof(float) * width);
	}


	unsigned char * out_img = new unsigned char[width_height_channel];
	for (int i = 0; i < width_height_channel; ++i)
		out_img[i] = static_cast<unsigned char>(out[i]);

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