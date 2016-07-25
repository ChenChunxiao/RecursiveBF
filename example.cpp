#include "recursive_bf.hpp"
#include <time.h>
#define LEN_MAX	256

int example_filtering(int argc, char*argv[]);
int example_detail_enhancement(int argc, char*argv[]);
int main(int argc, char*argv[])
{
	if (argc != 6)
	{
		printf("Usage:\n");
		printf("--------------------------------------------------------------------\n\n");
		printf("*.exe: filename_out filename_in (only ppm image) \n");
		printf("       sigma_spatial(e.g., 0.03) sigma_range(e.g., 0.1)\n");
		printf("       application_id (0 for GDBF; 1 for RBF)\n\n");
		printf("--------------------------------------------------------------------\n");
		return(-1);
	}
	else
	{
		return example_filtering(argc, argv);
	}
}

template <typename T>
inline T *** qx_alloc_3(int n, int r, int c, int padding = 10)
{
	int rc = r*c;
	int i, j;
	T *a = (T*)malloc(sizeof(T)*(n*rc + padding));
	if (a == NULL) { 
		printf("qx_allocu_3() fail, Memory is too huge, fail.\n"); 
		getchar(); 
		exit(0); 
	}
	T **p = (T**)malloc(sizeof(T*)*n*r);
	T ***pp = (T***)malloc(sizeof(T**)*n);
	for (i = 0; i<n; i++)
		for (j = 0; j<r; j++)
			p[i*r + j] = &a[i*rc + j*c];
	for (i = 0; i<n; i++)
		pp[i] = &p[i*r];
	return pp;
}
template <typename T>
inline T ** qx_alloc(int r, int c, int padding = 10)
{
	T *a = (T*)malloc(sizeof(T)*(r*c + padding));
	if (a == NULL) { 
		printf("qx_allocd() fail, Memory is too huge, fail.\n"); 
		getchar(); 
		exit(0); 
	}
	T **p = (T**)malloc(sizeof(T*)*r);
	for (int i = 0; i<r; i++) p[i] = &a[i*c];
	return p;
}

template <typename T>
inline void qx_free_3(T ***p)
{
	if (p != NULL)
	{
		free(p[0][0]);
		free(p[0]);
		free(p);
		p = NULL;
	}
}
template <typename T>
inline void qx_free(T **p)
{
	if (p != NULL)
	{
		free(p[0]);
		free(p);
		p = NULL;
	}
}

void qx_image_size(char* file_name, int &h, int &w, int *nr_channel = NULL)
{
	FILE * file_in; int nrc;
	char line[LEN_MAX];
	int	i;
	unsigned char *image = NULL;
	fopen_s(&file_in, file_name, "rb");
	//fopen_s(&file_in,file_name,"rb");
	if (!file_in)
	{
		printf("Please check input file_name: %s\n", file_name);
		getchar();
		exit(0);
	}
	if (fgetc(file_in) == 'P')
		//fscanf_s(file_in,"%d\n",&i);
		fscanf_s(file_in, "%d\n", &i);
	else
	{
		printf("Bad	header in ppm file.\n");
		getchar();
		exit(1);
	}
	while (fgets(line, LEN_MAX, file_in) != NULL)
	{
		if (line[0] == '#') continue;
		else
		{
			//sscanf_s(line, "%d %d\n",	&w,	&h);
			sscanf_s(line, "%d %d\n", &w, &h);
			if (i == 9 && nr_channel != NULL)
			{
				fscanf_s(file_in, "%d\n", &nrc);
				(*nr_channel) = nrc;
			}
			break;
		}
	}
}

int qx_loadimage(char* filename, unsigned char *image, int h, int w, int *nr_channel = NULL)
{
	FILE * file_in; int nrc;
	char line[LEN_MAX];
	int	i; int imax, hc, wc;
	unsigned char *image_ = image;
	//fopen_s(&file_in,filename,"rb");
	fopen_s(&file_in, filename, "rb");

	if (!file_in)
	{
		printf("Please check input filename: %s\n", filename);
		exit(0);
	}
	if (fgetc(file_in) == 'P')
		fscanf_s(file_in, "%d\n", &i);
	else
	{
		printf("Bad	header in ppm file.\n");
		exit(1);
	}
	while (fgets(line, LEN_MAX, file_in) != NULL)
	{
		if (line[0] == '#') continue;
		else
		{
			sscanf_s(line, "%d %d\n", &wc, &hc);
			break;
		}
	}
	char str_tmp[100];
	switch (i)
	{
	case 5:
		fgets(str_tmp, 100, file_in);
		imax = atoi(str_tmp);
		if (nr_channel != NULL) (*nr_channel) = 1;
		memset(image, 0, sizeof(unsigned char)*h*w);
		fread(image, sizeof(unsigned char), h*w, file_in);
		break;
	case 6:
		fgets(str_tmp, 100, file_in);
		imax = atoi(str_tmp);
		if (nr_channel != NULL) (*nr_channel) = 3;
		memset(image, 0, sizeof(unsigned char)*h*w * 3);
		fread(image, sizeof(unsigned char), h*w * 3, file_in);
		break;
	case 2:
		fgets(str_tmp, 100, file_in);
		imax = atoi(str_tmp);
		for (int y = 0; y<h; y++) for (int x = 0; x<w; x++)
		{
			//if(fscanf_s(file_in,"%d",&imax)!=1){printf("error in reading file.\n");getchar();exit(0);}
			fscanf_s(file_in, "%d", &imax);
			*image_++ = imax;
		}
		break;
	case 3:
		fgets(str_tmp, 100, file_in);
		imax = atoi(str_tmp);
		int cr, cg, cb;
		for (int y = 0; y<h; y++) for (int x = 0; x<w; x++)
		{
			//if(fscanf_s(file_in,"%d%d%d",&cr,&cg,&cb)!=3){printf("error in reading file.\n");getchar();exit(0);}
			fscanf_s(file_in, "%d%d%d", &cr, &cg, &cb);
			*image_++ = cr; *image_++ = cg; *image_++ = cb;
		}
		break;
	case 9:
		fgets(str_tmp, 100, file_in);
		nrc = atoi(str_tmp);
		fscanf_s(file_in, "%d\n", &nrc);
		if (nr_channel != NULL) (*nr_channel) = nrc;
		fgets(str_tmp, 100, file_in);
		imax = atoi(str_tmp);
		fread(image, sizeof(unsigned char), h*w*nrc, file_in);
		break;
	default:
		printf("Can not open image [%s]!!\n", filename);
		break;
	}
	fclose(file_in);
	return (0);
}

void qx_saveimage(char* filename, unsigned char *image, int h, int w, int channel)
{
	FILE* file_out; unsigned char maxx = 255;
	//fopen_s(&file_out,filename,"wb");
	fopen_s(&file_out, filename, "wb");

	if (channel == 1) fprintf(file_out, "P5\n%d %d\n%d\n", w, h, maxx);
	//else if(channel==3) fprintf(file_out,"P6\n%d %d\n%d\n",w,h,maxx);
	else if (channel == 3) fprintf(file_out, "P6\n%d %d\n%d\n", w, h, maxx);
	else fprintf(file_out, "P9\n%d %d\n%d\n%d\n", w, h, channel, maxx);
	fwrite(image, sizeof(unsigned char), h*w*channel, file_out);

	fclose(file_out);
}

class Timer {
private:
	unsigned long begTime;
public:
	void start() { begTime = clock(); }
	float elapsedTime() { return float((unsigned long)clock() - begTime) / CLOCKS_PER_SEC; }
};


int example_filtering(int argc, char*argv[])
{
	char*filename_out = argv[1];
	char*filename_in = argv[2];
	double sigma_spatial = atof(argv[3]);
	double sigma_range = atof(argv[4]);
	int filter_type = atoi(argv[5]);

	int h, w;//height and width of the input image
	qx_image_size(filename_in, h, w);//obtain the height and width of the input image
	unsigned char ***texture = qx_alloc_3<unsigned char>(h, w, 3);//allocate memory
	double ***image = qx_alloc_3<double>(h, w, 3);
	double ***image_filtered = qx_alloc_3<double>(h, w, 3);
	double ***temp = qx_alloc_3<double>(h, w, 3);
	double ***temp_2 = qx_alloc_3<double>(2, w, 3);

	qx_loadimage(filename_in, texture[0][0], h, w);//read an color image
	for (int y = 0; y<h; y++) for (int x = 0; x<w; x++) for (int c = 0; c<3; c++) image[y][x][c] = texture[y][x][c];//initialize the original color image

	Timer timer;
	timer.start();
	if (filter_type == 0)//GDBF
	{
		qx_gradient_domain_recursive_bilateral_filter(
			image_filtered, image, texture, 
			sigma_spatial, sigma_range, h, w, temp, temp_2);
	}
	else//RBF
	{
		double**temp_factor = qx_alloc<double>(h * 2 + 2, w);
		qx_recursive_bilateral_filter(
			image_filtered, image, texture, 
			sigma_spatial, sigma_range, h, w, temp, temp_2, 
			temp_factor, &(temp_factor[h]), &(temp_factor[h + h]));
		qx_free(temp_factor); 
		temp_factor = NULL;
	}
	printf("Elapsed time: %2.5f", timer.elapsedTime());

	for (int y = 0; y<h; y++) 
		for (int x = 0; x<w; x++) 
			for (int c = 0; c<3; c++) 
				texture[y][x][c] = (unsigned char)image_filtered[y][x][c];
	qx_saveimage(filename_out, texture[0][0], h, w, 3);

	qx_free_3(texture);
	qx_free_3(image);
	qx_free_3(image_filtered);
	qx_free_3(temp);
	qx_free_3(temp_2);
	return(0);
}