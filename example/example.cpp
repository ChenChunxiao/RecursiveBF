#include "../include/rbf.hpp"
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image.hpp"
#include "stb_image_write.hpp"
#include <stdio.h>
#include <time.h>
#define LEN_MAX	256


class Timer {
private:
	unsigned long begTime;
public:
	void start() { begTime = clock(); }
	float elapsedTime() { return float((unsigned long)clock() - begTime) / CLOCKS_PER_SEC; }
};

int main(int argc, char*argv[])
{
	if (argc != 5)
	{
		printf("Usage:\n");
		printf("--------------------------------------------------------------------\n\n");
		printf("*.exe: filename_out filename_in (only ppm image) \n");
		printf("       sigma_spatial(e.g., 0.03) sigma_range(e.g., 0.1)\n\n");
		printf("--------------------------------------------------------------------\n");
		return(-1);
	}
	else
	{
		const char * filename_out = argv[1];
		const char * filename_in = argv[2];
		double sigma_spatial = atof(argv[3]);
		double sigma_range = atof(argv[4]);

		int width, height, channel;
		unsigned char * img_in = stbi_load(filename_in, &width, &height, &channel, 0);

		Timer timer;
		timer.start();
		unsigned char * img_out = recursive_bf(
			img_in, sigma_spatial, sigma_range, width, height, channel);
		printf("Elapsed time = %2.5f secs", timer.elapsedTime());
		stbi_write_bmp(filename_out, width, height, channel, img_out);
		delete[] img_in;
		delete[] img_out;
	}
}
