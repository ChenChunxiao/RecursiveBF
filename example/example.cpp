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
		float sigma_spatial = static_cast<float>(atof(argv[3]));
		float sigma_range = static_cast<float>(atof(argv[4]));

		int width, height, channel;
		unsigned char * img = stbi_load(filename_in, &width, &height, &channel, 0);

		Timer timer;
		timer.start();
		recursive_bf(img, sigma_spatial, sigma_range, width, height, channel);
		printf("Elapsed time = %2.5f secs", timer.elapsedTime());

		stbi_write_bmp(filename_out, width, height, channel, img);
		delete[] img;
	}
}
