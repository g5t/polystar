
#include <cassert>
#include <fstream>
#include "bitmap.hpp"

using namespace polystar::bitmap;

const size_t BYTES_PER_PIXEL = 3; /// red, green, & blue
const size_t FILE_HEADER_SIZE = 14;
const size_t INFO_HEADER_SIZE = 40;

void write_header(std::ofstream & of, size_t height, size_t stride){
  auto fileSize = FILE_HEADER_SIZE + INFO_HEADER_SIZE + (stride * height);

  static unsigned char fileHeader[] = {
    0, 0,     /// signature
    0, 0, 0, 0, /// image file size in bytes
    0, 0, 0, 0, /// reserved
    0, 0, 0, 0, /// start of pixel array
  };

  fileHeader[0] = (unsigned char) ('B');
  fileHeader[1] = (unsigned char) ('M');
  fileHeader[2] = (unsigned char) (fileSize);
  fileHeader[3] = (unsigned char) (fileSize >> 8);
  fileHeader[4] = (unsigned char) (fileSize >> 16);
  fileHeader[5] = (unsigned char) (fileSize >> 24);
  fileHeader[10] = (unsigned char) (FILE_HEADER_SIZE + INFO_HEADER_SIZE);

  of.write((char *)(fileHeader), FILE_HEADER_SIZE);
}
void write_info_header(std::ofstream & of, size_t height, size_t width){
  static unsigned char infoHeader[] = {
    0, 0, 0, 0, /// header size
    0, 0, 0, 0, /// image width
    0, 0, 0, 0, /// image height
    0, 0,     /// number of color planes
    0, 0,     /// bits per pixel
    0, 0, 0, 0, /// compression
    0, 0, 0, 0, /// image size
    0, 0, 0, 0, /// horizontal resolution
    0, 0, 0, 0, /// vertical resolution
    0, 0, 0, 0, /// colors in color table
    0, 0, 0, 0, /// important color count
  };

  infoHeader[0] = (unsigned char) (INFO_HEADER_SIZE);
  infoHeader[4] = (unsigned char) (width);
  infoHeader[5] = (unsigned char) (width >> 8);
  infoHeader[6] = (unsigned char) (width >> 16);
  infoHeader[7] = (unsigned char) (width >> 24);
  infoHeader[8] = (unsigned char) (height);
  infoHeader[9] = (unsigned char) (height >> 8);
  infoHeader[10] = (unsigned char) (height >> 16);
  infoHeader[11] = (unsigned char) (height >> 24);
  infoHeader[12] = (unsigned char) (1);
  infoHeader[14] = (unsigned char) (BYTES_PER_PIXEL * 8);

  of.write((char *)(infoHeader), INFO_HEADER_SIZE);
}

void polystar::bitmap::write(const std::vector<std::vector<Color>> & image, const std::string & filename){
  auto height = image.size();
  auto width = image[0].size();

  auto width_bytes = width * BYTES_PER_PIXEL;
  unsigned char padding[3] = {0, 0, 0};
  auto paddingSize = (4 - (width_bytes) % 4) % 4;
  auto stride = width_bytes + paddingSize;

  auto of = std::ofstream(filename, std::ios::binary);

  write_header(of, height, stride);
  write_info_header(of, height, width);

  for (const auto & row: image){
    for (const auto & val: row) val.write(of);
    of.write((char *)padding, static_cast<std::streamsize>(paddingSize));
  }
}
