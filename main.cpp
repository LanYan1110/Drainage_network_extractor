#include <iostream>
#include <fstream>
#include <queue>
#include <cassert>
#include "gdal_priv.h"
#include "cpl_conv.h"
#include <algorithm>

// Storage and access of a raster of a given size
struct Raster {
  std::vector<int> pixels; // where everything is stored
  int max_x, max_y; // number of columns and rows
  
  // Initialise a raster with x columns and y rows
  Raster(int x, int y) {
    max_x = x;
    max_y = y;
    unsigned int total_pixels = x*y;
    pixels.reserve(total_pixels);
  }

    void output_direction(int& current_line, unsigned int* line)
    {
        for (int i = 0; i < max_x; ++i) line[i] =pixels[i + current_line*max_y];
    }


    // Fill values of an entire row
  void add_scanline(const int *line) {
    for (int i = 0; i < max_x -1; ++i) pixels.push_back(line[i]);
  }
  
  // Fill entire raster with zeros
  void fill() {
    unsigned int total_pixels = max_x*max_y;
    for (int i = 0; i < total_pixels; ++i) pixels.push_back(0);
  }
  
  // Access the value of a raster cell to read or write it
  int &operator()(int x, int y) {
    assert(x >= 0 && x < max_x);
    assert(y >= 0 && y < max_y);
    return pixels[x + y*max_x];
  }
  
  // Access the value of a raster cell to read it
  int operator()(int x, int y) const {
    assert(x >= 0 && x < max_x);
    assert(y >= 0 && y < max_y);
    return pixels[x + y*max_x];
  }
};


// A structure that links to a single cell in a Raster
struct RasterCell {
  int x, y; // row and column of the cell
  int elevation;
  int insertion_order;
  
  // Defines a new link to a cell
  RasterCell(int x, int y, int elevation, int insertion_order) {
    this->x = x;
    this->y = y;
    this->elevation = elevation;
    this->insertion_order = insertion_order;
  }
  
  // Define the order of the linked cells (to be used in a priority_queue)
  bool operator<(const RasterCell &other) const {
      if (this->elevation > other.elevation) {
          return true;
      }
      else if (this->elevation == other.elevation) {
          if (this->insertion_order < other.elevation)
              return true;
          else
              return false;
      }
      else
        return false;
  }
};

// Write the values in a linked raster cell (useful for debugging)
std::ostream& operator<<(std::ostream& os, const RasterCell& c) {
  os << "{h=" << c.elevation << ", o=" << c.insertion_order << ", x=" << c.x << ", y=" << c.y << "}";
  return os;
}

int insertion_order(){
        static int counter = 0;
        return ++counter;
}

void print_queue(std::priority_queue<RasterCell, std::deque<RasterCell>> q) {
    while (!q.empty()) {
        std::cout << q.top() << ' ';
        q.pop();
    }
    std::cout << '\n';
}

void process_points(int x, int y, std::priority_queue<RasterCell, std::deque<RasterCell>> priority_queue, Raster input_raster, Raster flow_direction, int scenario)
{
    if (flow_direction(x,y) ==0) {
        auto new_cell = RasterCell(x, y, input_raster(x,y), insertion_order());
        priority_queue.push(new_cell);
        flow_direction(x, y) = scenario;
    }
}

int main(int argc, const char * argv[]) {
  // Open dataset
  GDALDataset  *input_dataset;
  GDALAllRegister();
  input_dataset = (GDALDataset *)GDALOpen("S33W070.hgt", GA_ReadOnly);
  if (input_dataset == NULL) {
    std::cerr << "Couldn't open file" << std::endl;
    return 1;
  }
  
  // Print dataset info
  double geo_transform[6];
  std::cout << "Driver: " << input_dataset->GetDriver()->GetDescription() << "/" << input_dataset->GetDriver()->GetMetadataItem(GDAL_DMD_LONGNAME) << std::endl;;
  std::cout << "Size is " << input_dataset->GetRasterXSize() << "x" << input_dataset->GetRasterYSize() << "x" << input_dataset->GetRasterCount() << std::endl;
  if (input_dataset->GetProjectionRef() != NULL) std::cout << "Projection is '" << input_dataset->GetProjectionRef() << "'" << std::endl;
  if (input_dataset->GetGeoTransform(geo_transform) == CE_None) {
    std::cout << "Origin = (" << geo_transform[0] << ", " << geo_transform[3] << ")" << std::endl;
    std::cout << "Pixel Size = (" << geo_transform[1] << ", " << geo_transform[5] << ")" << std::endl;
  }
  
  // Print Band 1 info
  GDALRasterBand *input_band;
  int nBlockXSize, nBlockYSize;
  int bGotMin, bGotMax;
  double adfMinMax[2];
  input_band = input_dataset->GetRasterBand(1);
  input_band->GetBlockSize(&nBlockXSize, &nBlockYSize);
  std::cout << "Band 1 Block=" << nBlockXSize << "x" << nBlockYSize << " Type=" << GDALGetDataTypeName(input_band->GetRasterDataType()) << " ColorInterp=" << GDALGetColorInterpretationName(input_band->GetColorInterpretation()) << std::endl;
  adfMinMax[0] = input_band->GetMinimum(&bGotMin);
  adfMinMax[1] = input_band->GetMaximum(&bGotMax);
  if (!(bGotMin && bGotMax)) GDALComputeRasterMinMax((GDALRasterBandH)input_band, TRUE, adfMinMax);
  std::cout << "Min=" << adfMinMax[0] << " Max=" << adfMinMax[1] << std::endl;
  
  // Read Band 1 line by line
  int nXSize = input_band->GetXSize();
  int nYSize = input_band->GetYSize();
  Raster input_raster(nXSize, nYSize);
  for (int current_scanline = 0; current_scanline < nYSize; ++current_scanline) {
    int *scanline = (int *)CPLMalloc(sizeof(float)*nXSize);
    if (input_band->RasterIO(GF_Read, 0, current_scanline, nXSize, 1,
                         scanline, nXSize, 1, GDT_Int32,
                         0, 0) != CPLE_None) {
      std::cerr << "Couldn't read scanline " << current_scanline << std::endl;
      return 1;
    } input_raster.add_scanline(scanline);
    CPLFree(scanline);
  } std::cout << "Created raster: " << input_raster.max_x << "x" << input_raster.pixels.size()/input_raster.max_y << " = " << input_raster.pixels.size() << std::endl;
  
  // Flow direction
  Raster flow_direction(input_raster.max_x, input_raster.max_y);
    flow_direction.fill();
    std::priority_queue<RasterCell, std::deque<RasterCell>> priority_queue;
    for (int x=0; x < nXSize -1 ; x++){
        for (int y=0; y < nYSize -1 ; y++){
            if (x == 0 || x == nXSize -1 || y == 0 || y == nYSize -1){
                auto priority_cell = RasterCell(x, y, input_raster(x,y), insertion_order());
                priority_queue.push(priority_cell);
            }
        }

    }
    for (; !priority_queue.empty(); priority_queue.pop()){
        RasterCell c = priority_queue.top();
        if (c.x == 0 and c.y ==0 ){
            int x = c.x +1; int y = c.y +1; int scenario = 25;
            process_points(x, y, priority_queue, input_raster, flow_direction, scenario);
        }

        else if ( c.x ==0 and c.y == nYSize ){
            int x = c.x +1; int y = c.y -1; int scenario = 15;
            process_points(x, y, priority_queue, input_raster, flow_direction, scenario);
        }

        else if (c.y ==0 and c.x == nXSize){
            int x = c.x -1; int y = c.y +1; int scenario = 35;
            process_points(x, y, priority_queue, input_raster, flow_direction, scenario);
        }

        else if ( c.y == nYSize and c.x == nXSize){
            int x = c.x -1; int y = c.y -1; int scenario = 5;
            process_points(x, y, priority_queue, input_raster, flow_direction, scenario);
        }

        else if (c.y ==0 and c.x != 0 and c.x != nXSize){
            int x = c.x ; int y = c.y +1; int scenario = 20;
            process_points(x, y, priority_queue, input_raster, flow_direction, scenario);
        }

        else if (c.x ==0 and c.y != 0 and c.y != nYSize){
            int x = c.x +1; int y = c.y; int scenario = 30;
            process_points(x, y, priority_queue, input_raster, flow_direction, scenario);
        }
        else if (c.y == nYSize and c.x != 0 and c.x != nXSize){
            int x = c.x; int y = c.y -1; int scenario = 10;
            process_points(x, y, priority_queue, input_raster, flow_direction, scenario);
        }
        else if (c.x == nXSize and c.y != 0 and c.y != nXSize){
            int x = c.x -1; int y = c.y; int scenario = 40;
            process_points(x, y, priority_queue, input_raster, flow_direction, scenario);
        }

        else
        {
            int x_1 = c.x +1; int y_1 = c.y; int scenario_1 = 20;
            process_points(x_1, y_1, priority_queue, input_raster, flow_direction, scenario_1);
            int x_2 = c.x +1; int y_2 = c.y-1; int scenario_2 = 15;
            process_points(x_2, y_2, priority_queue, input_raster, flow_direction, scenario_2);
            int x_3 = c.x; int y_3 = c.y-1; int scenario_3 = 10;
            process_points(x_1, y_1, priority_queue, input_raster, flow_direction, scenario_3);
            int x_4 = c.x -1; int y_4 = c.y-1; int scenario_4 = 5;
            process_points(x_1, y_1, priority_queue, input_raster, flow_direction, scenario_4);
            int x_5 = c.x -1; int y_5 = c.y; int scenario_5 = 40;
            process_points(x_1, y_1, priority_queue, input_raster, flow_direction, scenario_5);
            int x_6 = c.x +1; int y_6 = c.y-1; int scenario_6 = 35;
            process_points(x_1, y_1, priority_queue, input_raster, flow_direction, scenario_6);
            int x_7 = c.x; int y_7 = c.y+1; int scenario_7 = 30;
            process_points(x_1, y_1, priority_queue, input_raster, flow_direction, scenario_7);
            int x_8 = c.x +1; int y_8 = c.y+1; int scenario_8 = 25;
            process_points(x_1, y_1, priority_queue, input_raster, flow_direction, scenario_8);
        }
    }

//   Write flow direction
    std::string tiffname("test.tif");
     // set the name and path of the flow_direction file
    GDALDriver* driverGeotiff(GetGDALDriverManager()->GetDriverByName("GTiff"));
    GDALDataset* geotiffDataset(driverGeotiff->Create(tiffname.c_str(), nXSize, nYSize, 1, GDT_Int32, NULL));

    geotiffDataset->SetGeoTransform(geo_transform);
    geotiffDataset->SetProjection(input_dataset->GetProjectionRef());

    auto* output_line((unsigned int*)CPLMalloc(sizeof(unsigned int)* nXSize));

    for (int current_scanline = 0; current_scanline != nYSize; ++current_scanline)
    {
        input_raster.output_direction(current_scanline, output_line); // similar to add_scanline line[i] = direction/accumulation of each cell

        if (geotiffDataset->GetRasterBand(1)->RasterIO(GF_Write, 0, current_scanline, nXSize, 1,
                                                       output_line, nXSize, 1, GDT_UInt32, 0, 0) != CPLE_None) // write the values into the raster file
        {
            std::cerr << "Couldn't load output_line " << current_scanline << '\n';
            return 1;
        }

    }

    CPLFree(output_line);

  // to do
  
  // Flow accumulation
//  Raster flow_accumulation(input_raster.max_x, input_raster.max_y);
  // to do
  
  // Write flow accumulation
  // to do

  // Close input dataset
  GDALClose(input_dataset);

  return 0;
}
