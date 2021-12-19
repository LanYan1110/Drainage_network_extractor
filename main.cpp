#include <iostream>
#include <fstream>
#include <queue>
#include <cassert>
#include <stack>
#include "gdal_priv.h"
#include "cpl_conv.h"

int insert_order = 0;
// Storage and access of a raster of a given size
struct Raster {
    std::vector<unsigned int> pixels; // where everything is stored
    int max_x, max_y; // number of columns and rows

    // Initialise a raster with x columns and y rows
    Raster(int x, int y) {
        max_x = x;
        max_y = y;
        unsigned int total_pixels = x * y;
        pixels.reserve(total_pixels);
    }

    // Fill values of an entire row
    void add_scanline(const int* line) {
        for (int i = 0; i < max_x; ++i) pixels.push_back(line[i]);
    }

    // Fill entire raster with zeros
    void fill() {
        unsigned int total_pixels = max_x * max_y;
        for (int i = 0; i < total_pixels; ++i) pixels.push_back(0);
    }

    void fill_1() {
        unsigned int total_pixels = max_x * max_y;
        for (int i = 0; i < total_pixels; ++i) pixels.push_back(1);
    }

    // Access the value of a raster cell to read or write it
    unsigned int& operator()(int x, int y) {
        assert(x >= 0 && x < max_x);
        assert(y >= 0 && y < max_y);
        return pixels[x + y * max_x];
    }

    // Access the value of a raster cell to read it
    unsigned int operator()(int x, int y) const {
        assert(x >= 0 && x < max_x);
        assert(y >= 0 && y < max_y);
        return pixels[x + y * max_x];
    }
};

// A structure that links to a single cell in a Raster
struct RasterCell {
    int x, y; // row and column of the cell
    unsigned int elevation;
    int insertion_order;
    // Defines a new link to a cell
    RasterCell(int x, int y, unsigned int elevation, int insertion_order) {
        this->x = x;
        this->y = y;
        this->elevation = elevation;
        this->insertion_order = insertion_order;
    }
    // Define the order of the linked cells (to be used in a priority_queue)
    bool operator<(const RasterCell& other) const {
        // to do with statements like if (this->elevation > other.elevation) return false/true;
        if (this->elevation > other.elevation) {
            return true;
        }
        else if (this->elevation == other.elevation){
            //(grid cells added earlier have higher precedence in case of equal elevation)
            if (this->insertion_order < other.insertion_order)
                return false;
            else
                return true;
        }
        else {
            return false;
        }
    }
};

// Write the values in a linked raster cell (useful for debugging)
std::ostream& operator<<(std::ostream& os, const RasterCell& c) {
    os << "{h=" << c.elevation << ", o=" << c.insertion_order << ", x=" << c.x << ", y=" << c.y << "}";
    return os;
}

//process neighbours that are not yet in the search list and not yet processed
void neighbour_processing(Raster& input_raster, int x, int y, Raster& flow_direction, int fd_value, std::priority_queue<RasterCell, std::deque<RasterCell>>& cells_to_process_flow, int nXSize, int nYSize)
{
    if (x >= 0 and x < nXSize and y >= 0 and y < nYSize and flow_direction(x, y) == 0){
        //add the neighbour into the prority queue
        auto height = input_raster(x,y);
        RasterCell n(x, y, height, insert_order++);
        cells_to_process_flow.push(n);
        //write the cell's direction in the flow direction raster
        flow_direction(x, y) = fd_value;
    }
    else if (x >= 0 and x < nXSize and y >= 0 and y < nYSize and flow_direction(x, y) == -1 ) {
        flow_direction(x, y) = fd_value;
    }
}

int main(int argc, const char* argv[]) {
    // Open dataset
    GDALDataset* input_dataset;
    GDALAllRegister();
    input_dataset = (GDALDataset *)GDALOpen("N74E018.hgt", GA_ReadOnly);
    if (input_dataset == NULL)
    {
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
    GDALRasterBand* input_band;
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
        int* scanline = (int*)CPLMalloc(sizeof(float) * nXSize);
        if (input_band->RasterIO(GF_Read, 0, current_scanline, nXSize, 1,
                                 scanline, nXSize, 1, GDT_Int32,
                                 0, 0) != CPLE_None) {
            std::cerr << "Couldn't read scanline " << current_scanline << std::endl;
            return 1;
        } input_raster.add_scanline(scanline);
        CPLFree(scanline);
    } std::cout << "Created raster: " << input_raster.max_x << "x" << input_raster.pixels.size() / input_raster.max_y << " = " << input_raster.pixels.size() << std::endl;

    // Flow direction
    Raster flow_direction(input_raster.max_x, input_raster.max_y);
    flow_direction.fill();
    std::priority_queue<RasterCell, std::deque<RasterCell>> cells_to_process_flow;
    std::stack<RasterCell> accumulation_stack;

    // Write flow direction
    //push boundary cells into priority queue (initial priority queue)
    for (int j = 0; j < nYSize; j++)
    {
        for (int i = 0; i < nXSize; i++)
        {
            if (i == 0 || j == 0 || i == nXSize - 1 || j == nYSize - 1)
            {
                RasterCell c(i, j, input_raster(i, j), insert_order++);
                cells_to_process_flow.push(c);
                flow_direction(i, j) = -1;
            }
        }
    }

    //now we have our initial priority queue, we begin to process it
    for (; !cells_to_process_flow.empty(); cells_to_process_flow.pop())
    {
        RasterCell c = cells_to_process_flow.top();
        accumulation_stack.push(c);
        neighbour_processing(input_raster, c.x + 1, c.y - 1, flow_direction, 25, cells_to_process_flow, nXSize, nYSize); //1
        neighbour_processing(input_raster, c.x, c.y - 1, flow_direction, 30, cells_to_process_flow, nXSize, nYSize);//2
        neighbour_processing(input_raster, c.x - 1, c.y - 1, flow_direction, 35, cells_to_process_flow, nXSize, nYSize);//3
        neighbour_processing(input_raster, c.x + 1, c.y, flow_direction, 20, cells_to_process_flow, nXSize, nYSize); //4
        neighbour_processing(input_raster, c.x - 1, c.y, flow_direction, 40, cells_to_process_flow, nXSize, nYSize); //5
        neighbour_processing(input_raster, c.x + 1, c.y + 1, flow_direction, 15, cells_to_process_flow, nXSize, nYSize); //6
        neighbour_processing(input_raster, c.x, c.y + 1, flow_direction, 10, cells_to_process_flow, nXSize, nYSize); //7
        neighbour_processing(input_raster, c.x - 1, c.y + 1, flow_direction, 5, cells_to_process_flow, nXSize, nYSize); //8+
    }

    //flow accumulation
    Raster flow_accumulation(input_raster.max_x, input_raster.max_y);
    flow_accumulation.fill();
    for (; !accumulation_stack.empty(); accumulation_stack.pop())
    {
        RasterCell c = accumulation_stack.top();
        //get its neighbours

        if (flow_direction(c.x, c.y) == 25) {
            flow_accumulation(c.x - 1, c.y + 1) = flow_accumulation(c.x - 1, c.y + 1) + flow_accumulation(c.x, c.y) +1;}
        if (flow_direction(c.x, c.y) == 30){
            flow_accumulation(c.x, c.y + 1) = flow_accumulation(c.x, c.y +1) + flow_accumulation(c.x, c.y) +1;}
        if ( flow_direction(c.x, c.y) == 35){
            flow_accumulation(c.x + 1, c.y + 1) =  flow_accumulation(c.x + 1, c.y + 1) + flow_accumulation(c.x, c.y) +1;}
        if (flow_direction(c.x, c.y) == 20){
            flow_accumulation(c.x - 1, c.y) = flow_accumulation(c.x - 1, c.y) + flow_accumulation(c.x, c.y) +1;}
        if (flow_direction(c.x, c.y) == 40){
            flow_accumulation(c.x + 1, c.y) = flow_accumulation(c.x + 1, c.y) + flow_accumulation(c.x, c.y) +1;}
        if (flow_direction(c.x, c.y) == 15){
            flow_accumulation(c.x - 1, c.y - 1) = flow_accumulation(c.x - 1, c.y - 1) + flow_accumulation(c.x, c.y) +1;}
        if (flow_direction(c.x, c.y) == 10){
            flow_accumulation(c.x, c.y - 1) = flow_accumulation(c.x, c.y - 1) + flow_accumulation(c.x, c.y) +1;}
        if (flow_direction(c.x, c.y) == 5){
            flow_accumulation(c.x + 1, c.y - 1) = flow_accumulation(c.x + 1, c.y - 1) + flow_accumulation(c.x, c.y) +1;}

    }

//        unsigned int value = flow_direction(c.x, c.y);
//        unsigned int add = flow_accumulation(c.x, c.y) +1;
//        switch (value){
//            case(5):{
//                auto flow_old = flow_accumulation(c.x -1, c.y -1);
//                flow_accumulation(c.x -1, c.y -1) = flow_old + add;
//                std::cout << "old value  " << flow_old << " addition = " << add << " new_flow: "<< flow_accumulation(c.x -1, c.y -1) << " case =5"<< std::endl;
//                break;
//            }
//            case(10):{
//                auto flow_old = flow_accumulation(c.x, c.y -1);
//                flow_accumulation(c.x, c.y -1) = flow_old + add;
//                std::cout << "old value  " << flow_old << " addition = " << add << " new_flow: "<< flow_accumulation(c.x, c.y -1) << " case =10"<< std::endl;
//                break;
//            }
//            case(15):{
//                auto flow_old = flow_accumulation(c.x +1, c.y-1);
//                flow_accumulation(c.x +1, c.y-1) =  flow_old + add;
//                std::cout << "old value  " << flow_old << " addition = " << add << " new_flow: "<< flow_accumulation(c.x +1, c.y-1) << " case =15"<< std::endl;
//                break;
//            }
//            case(40):{
//                auto flow_old = flow_accumulation(c.x - 1, c.y);
//                flow_accumulation(c.x - 1, c.y) =  flow_old + add;
//                std::cout << "old value  " << flow_old << " addition = " << add << " new_flow: "<< flow_accumulation(c.x - 1, c.y) << " case =40"<< std::endl;
//                break;
//            }
//            case(20):{
//                auto flow_old = flow_accumulation(c.x +1, c.y);
//                flow_accumulation(c.x +1, c.y)  = flow_old + add;
//                std::cout << "old value  " << flow_old << " addition = " << add << " new_flow: "<< flow_accumulation(c.x +1, c.y) << " case =20"<< std::endl;
//                break;
//            }
//            case(35):{
//                auto flow_old = flow_accumulation(c.x -1, c.y +1);
//                flow_accumulation(c.x -1, c.y +1) = flow_old + add;
//                std::cout << "old value  " << flow_old << " addition = " << add << " new_flow: "<< flow_accumulation(c.x -1, c.y +1) << " case =35"<< std::endl;
//                break;
//            }
//            case(30):{
//                auto flow_old = flow_accumulation(c.x , c.y +1);
//                flow_accumulation(c.x , c.y +1) = flow_old + add;
//                std::cout << "old value  " << flow_old << " addition = " << add << " new_flow: "<< flow_accumulation(c.x , c.y +1) << " case =30"<< std::endl;
//                break;
//            }
//            case(25):{
//                auto flow_old = flow_accumulation(c.x +1, c.y+1);
//                flow_accumulation(c.x +1, c.y+1) = flow_old + add;
//                std::cout << "old value  " << flow_old << " addition = " << add << " new_flow: "<< flow_accumulation(c.x +1, c.y+1) << " case =25"<< std::endl;
//                break;
//            }
//            default:{   break;
//            }
//        }
//    }

    GDALDataset* geotiffDataset;
    GDALDataset* geotiffDataset_2;
    GDALDriver* driverGeotiff;
    GDALDriver* driverGeotiff_2;
    GDALRasterBand* geotiffBand;
    GDALRasterBand* geotiffBand_2;
    input_dataset->GetGeoTransform(geo_transform);
    std::string extension(".tif"), tiffname;
    tiffname = (std::string)"flow_direction" + extension;
    std::string extension_2(".tif"), tiffname_2;
    tiffname_2 = (std::string)"flow_accumulation" + extension_2;
    driverGeotiff = GetGDALDriverManager()->GetDriverByName("GTiff");
    driverGeotiff_2 = GetGDALDriverManager()->GetDriverByName("GTiff");
    geotiffDataset = driverGeotiff->Create(tiffname.c_str(), nXSize, nYSize, 1, GDT_Int32, NULL);
    geotiffDataset_2 = driverGeotiff_2->Create(tiffname_2.c_str(), nXSize, nYSize, 1, GDT_Int32 , NULL);
    geotiffDataset->SetGeoTransform(geo_transform);
    geotiffDataset_2->SetGeoTransform(geo_transform);
    geotiffDataset->SetProjection(input_dataset->GetProjectionRef());
    geotiffDataset_2->SetProjection(input_dataset->GetProjectionRef());
    auto* rowBuff = (unsigned int*)CPLMalloc(sizeof(unsigned int) * nYSize);
    auto* rowBuff_2 = (unsigned int*)CPLMalloc(sizeof(unsigned int) * nYSize);
    for (int row = 0; row < nYSize; row++)
    {
        for (int col = 0; col < nXSize; col++)
        {
            rowBuff[col] = flow_direction(col, row);
        }
        geotiffDataset->GetRasterBand(1)->RasterIO(GF_Write, 0, row, nYSize, 1, rowBuff, nYSize, 1, GDT_UInt32, 0, 0);
    }
    for (int row = 0; row < nYSize; row++)
    {
        for (int col = 0; col < nXSize; col++){
            rowBuff_2[col] = flow_accumulation(col, row);
        }
        geotiffDataset_2->GetRasterBand(1)->RasterIO(GF_Write, 0, row, nYSize, 1, rowBuff_2, nYSize, 1, GDT_UInt32, 0, 0);
    }
    GDALClose(input_dataset);
    GDALClose(geotiffDataset);
    GDALClose(geotiffDataset_2);
    GDALDestroyDriverManager();
    return 0;
}