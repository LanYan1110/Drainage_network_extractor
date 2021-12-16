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
    std::vector<int> pixels; // where everything is stored
    std::vector<int> is_in_lists;//initially 0, this cell is put in the search list then set to 1
    std::vector<int> is_processed;//initially 0, processed then set to 1
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

    void fillmarkers()
    {
        unsigned int total_pixels = max_x * max_y;
        for (int i = 0; i < total_pixels; ++i) is_in_lists.push_back(0);
        for (int i = 0; i < total_pixels; ++i) is_processed.push_back(0);
    }

    // Access the value of a raster cell to read or write it
    int& operator()(int x, int y) {
        assert(x >= 0 && x < max_x);
        assert(y >= 0 && y < max_y);
        return pixels[x + y * max_x];
    }

    // Access the value of a raster cell to read it
    int operator()(int x, int y) const {
        assert(x >= 0 && x < max_x);
        assert(y >= 0 && y < max_y);
        return pixels[x + y * max_x];
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
    bool operator<(const RasterCell& other) const {
        // to do with statements like if (this->elevation > other.elevation) return false/true;
        if (this->elevation < other.elevation)
            return false;
        else if (this->elevation > other.elevation)
            return true;
        else if (this->elevation == other.elevation)
        {
            //(grid cells added earlier have higher precedence in case of equal elevation)
            if (this->insertion_order < other.insertion_order)
                return false;
            else
                return true;
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
    if ((x >= 0 and x < nXSize and y >= 0 and y < nYSize and input_raster.is_in_lists[(x + y * input_raster.max_x)] == 0)
        and (flow_direction(x, y) == 0))
    {
        //add the neighbour into the prority queue
        RasterCell n(x, y, input_raster(x, y), insert_order++);
        cells_to_process_flow.push(n);
        input_raster.is_in_lists[(x + y * input_raster.max_x)] = 1;
        //write the cell's direction in the flow direction raster
        flow_direction.pixels[(x + y * input_raster.max_x)] = fd_value;
    }
    else if ((x >= 0 and x < nXSize and y >= 0 and y < nYSize and flow_direction(x, y) == 0 and input_raster.is_in_lists[x + y * input_raster.max_x] == -1)) {
        flow_direction.pixels[(x + y * input_raster.max_x)] = fd_value;
        input_raster.is_in_lists[(x + y * input_raster.max_x)] = 1;
    }

}

int flow_processing(Raster& input_raster, int x, int y, Raster& flow_direction, int fd_value, std::stack<RasterCell>& a, int nXSize, int nYSize)
{
    if (x >= 0 and x < nXSize and y >= 0 and y < nYSize)
    {
        if (flow_direction(x, y) == fd_value)  return 1;
        else return 0;
    }
}

int main(int argc, const char* argv[]) {

    // Open dataset
    GDALDataset* input_dataset;
    GDALAllRegister();
    input_dataset = (GDALDataset*)GDALOpen("E:\\Geomatics q2\\GEO1015\\runoffmodelling\\N36E076.hgt", GA_ReadOnly);
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

    input_raster.fillmarkers();
    // Flow direction
    Raster flow_direction(input_raster.max_x, input_raster.max_y);
    flow_direction.fill();
    std::priority_queue<RasterCell, std::deque<RasterCell>> cells_to_process_flow;
    std::stack<RasterCell> accumulation_stack;
    // Write flow direction
    //push boundary cells into priority queue (initial priority queue)
    for (int i = 0; i < input_raster.max_x; i++)
    {
        for (int j = 0; j < input_raster.max_y; j++)
        {
            if (i == 0 || j == 0 || i == input_raster.max_x - 1 || j == input_raster.max_y - 1)
            {
                RasterCell c(i, j, input_raster(i, j), insert_order++);
                cells_to_process_flow.push(c);
                input_raster.is_in_lists[i + j * input_raster.max_x] = -1;
            }
        };
    };
    //now we have our initial priority queue, we begin to process it
    for (; !cells_to_process_flow.empty(); cells_to_process_flow.pop())
    {
        RasterCell c = cells_to_process_flow.top();
        accumulation_stack.push(c);
        input_raster.is_processed[c.x * input_raster.max_x + c.y] = 1;
        //get its neighbours
        //if they are four corners
        neighbour_processing(input_raster, c.x - 1, c.y - 1, flow_direction, 25, cells_to_process_flow, nXSize, nYSize); //1
        neighbour_processing(input_raster, c.x, c.y - 1, flow_direction, 20, cells_to_process_flow, nXSize, nYSize);//2
        neighbour_processing(input_raster, c.x + 1, c.y - 1, flow_direction, 35, cells_to_process_flow, nXSize, nYSize);//3
        neighbour_processing(input_raster, c.x - 1, c.y, flow_direction, 30, cells_to_process_flow, nXSize, nYSize); //4
        neighbour_processing(input_raster, c.x + 1, c.y, flow_direction, 40, cells_to_process_flow, nXSize, nYSize); //5
        neighbour_processing(input_raster, c.x - 1, c.y + 1, flow_direction, 15, cells_to_process_flow, nXSize, nYSize); //6
        neighbour_processing(input_raster, c.x, c.y + 1, flow_direction, 10, cells_to_process_flow, nXSize, nYSize); //7
        neighbour_processing(input_raster, c.x + 1, c.y + 1, flow_direction, 5, cells_to_process_flow, nXSize, nYSize); //8+
    }

    //flow accumulation
    Raster flow_accumulation(input_raster.max_x, input_raster.max_y);
    flow_accumulation.fill_1();
    for (; !accumulation_stack.empty(); accumulation_stack.pop())
    {
        RasterCell c = accumulation_stack.top();
        //get its neighbours
        if (c.x - 1 >= 0 and c.x - 1 < nXSize and c.y - 1 >= 0 and c.y - 1 < nYSize and flow_direction(c.x, c.y - 1) == 25)
            flow_accumulation.pixels[(c.x + c.y * input_raster.max_x)] += flow_accumulation(c.x - 1, c.y - 1);
        if (c.x >= 0 and c.x < nXSize and c.y-1 >= 0 and c.y-1 < nYSize and flow_direction(c.x, c.y-1) == 20)
            flow_accumulation(c.x, c.y) += flow_accumulation(c.x, c.y-1);
        if (c.x + 1 >= 0 and c.x + 1 < nXSize and c.y - 1 >= 0 and c.y - 1 < nYSize and flow_direction(c.x + 1, c.y - 1) == 35)
            flow_accumulation(c.x, c.y) += flow_accumulation(c.x + 1, c.y - 1);
        if (c.x - 1 >= 0 and c.x - 1 < nXSize and c.y >= 0 and c.y < nYSize and flow_direction(c.x - 1, c.y) == 30)
            flow_accumulation(c.x, c.y) += flow_accumulation(c.x - 1, c.y);
        if (c.x + 1 >= 0 and c.x + 1 < nXSize and c.y >= 0 and c.y < nYSize and flow_direction(c.x + 1, c.y) == 40)
            flow_accumulation(c.x, c.y) += flow_accumulation(c.x + 1, c.y);
        if (c.x - 1 >= 0 and c.x - 1 < nXSize and c.y + 1 >= 0 and c.y + 1 < nYSize and flow_direction(c.x - 1, c.y + 1) == 15)
            flow_accumulation(c.x, c.y) += flow_accumulation(c.x - 1, c.y + 1);
        if (c.x >= 0 and c.x < nXSize and c.y + 1 >= 0 and c.y + 1 < nYSize and flow_direction(c.x, c.y + 1) == 10)
            flow_accumulation(c.x, c.y) += flow_accumulation(c.x, c.y + 1);
        if (c.x + 1 >= 0 and c.x + 1 < nXSize and c.y +1 >= 0 and c.y +1 < nYSize and flow_direction(c.x + 1, c.y + 1) == 5)
            flow_accumulation(c.x, c.y) += flow_accumulation(c.x + 1, c.y + 1);
    }


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
    geotiffDataset = driverGeotiff->Create(tiffname.c_str(), nXSize, nYSize, 1, GDT_Float32, NULL);
    geotiffDataset_2 = driverGeotiff_2->Create(tiffname_2.c_str(), nXSize, nYSize, 1, GDT_Float32, NULL);
    geotiffDataset->SetGeoTransform(geo_transform);
    geotiffDataset_2->SetGeoTransform(geo_transform);
    geotiffDataset->SetProjection(input_dataset->GetProjectionRef());
    geotiffDataset_2->SetProjection(input_dataset->GetProjectionRef());
    int* rowBuff = (int*)CPLMalloc(sizeof(float) * nXSize);
    int* rowBuff_2 = (int*)CPLMalloc(sizeof(float) * nXSize);
    for (int row = 0; row < nYSize; row++)
    {
        for (int col = 0; col < nXSize; col++)
        {
            rowBuff[col] = flow_accumulation(col, row);
        }
        geotiffDataset->GetRasterBand(1)->RasterIO(GF_Write, 0, row, nXSize, 1, rowBuff, nXSize, 1, GDT_Int32, 0, 0);
    }
    for (int row = 0; row < nYSize; row++)
    {
        for (int col = 0; col < nXSize; col++)
        {
            rowBuff_2[col] = flow_direction(col, row);
        }
        geotiffDataset_2->GetRasterBand(1)->RasterIO(GF_Write, 0, row, nXSize, 1, rowBuff_2, nXSize, 1, GDT_Int32, 0, 0);
    }
    GDALClose(input_dataset);
    GDALClose(geotiffDataset);
    GDALClose(geotiffDataset_2);
    GDALDestroyDriverManager();
    return 0;
}