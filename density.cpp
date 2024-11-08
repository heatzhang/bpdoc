//
// Created by zhangzhijie on 2024/10/14.
//

#include "density.h"

namespace LayoutAdvance::Placer {
    DensityCalculator::DensityCalculator(Context *ctx)
            : context_(ctx) {
        initializeBins();
//        expandSmallModulesToBinSize();
        initializeFFT();
    }

    void DensityCalculator::densityGradientUpdate() {
        // in density gradient update, filler cells are treated like normal modules

        updateBinDensity();

        executeFFT();

        propagateElectroForceToCells();
    }

    void DensityCalculator::updateBinDensity() {
        // first we need to clear the node and filler density of all bins, they are outdated
        for (int i = 0; i < context_->binDimensionX; i++) {
            for (int j = 0; j < context_->binDimensionY; j++) {
                context_->bins[i][j]->nodeDensity = 0.0;
                context_->bins[i][j]->fillerDensity = 0.0;
            }
        }

        // second we iterate over all modules and calculate their density contribution to each bin
        for (Cell *curModule: context_->modules) {
            // we need to calculate the touched bins of the module
            auto trueLL = context_->binStep;
            trueLL.x *= context_->ringThicknessX;
            trueLL.y *= context_->ringThicknessY;
            trueLL = context_->placeRegionLL - trueLL;
            curModule->updateBinIdx(context_->binDimensionX,
                                    context_->binDimensionY,
                                    context_->binStep,
                                    trueLL);

            for (int i = curModule->binStartIdx_x; i <= curModule->binEndIdx_x; i++) {
                for (int j = curModule->binStartIdx_y; j <= curModule->binEndIdx_y; j++) {
                    const auto &bin = context_->bins[i][j];
                    double overlapArea = curModule->getOverlapArea(bin.get());
                    overlapArea *= curModule->localSmoothLengthScale.x *
                                   curModule->localSmoothLengthScale.y;
                    if (curModule->isFiller) {
                        context_->bins[i][j]->fillerDensity += overlapArea;
                    } else {
                        context_->bins[i][j]->nodeDensity += overlapArea;
                    }
                }
            }
        }
    }

    double DensityCalculator::calculateUnionAreaOfModules() {
        // // 使用布局区域的边界
        // double xmin = context_->placeRegionLL.x, ymin = context_->placeRegionLL.y;
        // double xmax = context_->placeRegionUR.x, ymax = context_->placeRegionUR.y;
        // for (const Cell* module : context_->modules) {
        //     xmin = std::min(xmin, module->coordinate.x);
        //     ymin = std::min(ymin, module->coordinate.y);
        //     xmax = std::max(xmax, module->coordinate.x + module->getWidth());
        //     ymax = std::max(ymax, module->coordinate.y + module->getHeight());
        // }

        // // 使用bin的大小的10倍作为网格大小，提高精度
        // int ratioX = 200;
        // int ratioY = 3;
        // int grid_size_x = context_->binDimensionX * ratioX;
        // int grid_size_y = context_->binDimensionY * ratioY;
        // double dx = (xmax - xmin) / grid_size_x;
        // double dy = (ymax - ymin) / grid_size_y;
        

        // // 初始化网格
        // std::vector<std::vector<bool>> grid(grid_size_x, std::vector<bool>(grid_size_y, false));

        // // 标记被模块覆盖的网格
        // for (const Cell* module : context_->modules) {
        //     if (module->isFiller) continue; // 跳过填充单元

        //     int x_start = std::max(0, int((module->coordinate.x - xmin) / dx));
        //     int x_end = std::min(grid_size_x - 1, int((module->coordinate.x + module->getWidth() - xmin) / dx));
        //     int y_start = std::max(0, int((module->coordinate.y - ymin) / dy));
        //     int y_end = std::min(grid_size_y - 1, int((module->coordinate.y + module->getHeight() - ymin) / dy));

        //     for (int x = x_start; x <= x_end; ++x) {
        //         for (int y = y_start; y <= y_end; ++y) {
        //             grid[x][y] = true;
        //         }
        //     }
        // }

        // // 计算被覆盖的网格总数
        // int count = 0;
        // std::vector<std::vector<int>> bin_count(context_->binDimensionX, std::vector<int>(context_->binDimensionY, 0));
        // for (int x = 0; x < grid_size_x; ++x) {
        //     for (int y = 0; y < grid_size_y; ++y) {
        //         if (grid[x][y]) {
        //             count++;
        //             bin_count[x/ratioX][y/ratioY]++;
        //         }
        //     }
        // }

        
        // auto binStepX = context_->binStep.x;
        // auto binStepY = context_->binStep.y;
        // double globalOverflowArea = 0.0;
        // double binArea = binStepX * binStepY;
        // double scaledBinArea = context_->ps->targetDensity * binArea;
        // double invertedBinArea = 1.0 / binArea;
        // double placementRegionArea = (context_->placeRegionUR.x - context_->placeRegionLL.x) * (context_->placeRegionUR.y - context_->placeRegionLL.y);

        // double estimatedStdCellsAreaRatio = context_->fillerCellsTotalArea / std::max(1e-6, placementRegionArea - count * dx * dy);
        // // estimatedStdCellsAreaRatio = std::min(estimatedStdCellsAreaRatio, 1.0);
        // Logger::info("Estimated std cells area ratio =", estimatedStdCellsAreaRatio);
        // for (int i = 0; i < context_->binDimensionX; i++) {
        //     for (int j = 0; j < context_->binDimensionY; j++) {
        //         auto& bin = context_->bins[i][j];
        //         double occupiedArea = double(bin_count[i][j]) / double(ratioX * ratioY) * binArea;
        //         bin->fillerDensity = std::max(0.0, (binArea - occupiedArea) * estimatedStdCellsAreaRatio);
        //     }
        // }
        // // 返回近似并集面积
        // return count * dx * dy;
        return 0;
    }

    void DensityCalculator::executeFFT() {
        if (!context_->ps->addFillers && context_->ps->smartFiller) {
            calculateUnionAreaOfModules();
        }
        auto binStepX = context_->binStep.x;
        auto binStepY = context_->binStep.y;
        double globalOverflowArea = 0.0;
        double binArea = binStepX * binStepY;
        double scaledBinArea = context_->ps->targetDensity * binArea;
        double invertedBinArea = 1.0 / binArea;
        double placementRegionArea = (context_->placeRegionUR.x - context_->placeRegionLL.x) * (context_->placeRegionUR.y - context_->placeRegionLL.y);
        for (int i = 0; i < context_->bins.size(); i++) {
            for (int j = 0; j < context_->bins[i].size(); j++) {
                double eDensity =
                        context_->bins[i][j]->nodeDensity +
                        context_->bins[i][j]->fillerDensity +
                        context_->bins[i][j]->terminalDensity;
                if (context_->bins[i][j]->inPlaceRegion) {// note: overflow calculation do not consider filler cells
                    globalOverflowArea += std::max(0.0, eDensity -
                                                   scaledBinArea);
                }
                eDensity += context_->bins[i][j]->baseDensity;
                eDensity *= invertedBinArea;
                fft->updateDensity(i, j, eDensity);
            }
        }

        context_->ps->globalDensityOverflow = globalOverflowArea / context_->stdCellsTotalArea;
        fft->doFFT();
        for (int i = 0; i < context_->binDimensionX; i++) {
            for (int j = 0; j < context_->binDimensionY; j++) {
                auto eForcePair = fft->getElectroForce(i, j);
                context_->bins[i][j]->E.x = eForcePair.first;
                context_->bins[i][j]->E.y = eForcePair.second;
            }
        }
    }

    void DensityCalculator::propagateElectroForceToCells() {
        int index = 0;
        for (Cell *curModule: context_->modules) {
            assert(index == curModule->idx);
            double ratio =
                    curModule->localSmoothLengthScale.x * curModule->localSmoothLengthScale.y;
            for (int i = curModule->binStartIdx_x; i <= curModule->binEndIdx_x; i++) {
                for (int j = curModule->binStartIdx_y; j <= curModule->binEndIdx_y; j++) {
                    double overlapArea = curModule->getOverlapArea(context_->bins[i][j].get());
                    overlapArea *= ratio;
                    context_->densityGradientX[index] += overlapArea * context_->bins[i][j]->E.x;
                    context_->densityGradientY[index] += overlapArea * context_->bins[i][j]->E.y;
                }
            }
            double leftArea = std::max(0.0, context_->placeRegionLL.x - curModule->coordinate.x) *
                              curModule->getHeight();
            double rightArea = std::max(0.0, curModule->coordinate.x + curModule->getWidth() -
                    context_->placeRegionUR.x) * curModule->getHeight();
            double lowerArea = std::max(0.0, context_->placeRegionLL.y - curModule->coordinate.y) *
                               curModule->getWidth();
            double upperArea = std::max(0.0, curModule->coordinate.y + curModule->getHeight() -
                    context_->placeRegionUR.y) * curModule->getWidth();
            leftArea *= ratio;
            rightArea *= ratio;
            lowerArea *= ratio;
            upperArea *= ratio;
            if (curModule->isFiller) {
//                context_->densityGradientX[index] += (leftArea - rightArea) * 0.01;
//                context_->densityGradientY[index] += (lowerArea - upperArea) * 0.01;
            } else {
//                context_->densityGradientX[index] += (leftArea - rightArea) * 0.1;
//                context_->densityGradientY[index] += (lowerArea - upperArea) * 0.1;
            }
            index++;
        }
    }

    void DensityCalculator::initializeFFT() {
        fft = std::make_unique<FFTW3Wrapper>(context_->binDimensionX, context_->binDimensionY, context_->binStep.x,
                                                context_->binStep.y);
        fft->createPlans(context_->ps->useWisdom);
    }

    void DensityCalculator::expandSmallModulesToBinSize() {
        for (Cell *curModule: context_->modules) {
            // some modules may be too small, we need to expand their size to at least the bin size
            curModule->reshape(context_->binStep);
        }
    }

    void DensityCalculator::initializeBins() {
        auto& modules = context_->modules;
        auto& placeRegionUR = context_->placeRegionUR;
        auto& placeRegionLL = context_->placeRegionLL;
        auto& binDimensionX = context_->binDimensionX;
        auto& binDimensionY = context_->binDimensionY;
        auto& binStep = context_->binStep;
        auto& bins = context_->bins;

        auto smallest_module_x = std::min_element(modules.begin(), modules.end(), [](Cell* a, Cell* b) {
            return a->getWidth() < b->getWidth();
        });
        auto smallest_module_y = std::min_element(modules.begin(), modules.end(), [](Cell* a, Cell* b) {
            return a->getHeight() < b->getHeight();
        });

        if (smallest_module_x != modules.end() && smallest_module_y != modules.end()) {
            Cell* minModuleX = *smallest_module_x;
            Cell* minModuleY = *smallest_module_y;
            int idealBinCntX = std::ceil(BIN_RATIO*(placeRegionUR.x - placeRegionLL.x) / minModuleX->getWidth());
            int idealBinCntY = std::ceil(BIN_RATIO*(placeRegionUR.y - placeRegionLL.y) / minModuleY->getHeight());
            double idealBinWidth = (placeRegionUR.x - placeRegionLL.x) / idealBinCntX;
            double idealBinHeight = (placeRegionUR.y - placeRegionLL.y) / idealBinCntY;
            double bestBinEdgeLength = std::min(idealBinWidth, idealBinHeight);
            binDimensionX = std::ceil((placeRegionUR.x - placeRegionLL.x) / bestBinEdgeLength);
            binDimensionY = std::ceil((placeRegionUR.y - placeRegionLL.y) / bestBinEdgeLength);
            int cntX = 1;
            int cntY = 1;
            while ((2 << cntX) < binDimensionX && cntX <= 10) {
                cntX++;
            }
            while ((2 << cntY) < binDimensionY && cntY <= 10) {
                cntY++;
            }
            binDimensionX = binDimensionY = 2 << std::max(cntY, cntX);
        } else {
            std::cerr << "No modules found." << std::endl;
            binDimensionX = binDimensionY = 1024;
        }

        binStep.x = (placeRegionUR.x - placeRegionLL.x) / binDimensionX;
        binStep.y = (placeRegionUR.y - placeRegionLL.y) / binDimensionY;
        Logger::info("Bin dimension is", binDimensionX, "*", binDimensionY, "size =", binStep);

        bins.resize(binDimensionX);
        for (int i = 0; i < binDimensionX; i++) {
            bins[i].resize(binDimensionY);
            for (int j = 0; j < binDimensionY; j++) {
                bins[i][j] = std::make_unique<Bin>();
                auto& curBin = bins[i][j];
                curBin->ll.x = i * binStep.x + placeRegionLL.x;
                curBin->ll.y = j * binStep.y + placeRegionLL.y;
                curBin->width = binStep.x;
                curBin->height = binStep.y;
                curBin->ur.x = curBin->ll.x + curBin->width;
                curBin->ur.y = curBin->ll.y + curBin->height;
            }
        }
    }
} // LayoutAdvance::Placer