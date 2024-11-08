//
// Created by zhangzhijie on 2024/10/15.
//

#include "backtracking.h"

namespace LayoutAdvance::Placer {

    void BTOptimizer::calculateNewSolutions() {
        std::vector<VECTOR_2D> new_gradient;

        size_t moduleCnt = curReferenceSolution.size();
        new_main_solution.resize(moduleCnt);
        new_reference_solution.resize(moduleCnt);

        // 计算初始步长
        double step_size = (ps->iterCount == 0) ? 1.0 : 1.0 / calcLipschitzConstant(curReferenceSolution,
                                              lastReferenceSolution,
                                              curGradient, lastGradient);

        double new_NS_opt_param = (1 + std::sqrt(4 * (ps->nesterovPara * ps->nesterovPara) + 1)) / 2;

        while (true) {
            // 移动单元
            for (int i = 0; i < moduleCnt; ++i) {
                VECTOR_2D& new_pos = new_main_solution[i];
                VECTOR_2D& new_ref_pos = new_reference_solution[i];

                VECTOR_2D gradient = curGradient[i];
                VECTOR_2D cur_pos = curMainSolution[i];
                VECTOR_2D cur_ref_pos = curReferenceSolution[i];

                new_pos.x = cur_ref_pos.x + gradient.x * step_size;
                new_pos.y = cur_ref_pos.y + gradient.y * step_size;

                new_pos.x = std::min(std::max(new_pos.x, context->placeRegionLL.x + context->modules[i]->width * 0.5), context->placeRegionUR.x - context->modules[i]->width * 0.5);
                new_pos.y = std::min(std::max(new_pos.y, context->placeRegionLL.y + context->modules[i]->height * 0.5), context->placeRegionUR.y - context->modules[i]->height * 0.5);

                new_ref_pos.x = new_pos.x + (new_pos.x - cur_pos.x) * ((ps->nesterovPara - 1) / new_NS_opt_param);
                new_ref_pos.y = new_pos.y + (new_pos.y - cur_pos.y) * ((ps->nesterovPara - 1) / new_NS_opt_param);
            }

            context->setPositionFrom(new_reference_solution);
            model->totalGradientUpdate();
            context->writeGradientsTo(new_gradient);
            double new_lipschitz_constant = calcLipschitzConstant(
                new_reference_solution, curReferenceSolution, new_gradient,
                curGradient);
            double new_step_size = 1.0 / new_lipschitz_constant;

            if (BACKTRACK_EPS * step_size <= new_step_size) {
                break;
            }
            step_size = new_step_size;
        }
        ps->nesterovPara = new_NS_opt_param;
    }
} // LayoutAdvance::Placer