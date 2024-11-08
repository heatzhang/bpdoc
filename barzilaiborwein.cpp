//
// Created by zhangzhijie on 2024/10/15.
//

#include "barzilaiborwein.h"

namespace LayoutAdvance::Placer {

    void BBOptimizer::calculateNewSolutions() {
        size_t moduleCnt = curReferenceSolution.size();
        new_main_solution.resize(moduleCnt);
        new_reference_solution.resize(moduleCnt);

        // 计算步长
        double step_size;
        if (ps->iterCount == 0) {
            step_size = 1.0;
        } else {
            double lipschitz_constant = calcLipschitzConstant(
              curReferenceSolution, lastReferenceSolution, curGradient,
              lastGradient);
            double lip_step = 1.0 / lipschitz_constant;

            // 计算s和y向量
            double s_total = 0.0;
            double y_total = 0.0;
            double s_dot_y = 0.0;
            for (int i = 0; i < moduleCnt; i++) {
                VECTOR_2D s =
                  curReferenceSolution[i] - lastReferenceSolution[i];
                VECTOR_2D y = curGradient[i] - lastGradient[i];
                s_total += s.x * s.x + s.y * s.y;
                y_total += y.x * y.x + y.y * y.y;
                s_dot_y += s.x * y.x + s.y * y.y;
            }
            double s_norm_sq = s_total;
            double y_norm_sq = y_total;

            double bb_long_step = s_norm_sq / (s_dot_y > 0 ? s_dot_y : 1e-12);
            double bb_short_step = s_dot_y / (y_norm_sq > 0 ? y_norm_sq : 1e-12);

            if (bb_short_step > 0) {
                step_size = bb_short_step;
            } else {
                step_size = std::min(std::sqrt(s_norm_sq) / (std::sqrt(y_norm_sq) + 1e-12), lip_step);
            }
        }

        double new_NS_opt_param = (1 + std::sqrt(4 * (ps->nesterovPara * ps->nesterovPara) + 1)) / 2;

        // 执行优化步骤
        for (int i = 0; i < moduleCnt; i++) {
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
        ps->nesterovPara = new_NS_opt_param;
    }
} // LayoutAdvance::Placer