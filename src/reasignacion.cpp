//
// Created by oscar on 13/10/25.
//

#include "reasignacion.hpp"

void reasignar
(
    Campo::Presion  & presion,
    Campo::Momentum & velU,
    Campo::MDotStar & mdotstar,
    Campo::velFace  & velface
)
{

    // Campos nuevos
    const auto& u_star     = velU.u_star;
    const auto& v_star     = velU.v_star;
    const auto& P_star     = presion.P_star;
    const auto& mdotstar_x = mdotstar.mDotStar_x;
    const auto& mdotstar_y = mdotstar.mDotStar_y;
    const auto& ue         = velface.uFace_x;
    const auto& vn         = velface.vFace_y;
    const auto& ue_i       = velface.uFace_x_interp;
    const auto& vn_i       = velface.vFace_y_interp;

    // Campos viejos
    auto& u_old          = velU.u_old;
    auto& v_old          = velU.v_old;
    auto& P_old          = presion.P_old;
    auto& mdotstar_x_old = mdotstar.mDotStar_x_old;
    auto& mdotstar_y_old = mdotstar.mDotStar_y_old;
    auto& ue_n           = velface.uFace_x_n;
    auto& vn_n           = velface.vFace_y_n;
    auto& ue_i_n         = velface.uFace_x_interp_n;
    auto& vn_i_n         = velface.vFace_y_interp_n;

    // Reasignacion
    P_old          = P_star;
    u_old          = u_star;
    v_old          = v_star;
    mdotstar_x_old = mdotstar_x;
    mdotstar_y_old = mdotstar_y;
    ue_n           = ue;
    vn_n           = vn;
    ue_i_n         = ue_i;
    vn_i_n         = vn_i;

}
