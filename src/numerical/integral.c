/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/numerical/integral.h>
#include <gemstore/model/model.h>
#include <gemstore/basis/orbit.h>

#define OHP 50          /* order of half-range Hermite polynomials */

/* weights of half-range Hermite polynomials */
static const double weights[OHP] = {
    9.49074370469973147082E-003, 2.20248042432282553631E-002, 3.43744122779471227909E-002, 4.62954423315795433979E-002, 5.74270673524102602228E-002,
    6.72616743347232748312E-002, 7.51655105427606549250E-002, 8.04487210189303602989E-002, 8.24866653780293849459E-002, 8.08772841479732588309E-002,
    7.55968131214373946840E-002, 6.70988969096135482036E-002, 5.63037514101928781474E-002, 4.44533647270852753206E-002, 3.28602655529441186895E-002,
    2.26281047022441041912E-002, 1.44421205638529694046E-002, 8.49988917441698063564E-003, 4.58981582947691495106E-003, 2.26248032155746300357E-003,
    1.01296345857313001359E-003, 4.09855051025779875808E-004, 1.49103396554347114376E-004, 4.85209725035066145407E-005, 1.40499187954310437190E-005,
    3.60057445584129230942E-006, 8.12054885660541412344E-007, 1.60236226323029379083E-007, 2.74914498210973363181E-008, 4.07393870688604362262E-009,
    5.17733909405426558330E-010, 5.59883137229071183484E-011, 5.10833440280528370032E-012, 3.89531180860791485073E-013, 2.45631708315702786126E-014,
    1.26561330931793734941E-015, 5.25582653333676002424E-017, 1.73149373672351779729E-018, 4.44198240366431239009E-020, 8.68013900359727462358E-022,
    1.25805027635501783452E-023, 1.30872128650771015724E-025, 9.37652553566303622113E-028, 4.38604154939073732093E-030, 1.24690572705593213602E-032,
    1.94847016939779663560E-035, 1.44052991578865261374E-038, 3.94314802874471338715E-042, 2.51219924772053768487E-046, 1.15324953442160894479E-051
};

/* nodes of half-range Hermite polynomials */
static const double nodes[OHP] = {
    3.69941668941189387078E-003, 1.94655785593001262377E-002, 4.77232575110775569884E-002, 8.82994741423326136810E-002, 1.40944249831512837831E-001,
    2.05343824732975154536E-001, 2.81130944731155938791E-001, 3.67895701867100816092E-001, 4.65196632382567480466E-001, 5.72571575839739193266E-001,
    6.89547864293400636314E-001, 8.15651524129035456870E-001, 9.50415293767051842487E-001, 1.09338537085220496150E+000, 1.24412689394047714663E+000,
    1.40222823252836353434E+000, 1.56730420579590844586E+000, 1.73899837732924138424E+000, 1.91698458425367595630E+000, 2.10096785886658098596E+000,
    2.29068489295538163123E+000, 2.48590418287233886279E+000, 2.68642597976929227008E+000, 2.89208215616667258410E+000, 3.10273608869143607692E+000,
    3.31828264842635672287E+000, 3.53864838570112628139E+000, 3.76379199610129929525E+000, 3.99370515987453091261E+000, 4.22841385900133938902E+000,
    4.46798029678176976527E+000, 4.71250557663335036490E+000, 4.96213334417320913961E+000, 5.21705466625021294844E+000, 5.47751452299329293275E+000,
    5.74382044124152601111E+000, 6.01635402812180604868E+000, 6.29558651982983683470E+000, 6.58210002638651010927E+000, 6.87661707965657488199E+000,
    7.18004266516692291795E+000, 7.49352570514240963406E+000, 7.81855215039574067889E+000, 8.15709210450168878089E+000, 8.51184526474371849140E+000,
    8.88668005924412897894E+000, 9.28749674141648604654E+000, 9.72416586588463146083E+000, 1.02158862585784281522E+001, 1.08129860729453608573E+001
};

double integral_wfn_overlap(
    orbit_wfn_t wfn,
    double node_factor,
    const argsOrbit_t *args_bra,
    const argsOrbit_t *args_ket)
{
    double sum = 0.0;

    for (int i = 0; i < OHP; i++) {
        sum += node_factor * weights[i]
             * wfn(node_factor * nodes[i], args_bra->n, args_bra->l, args_bra->scale)
             * wfn(node_factor * nodes[i], args_ket->n, args_ket->l, args_ket->scale)
             * node_factor * node_factor * nodes[i] * nodes[i];
    }

    return sum;
}

double integral_wfn_overlap_complex(
    orbit_wfn_complex_t wfn,
    double node_factor,
    const argsOrbit_t *args_bra,
    const argsOrbit_t *args_ket)
{
    double sum = 0.0;

    for (int i = 0; i < OHP; i++) {
        sum += node_factor * weights[i]
             * wfn(node_factor * nodes[i], args_bra->n, args_bra->l, args_bra->scale)
             * wfn(node_factor * nodes[i], args_ket->n, args_ket->l, args_ket->scale)
             * node_factor * node_factor * nodes[i] * nodes[i];
    }

    return sum;
}

double integral_matrix_element(
    orbit_wfn_t wfn,
    potential_t pot,
    double node_factor,
    const argsOrbit_t *args_bra,
    const argsOrbit_t *args_ket,
    const argsModel_t *args_model,
    const argsModelDy_t *args_dynmc)
{
    double sum = 0.0;

    for (int i = 0; i < OHP; i++) {
        sum += node_factor * weights[i]
             * wfn(node_factor * nodes[i], args_bra->n, args_bra->l, args_bra->scale)
             * wfn(node_factor * nodes[i], args_ket->n, args_ket->l, args_ket->scale)
             * pot(node_factor * nodes[i], args_model, args_dynmc)
             * node_factor * node_factor * nodes[i] * nodes[i];
    }

    return sum;
}

double integral_matrix_element_complex(
    orbit_wfn_complex_t wfn,
    potential_t pot,
    double node_factor,
    const argsOrbit_t *args_bra,
    const argsOrbit_t *args_ket,
    const argsModel_t *args_model,
    const argsModelDy_t *args_dynmc)
{
    double sum = 0.0;

    for (int i = 0; i < OHP; i++) {
        sum += node_factor * weights[i]
             * wfn(node_factor * nodes[i], args_bra->n, args_bra->l, args_bra->scale)
             * wfn(node_factor * nodes[i], args_ket->n, args_ket->l, args_ket->scale)
             * pot(node_factor * nodes[i], args_model, args_dynmc)
             * node_factor * node_factor * nodes[i] * nodes[i];
    }

    return sum;
}