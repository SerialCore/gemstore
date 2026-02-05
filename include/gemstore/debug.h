/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_DEBUG
#define GEMSTORE_DEBUG

/* Debug su3 tensor product */
void debug_su3_product();

/* Debug spin-orbit coupling operators */
void debug_soc_operator();

/* Debug color wavefunctions */
void debug_color_wfn();

/* Debug spin wavefunctions */
void debug_spin_wfn();

/* Debug isospin wavefunctions */
void debug_isospin_wfn();

/* Debug orbital wavefunctions */
void debug_orbit_wfn();

/* Debug potential integral */
void debug_potential_model();

/* Debug eigen system */
void debug_eigen_system();

#endif