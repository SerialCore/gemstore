/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_PARAM_TYPECC
#define GEMSTORE_PARAM_TYPECC

#include <iostream>
#include <fstream>

class DualStream {
public:
    DualStream(const std::string& filename = "Fitting.out")
        : console(std::cout), file(filename, std::ios::out | std::ios::app) {
        if (!file.is_open()) {
            std::cerr << "Cannot open log file: " << filename << std::endl;
        }
    }

    DualStream(const DualStream&) = delete;
    DualStream& operator=(const DualStream&) = delete;
    DualStream(DualStream&&) = delete;
    DualStream& operator=(DualStream&&) = delete;

    ~DualStream() {
        if (file.is_open()) {
            file << std::endl;
            file.close();
        }
    }

    template<typename T>
    DualStream& operator<<(const T& value) {
        console << value;
        if (file.is_open()) {
            file << value;
        }
        return *this;
    }

    DualStream& operator<<(std::ostream& (*manip)(std::ostream&)) {
        console << manip;
        if (file.is_open()) {
            file << manip;
        }
        return *this;
    }
    
private:
    std::ostream& console;
    std::ofstream file;
};

struct State {
    int f1, f2, N, S, L, J;		/* quantum numbers */
	double exp_mass;			/* experimental mass */
	double exp_error;			/* experimental error */
};

#endif