/*!
 * 
 *
 * \file
 *
 *
 * \par Copyright 1995-2017 Shark Development Team
 * 
 * <BR><HR>
 * This file is part of Shark.
 * <http://shark-ml.org/>
 * 
 * Shark is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published 
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Shark is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with Shark.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
//#include <shark/ObjectiveFunctions/Benchmarks/cec09.h>
#include <shark/ObjectiveFunctions/Benchmarks/ZDT6.h>
#include <shark/ObjectiveFunctions/Benchmarks/ZDT4.h>
#include <shark/ObjectiveFunctions/Benchmarks/ZDT3.h>
#include <shark/ObjectiveFunctions/Benchmarks/ZDT2.h>
#include <shark/ObjectiveFunctions/Benchmarks/ZDT1.h>
#include <shark/ObjectiveFunctions/Benchmarks/Himmelblau.h>
#include <shark/ObjectiveFunctions/Benchmarks/Sphere.h>
#include <shark/ObjectiveFunctions/Benchmarks/Rosenbrock.h>
#include <shark/ObjectiveFunctions/Benchmarks/LZ9.h>
#include <shark/ObjectiveFunctions/Benchmarks/LZ8.h>
#include <shark/ObjectiveFunctions/Benchmarks/LZ7.h>
#include <shark/ObjectiveFunctions/Benchmarks/LZ6.h>
#include <shark/ObjectiveFunctions/Benchmarks/LZ5.h>
#include <shark/ObjectiveFunctions/Benchmarks/LZ4.h>
#include <shark/ObjectiveFunctions/Benchmarks/LZ3.h>
#include <shark/ObjectiveFunctions/Benchmarks/LZ2.h>
#include <shark/ObjectiveFunctions/Benchmarks/LZ1.h>
#include <shark/ObjectiveFunctions/Benchmarks/IHR6.h>
#include <shark/ObjectiveFunctions/Benchmarks/IHR4.h>
#include <shark/ObjectiveFunctions/Benchmarks/IHR3.h>
#include <shark/ObjectiveFunctions/Benchmarks/IHR2.h>
#include <shark/ObjectiveFunctions/Benchmarks/IHR1.h>
#include <shark/ObjectiveFunctions/Benchmarks/GSP.h>
#include <shark/ObjectiveFunctions/Benchmarks/Fonseca.h>
#include <shark/ObjectiveFunctions/Benchmarks/Ellipsoid.h>
#include <shark/ObjectiveFunctions/Benchmarks/ELLI2.h>
#include <shark/ObjectiveFunctions/Benchmarks/ELLI1.h>
#include <shark/ObjectiveFunctions/Benchmarks/DTLZ7.h>
#include <shark/ObjectiveFunctions/Benchmarks/DTLZ6.h>
#include <shark/ObjectiveFunctions/Benchmarks/DTLZ5.h>
#include <shark/ObjectiveFunctions/Benchmarks/DTLZ4.h>
#include <shark/ObjectiveFunctions/Benchmarks/DTLZ3.h>
#include <shark/ObjectiveFunctions/Benchmarks/DTLZ2.h>
#include <shark/ObjectiveFunctions/Benchmarks/DTLZ1.h>
#include <shark/ObjectiveFunctions/Benchmarks/Discus.h>
#include <shark/ObjectiveFunctions/Benchmarks/DiffPowers.h>
#include <shark/ObjectiveFunctions/Benchmarks/CIGTAB2.h>
#include <shark/ObjectiveFunctions/Benchmarks/CIGTAB1.h>
#include <shark/ObjectiveFunctions/Benchmarks/CigarDiscus.h>
#include <shark/ObjectiveFunctions/Benchmarks/Cigar.h>
#include <shark/ObjectiveFunctions/Benchmarks/Ackley.h>
#include <shark/ObjectiveFunctions/Benchmarks/UF1.h>
#include <shark/ObjectiveFunctions/Benchmarks/UF2.h>
#include <shark/ObjectiveFunctions/Benchmarks/UF3.h>
#include <shark/ObjectiveFunctions/Benchmarks/UF4.h>
#include <shark/ObjectiveFunctions/Benchmarks/UF5.h>
#include <shark/ObjectiveFunctions/Benchmarks/UF6.h>
#include <shark/ObjectiveFunctions/Benchmarks/UF7.h>
#include <shark/ObjectiveFunctions/Benchmarks/UF8.h>
#include <shark/ObjectiveFunctions/Benchmarks/UF9.h>
#include <shark/ObjectiveFunctions/Benchmarks/UF10.h>
#include <shark/ObjectiveFunctions/Benchmarks/WFG1.h>
#include <shark/ObjectiveFunctions/Benchmarks/WFG2.h>
#include <shark/ObjectiveFunctions/Benchmarks/WFG3.h>
#include <shark/ObjectiveFunctions/Benchmarks/WFG4.h>
#include <shark/ObjectiveFunctions/Benchmarks/WFG5.h>
#include <shark/ObjectiveFunctions/Benchmarks/WFG6.h>
#include <shark/ObjectiveFunctions/Benchmarks/WFG7.h>
#include <shark/ObjectiveFunctions/Benchmarks/WFG8.h>
#include <shark/ObjectiveFunctions/Benchmarks/WFG9.h>
#include <shark/ObjectiveFunctions/Benchmarks/BT1.h>
#include <shark/ObjectiveFunctions/Benchmarks/BT2.h>
#include <shark/ObjectiveFunctions/Benchmarks/BT3.h>
#include <shark/ObjectiveFunctions/Benchmarks/BT4.h>
#include <shark/ObjectiveFunctions/Benchmarks/BT5.h>
#include <shark/ObjectiveFunctions/Benchmarks/BT6.h>
#include <shark/ObjectiveFunctions/Benchmarks/BT7.h>
#include <shark/ObjectiveFunctions/Benchmarks/BT8.h>
#include <shark/ObjectiveFunctions/Benchmarks/BT9.h>
#include <shark/ObjectiveFunctions/Benchmarks/IMB1.h>
#include <shark/ObjectiveFunctions/Benchmarks/IMB2.h>
#include <shark/ObjectiveFunctions/Benchmarks/IMB3.h>
#include <shark/ObjectiveFunctions/Benchmarks/IMB4.h>
#include <shark/ObjectiveFunctions/Benchmarks/IMB5.h>
#include <shark/ObjectiveFunctions/Benchmarks/IMB6.h>
#include <shark/ObjectiveFunctions/Benchmarks/IMB7.h>
#include <shark/ObjectiveFunctions/Benchmarks/IMB8.h>
#include <shark/ObjectiveFunctions/Benchmarks/IMB9.h>
#include <shark/ObjectiveFunctions/Benchmarks/IMB10.h>

/// \defgroup benchmarks Benchmark functions
/// \ingroup objfunctions
///
/// A large set of benchmarks for single and multi-objective optimization.
/// Benchmarks are in the namespace shark::benchmark

