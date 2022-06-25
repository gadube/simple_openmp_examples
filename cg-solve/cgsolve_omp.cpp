// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER

/*
  Adapted from the Mantevo Miniapp Suite.
  https://mantevo.github.io/pdfs/MantevoOverview.pdf
*/

#include <cmath>
#include <limits>
#include <chrono>
#include <iostream>

#include "generate_matrix.hpp"

template <typename T>
void print_vec(T vec)
{
  for(size_t i=0; i<vec.extent(0); ++i) {
    std::cout << vec.name() << "["<<i<<"]: " << vec[i] << std::endl;
  }
}

template <class YType, class AType, class XType>
void spmv(YType y, AType A, XType x) {
        const size_t rows_size     = A.num_rows();
        const double* const xcoefs            = &x[0];
        double* ycoefs                        = &y[0];
        const double beta                     = 0;

        #pragma omp parallel for
        for(size_t row = 0; row < rows_size; ++row) {
                const size_t row_start = A.row_ptr[row];
                const size_t row_end   = A.row_ptr[row + 1];

                double sum = 0;

                #pragma omp simd reduction(+:sum)
                for(size_t i = row_start; i < row_end; ++i) {
                        sum += A.values(i) * x(A.col_idx(i));
                }

                ycoefs[row] = sum;
        }
}

template <class YType, class XType>
double dot(YType y, XType x) {
  double result = 0.0;
  const size_t n = y.extent(0);

  #pragma omp parallel for reduction(+:result)
  for(int i=0; i<n; ++i) {
    result += x[i] * y[i];
  }

  return result;
}

template <class ZType, class YType, class XType>
void axpby(ZType z, double alpha, XType x, double beta, YType y) {
    #ifdef MINIFE_DEBUG_OPENMP
      std::cout << "Starting WAXPBY..." << std::endl;
    #endif

      const INT_TYPE n = z.extent(0);
      const double*  xcoefs = &x[0];
      const double*  ycoefs = &y[0];
            double*  wcoefs = &z[0];

      if(beta == 0.0) {
    	if(alpha == 1.0) {
      		#pragma omp parallel for
      		for(int i=0; i<n; ++i) {
        			wcoefs[i] = xcoefs[i];
      		}
      	} else {
      		#pragma omp parallel for
      		for(int i=0; i<n; ++i) {
        			wcoefs[i] = alpha * xcoefs[i];
      		}
      	}
      } else {
    	if(alpha == 1.0) {
      		#pragma omp parallel for
      		for(int i=0; i<n; ++i) {
        			wcoefs[i] = xcoefs[i] + beta * ycoefs[i];
      		}
      	} else {
      		#pragma omp parallel for
      		for(int i=0; i<n; ++i) {
        			wcoefs[i] = alpha * xcoefs[i] + beta * ycoefs[i];
      		}
      	}
      }

    #ifdef MINIFE_DEBUG_OPENMP
      std::cout << "Finished WAXPBY." << std::endl;
    #endif
}

template <class VType, class AType>
int cg_solve(VType y, AType A, VType b, int max_iter, double tolerance) {
  int myproc    = 0;
  int num_iters = 0;

  double normr     = 0;
  double rtrans    = 0;
  double oldrtrans = 0;

  INT_TYPE print_freq = max_iter / 10;
  if (print_freq > 50) print_freq = 50;
  if (print_freq < 1) print_freq = 1;
  VType x("x", b.extent(0));
  VType r("r", x.extent(0));
  VType p("r", x.extent(0));
  VType Ap("r", x.extent(0));
  double one  = 1.0;
  double zero = 0.0;

  axpby(p, one, x, zero, x);
  spmv(Ap, A, p);
  axpby(r, one, b, -one, Ap);

  rtrans = dot(r, r);

  normr = std::sqrt(rtrans);

  if (myproc == 0) {
    std::cout << "Initial Residual = " << normr << std::endl;
  }

  double brkdown_tol = std::numeric_limits<double>::epsilon();

  for (INT_TYPE k = 1; k <= max_iter && normr > tolerance; ++k) {
    if (k == 1) {
      axpby(p, one, r, zero, r);
    } else {
      oldrtrans   = rtrans;
      rtrans      = dot(r, r);
      double beta = rtrans / oldrtrans;
      axpby(p, one, r, beta, p);
    }

    normr = std::sqrt(rtrans);

    if (myproc == 0 && (k % print_freq == 0 || k == max_iter)) {
      std::cout << "Iteration = " << k << "   Residual = " << normr
                << std::endl;
    }

    double alpha    = 0;
    double p_ap_dot = 0;

    spmv(Ap, A, p);

    p_ap_dot = dot(Ap, p);

    if (p_ap_dot < brkdown_tol) {
      if (p_ap_dot < 0) {
        std::cerr << "miniFE::cg_solve ERROR, numerical breakdown!"
                  << std::endl;
        return num_iters;
      } else
        brkdown_tol = 0.1 * p_ap_dot;
    }
    alpha = rtrans / p_ap_dot;

    axpby(x, one, x, alpha, p);
    axpby(r, one, r, -alpha, Ap);
    num_iters = k;
  }
  return num_iters;
}

int main(int argc, char* argv[]) {
  int N            = argc > 1 ? atoi(argv[1]) : 100;
  int max_iter     = argc > 2 ? atoi(argv[2]) : 200;
  double tolerance = argc > 3 ? atoi(argv[3]) : 0;

  CrsMatrix A = Impl::generate_miniFE_matrix(N);
  Kokkos::View<double> x = Impl::generate_miniFE_vector(N);

  Kokkos::View<INT_TYPE> row_ptr("row_ptr", A.row_ptr.extent(0));
  Kokkos::View<INT_TYPE> col_idx("col_idx", A.col_idx.extent(0));
  Kokkos::View<double> values("values", A.values.extent(0));

  Kokkos::View<double> y("Y", x.extent(0));

  std::cout << "============" << std::endl;
  std::cout << "WarmUp Solve" << std::endl;
  std::cout << "============" << std::endl << std::endl;
  int num_iters = cg_solve(y, A, x, 20, tolerance);

  std::cout << std::endl;
  std::cout << "============" << std::endl;
  std::cout << "Timing Solve" << std::endl;
  std::cout << "============" << std::endl << std::endl;
  auto start = std::chrono::high_resolution_clock::now();
  num_iters   = cg_solve(y, A, x, max_iter, tolerance);
  auto end = std::chrono::high_resolution_clock::now();
  double time = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();

  // Compute Bytes and Flops
  double spmv_bytes = A.num_rows() * sizeof(INT_TYPE) +
                      A.nnz() * sizeof(INT_TYPE) + A.nnz() * sizeof(double) +
                      A.nnz() * sizeof(double) + A.num_rows() * sizeof(double);

  double dot_bytes   = x.extent(0) * sizeof(double) * 2;
  double axpby_bytes = x.extent(0) * sizeof(double) * 3;

  double spmv_flops  = A.nnz() * 2;
  double dot_flops   = x.extent(0) * 2;
  double axpby_flops = x.extent(0) * 3;

  int spmv_calls  = 1 + num_iters;
  int dot_calls   = num_iters * 2;
  int axpby_calls = 2 + num_iters * 3;

  printf("CGSolve for 3D (%i %i %i); %i iterations; %lf time(ms)\n", N, N, N,
         num_iters, time);
  printf(
      "Performance: %lf GFlop/s %lf GB/s (Calls SPMV: %i Dot: %i AXPBY: %i\n",
      1e-9 *
          (spmv_flops * spmv_calls + dot_flops * dot_calls +
           axpby_flops * axpby_calls) /
          time,
      (1.0 / 1024 / 1024 / 1024) *
          (spmv_bytes * spmv_calls + dot_bytes * dot_calls +
           axpby_bytes * axpby_calls) /
          time,
      spmv_calls, dot_calls, axpby_calls);
}
